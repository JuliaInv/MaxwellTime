#Test sensitivities for backward Euler forward code with OcTree Mesh and
#grounded source. Single source and multiple source configurations tested.

using MaxwellTime
using JOcTree
using Base.Test
using jInv.LinearSolvers
using jInv.Utils

@testset "Grd src BE time-stepping" begin
println("Testing single source")

L  = [4096, 4096, 2048.0]
x0 = [0.0,  0.0,  0.0]
n  = [16, 16, 16]
h  = L./n

xmin = 1548.0; xmax = 2548.0
ymin = 1548.0; ymax = 2548.0
zmin =  524.0; zmax = 1524.0


S    = createOcTreeFromBox(x0[1],x0[2],x0[3],n[1],n[2],n[3],h[1],h[2],h[3],xmin,xmax,ymin,ymax,zmin,zmax,1,1)
Msh  = getOcTreeMeshFV(S,h;x0=x0)
ne   = sum(Msh.ne)
Xn   = getNodalGrid(Msh)
nn   = size(Xn,1)

#First test that sensitivities from getSensMatVec and getSensTMatVec
#Match those computed from explicitly constructing J.
Sources = zeros(ne)
Tx = [1792.0 2048.0 1024.0;
      2048.0 2048.0 1024.0;]

Sources = getEdgeIntegralOfPolygonalChain(Msh,Tx)

#Setup receivers
Iobs = round.(Integer,ceil.(Msh.ne[1]*rand(16)))
P    = speye(ne)
P    = P[Iobs,:]

#Time stepping
dt      = [1e-4;1.5e-4;3e-4;4e-4]
t       = [0.0;cumsum(dt)]
wave    = zeros(length(dt)+1)
wave[1] = 1.0

#Random model
a     = rand(Msh.nc)
b     = 10*rand(Msh.nc) - 4
sigma = a.^b
mu0   = 4*pi*1e-7  # magnetic permeability
#chi   = 100*rand(Msh.nc)
chi   = zeros(Msh.nc)
mu    = mu0*(1+chi)

#Get data at initial model
sourceType = :Galvanic
obsTimes   = [0;cumsum(dt)]
t0         = 0.0
pFor       = getMaxwellTimeParam(Msh,Sources,P',obsTimes,t0,dt,wave,sourceType)
println("Getting data")
d,pFor = getData(sigma,pFor)
ew     = pFor.fields

#Form sensitivity matrix and related quantities
println("Forming sensitivity matrix")
Curl   = getCurlMatrix(Msh)
G      = getNodalGradientMatrix(Msh)
Msig   = getEdgeMassMatrix(Msh,vec(sigma))
Mmu    = getFaceMassMatrix(Msh,1./mu)
Ne,Qe, = getEdgeConstraints(Msh)
Nn,Qn, = getNodalConstraints(Msh)
Nf,Qf, = getFaceConstraints(Msh)

Curl = Qf*Curl*Ne
Msig = Ne'*Msig*Ne
Mmu   = Nf'*Mmu*Nf
G    = Qe*G*Nn
#P    = Ne'*P
nen  = size(Ne,2)
nnn  = size(Nn,2)

dCdm   = [G'*Ne'*getdEdgeMassMatrix(Msh,-Ne*ew[:,1,1]);
          Ne'*getdEdgeMassMatrix(Msh,1/dt[1]*Ne*(ew[:,1,2]-ew[:,1,1]));
          Ne'*getdEdgeMassMatrix(Msh,1/dt[2]*Ne*(ew[:,1,3]-ew[:,1,2]));
          Ne'*getdEdgeMassMatrix(Msh,1/dt[3]*Ne*(ew[:,1,4]-ew[:,1,3]));
          Ne'*getdEdgeMassMatrix(Msh,1/dt[4]*Ne*(ew[:,1,5]-ew[:,1,4]));]

blnkn = spzeros(nnn,4*nen)
blnkee = spzeros(nen,nen)
blnken = spzeros(nen,nnn)
K     = Curl'*Mmu*Curl
dCdu  = [G'*Msig*G blnkn;
         1/dt[1]*Msig*G K+1/dt[1]*Msig blnkee blnkee blnkee;
         blnken -1/dt[2]*Msig K+1/dt[2]*Msig blnkee blnkee;
         blnken blnkee -1/dt[3]*Msig K+1/dt[3]*Msig blnkee;
         blnken blnkee blnkee -1/dt[4]*Msig K+1/dt[4]*Msig;]


Pj     = blkdiag(-P*Ne*G,P*Ne,P*Ne,P*Ne,P*Ne)
z      = rand(Msh.nc)
JzMat  = -Pj*(dCdu\(dCdm*z))
println("Getting MaxwellTime sens mat vec product")
JzStep = getSensMatVec(z,sigma,pFor)
errInf = norm(JzMat-JzStep,Inf)/norm(JzMat,Inf)
errSL2 = norm(JzMat-JzStep)/norm(JzMat)
@test errSL2 < 1e-9
println("Relative Inf and L2 norm errors in Jz are $errInf and $errSL2")

#Transpose
x       = rand(5*size(P,1))
tmp     = dCdu'\(Pj'*x)
JtxMat  = -dCdm'*tmp
println("Getting MaxwellTime sens transpose mat vec product")
JtxStep = getSensTMatVec(x,sigma,pFor)
errInfT = norm(JtxMat-JtxStep,Inf)/norm(JtxMat,Inf)
errSL2T = norm(JtxMat-JtxStep)/norm(JtxMat)
@test errSL2T < 1e-9
println("Relative Inf and L2 norm errors in Jtx are $errInfT and $errSL2T")

#--------------------------------------------------------------------------

#Test multiple source survey
println("Testing multiple sources")

#setup e-dipole sources
ns = 2
Tx = zeros(2,3,ns)
Sources = zeros(ne,ns)
Tx[:,:,1] = [1792.0 2048.0 1024.0;
             2048.0 2048.0 1024.0;]
Tx[:,:,2] = [2048.0 2048.0 1024.0;
             2304.0 2048.0 1024.0;]

for i=1:ns
  Sources[:,i] = getEdgeIntegralOfPolygonalChain(Msh,Tx[:,:,i])
end

#Get data at initial model
obsTimes = [0.0;1.5e-4;2.5e-4;5e-4]
pFor     = getMaxwellTimeParam(Msh,Sources,P',obsTimes,t0,dt,wave,sourceType)
println("Getting data")
D,pFor   = getData(sigma,pFor)


println(" ")
println("==========  Test sigma inversion ======================")
println(" ")

println(" ")
println("==========  Derivative Test ======================")
println(" ")

function f(sigdum)
  d, = getData(sigdum,pFor)
  return d
end

df(zdum,sigdum) = getSensMatVec(zdum,sigdum,pFor)
pass,Error,Order = checkDerivativeMax(f,df,sigma;nSuccess=4)
@test pass


println(" ")
println("==========  Adjoint Test ======================")
println(" ")

v  = randn(prod(size(D)))
tic()
Jz = getSensMatVec(z,sigma,pFor)
toc()
I1 = dot(v,Jz)

tic()
JTz = getSensTMatVec(v,sigma,pFor)
toc()

I2 = dot(JTz,z)

println(I1,"      ",I2)
println("Relative error:",abs(I1-I2)/abs(I1))
@test abs(I1-I2)/abs(I1) < 1e-10



println(" ")
println("==========  Test mu inversion ======================")
println(" ")

m  = MaxwellTimeModel(Dict("sigmaCell"=>sigma,"muCell"=>chi),
                      ["muCell"])
mInit  = MaxwellTimeModel(Dict("sigmaCell"=>sigma,"muCell"=>mu),
                        ["muCell"])
m0 = MaxwellTimeModel(Dict("sigmaCell"=>sigma),Array{String,1}())
obsTimes   = cumsum(dt[2:end])
pFor       = getMaxwellTimeParam(Msh,Sources,P',obsTimes,t0,dt,wave,sourceType)
d,pFor = getData(mInit,pFor)

function f2(sigdum)
  d, = getData(sigdum,pFor)
  return d
end

df2(zdum,sigdum) = getSensMatVec(zdum,sigdum,pFor)

println(" ")
println("==========  Derivative Test ======================")
println(" ")

pass,Error,Order = checkDerivativeMax(f2,df2,m,m0;nSuccess=4)
@test pass


println(" ")
println("==========  Adjoint Test ======================")
println(" ")

v     = randn(prod(size(d)));
dmudm = spdiagm(fill(pi*4e-7,length(mu)));
dm    = MaxwellTimeModelDerivative(Dict(
        "muCell"=>dmudm),["muCell"])

zmod = dm*z
# zmod = MaxwellTimeModel(Dict("sigmaCell"=>sigma,
#                        "muCell"=>dmudm*z),["muCell"])
m2   = MaxwellTimeModel(Dict("sigmaCell"=>sigma,
                        "muCell"=>mu0*(1+chi)),["muCell"])
tic()
Jz = getSensMatVec(zmod,m2,pFor)
toc()
I1 = dot(v,Jz)

tic()
JTzLoc = getSensTMatVec(v,m2,pFor)
JTz    = dm'*JTzLoc
toc()

I2 = dot(JTz,z)

println(I1,"      ",I2)
println("Relative error:",abs(I1-I2)/abs(I1))
@test abs(I1-I2)/abs(I1) < 1e-10

println(" ")
println("==========  Test simultaneous sigma and mu inversion ======================")
println(" ")

activeProps = ["sigmaCell", "muCell"]
m    = MaxwellTimeModel(Dict("sigmaCell"=>log.(sigma),
                        "muCell"=>chi),activeProps)
m0   = MaxwellTimeModel()
pFor = getMaxwellTimeParam(Msh,Sources,P',obsTimes,t0,dt,wave,sourceType)

function f3(sigdum)
  d, = getData(sigdum,pFor)
  return d
end

df3(zdum,sigdum) = getSensMatVec(zdum,sigdum,pFor)

println(" ")
println("==========  Derivative Test ======================")
println(" ")


pass,Error,Order = checkDerivativeMax(f3,df3,m,m0;nSuccess=4)
@test pass


println(" ")
println("==========  Adjoint Test ======================")
println(" ")

dsigdm = [spdiagm(sigma) spzeros(Msh.nc,Msh.nc)]
dmudm  = [spzeros(Msh.nc,Msh.nc) spdiagm(fill(pi*4e-7,length(mu)))]
dm     = MaxwellTimeModelDerivative(Dict(
           "sigmaCell"=>dsigdm,"muCell"=>dmudm),activeProps)

z2 = rand(2*Msh.nc)

zmod = dm*z2
# zmod = MaxwellTimeModel(Dict("sigmaCell"=>dsigdm*z2,
#          "muCell"=>dmudm*z2),activeProps)
m2   = MaxwellTimeModel(Dict("sigmaCell"=>sigma,
         "muCell"=>mu0*(1+chi)),activeProps)
tic()
Jz = getSensMatVec(zmod,m2,pFor)
toc()
I1 = dot(v,Jz)

tic()
JTzLoc = getSensTMatVec(v,m2,pFor)
JTz    = dm'*JTzLoc
toc()

I2 = dot(JTz,z2)

println(I1,"      ",I2)
println("Relative error:",abs(I1-I2)/abs(I1))
@test abs(I1-I2)/abs(I1) < 1e-10

end
