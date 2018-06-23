#Test sensitivities for BDF2 forward code with OcTree Mesh and grounded
#source. Single source and multiple source configurations tested. Constant
#time-step size used.

@testset "Grd src BDF2 variable step-size" begin

using MaxwellTime
using JOcTree
using Base.Test
using jInv.LinearSolvers
using jInv.Utils
include("getC.jl")
# println("Testing single source")

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
#
# #First test that sensitivities from getSensMatVec and getSensTMatVec
# #Match those computed from explicitly constructing J.
Sources = zeros(ne)
Tx = [1792.0 2048.0 1024.0;
      2048.0 2048.0 1024.0;]

Sources = getEdgeIntegralOfPolygonalChain(Msh,Tx)

#Setup receivers
Iobs = round.(Integer,ceil.(Msh.ne[1]*rand(16)))
P    = speye(ne)
P    = P[Iobs,:]

#Time stepping
dt0     = 1e-4
nt      = 4
#dt      = dt0*ones(nt)
dt      = dt0*cumprod([1.0;1.25*ones(nt-1)])
# dt      = [dt[1]*ones(3);dt[2]*ones(3);dt[3]*ones(3);dt[4]*ones(3);]
#dt       = dt0*[1.0;1.0;1.0;1.25^2]
t       = [0;cumsum(dt)]
t0      = t[1]
wave    = zeros(length(dt)+1)
wave[1] = 1.0
#
# #Random model
a     = rand(Msh.nc)
b     = 10*rand(Msh.nc) - 4
sigma = a.^b

#Get data at initial model
sourceType            = :Galvanic
timeIntegrationMethod = :BDF2
obsTimes              = t
pFor                  = getMaxwellTimeParam(Msh,Sources,P',obsTimes,t0,dt,wave,sourceType,
                                            timeIntegrationMethod=timeIntegrationMethod,
                                            EMsolverType=:Pardiso)
pFor.cgTol = 1e-15
println("Getting data")
d,pFor = getData(sigma,pFor)
ew     = pFor.fields
ehat   = pFor.AuxFields

#Form sensitivity matrix and related quantities
println("Forming sensitivity matrix")
mu     = 4*pi*1e-7  # magnetic permeability
Curl   = getCurlMatrix(Msh)
G      = getNodalGradientMatrix(Msh)
Msig   = getEdgeMassMatrix(Msh,vec(sigma))
Mmu    = getFaceMassMatrix(Msh,vec(zeros(size(sigma)).+1/mu))
Ne,Qe, = getEdgeConstraints(Msh)
Nn,Qn, = getNodalConstraints(Msh)
Nf,Qf, = getFaceConstraints(Msh)

Curl = Qf*Curl*Ne
Msig = Ne'*Msig*Ne
Mmu  = Nf'*Mmu*Nf
G    = Qe*G*Nn
nen  = size(Ne,2)
nnn  = size(Nn,2)

blnkn = spzeros(nnn,5*nen)
blnkee = spzeros(nen,nen)
blnken = spzeros(nen,nnn)
K     = Curl'*Mmu*Curl
tau   = dt[2:end]./dt[1:end-1]
g12   = (1 + 2*tau[1])/(1+tau[1])
g13   = (1 + 2*tau[2])/(1+tau[2])
g14   = (1 + 2*tau[3])/(1+tau[3])
g22   = 1 + tau[1]
g23   = 1 + tau[2]
g24   = 1 + tau[3]
g32   = (tau[1]^2)/(1+tau[1])
g33   = (tau[2]^2)/(1+tau[2])
g34   = (tau[3]^2)/(1+tau[3])
A1    = K + 3/(2*dt[1])*Msig
A2    = K + g12/dt[2]*Msig
A3    = K + g13/dt[3]*Msig
A4    = K + g14/dt[4]*Msig
dCdu  = [G'*Msig*G blnkn;
         3/(2*dt[1])*Msig*G A1 blnkee blnkee blnkee blnkee;
         3/(4*dt[1])*Msig*G -3/(4*dt[1])*Msig A1 blnkee blnkee blnkee;
         -g32/(dt[2])*Msig*G blnkee -g22/dt[2]*Msig A2 blnkee blnkee;
         blnken blnkee g33/(dt[3])*Msig -g23/dt[3]*Msig A3 blnkee;
         blnken blnkee blnkee g34/(dt[4])*Msig -g24/dt[4]*Msig A4]

dCdm   = [G'*Ne'*getdEdgeMassMatrix(Msh,-Ne*ew[:,1,1]);
          Ne'*getdEdgeMassMatrix(Msh,3/(2*dt[1])*Ne*(ehat[:,1]-ew[:,1,1]));
          Ne'*getdEdgeMassMatrix(Msh,1/dt[1]*Ne*(1.5*ew[:,1,2]-0.75*ehat[:,1]-0.75*ew[:,1,1]));
          Ne'*getdEdgeMassMatrix(Msh,1/dt[2]*Ne*(g12*ew[:,1,3]-g22*ew[:,1,2]+g32*ew[:,1,1]));
          Ne'*getdEdgeMassMatrix(Msh,1/dt[3]*Ne*(g13*ew[:,1,4]-g23*ew[:,1,3]+g33*ew[:,1,2]));
          Ne'*getdEdgeMassMatrix(Msh,1/dt[4]*Ne*(g14*ew[:,1,5]-g24*ew[:,1,4]+g34*ew[:,1,3]));]

phi0 = (G'*Msig*G)\(G'*Ne'*Sources)
u0   = [phi0;ehat;vec(ew[:,1,2:end])];
q = [G'*Ne'*Sources;zeros(5*size(ew,1))]



# println("--------test dCdm---------------")
# v = rand(size(sigma))*1e-2; v[sigma.<0.009] = 0;
# f = C(u0,sigma,pFor);
# for i=1:10
#   h = 10.0^(-i);
#   fp = C(u0,sigma+h*v,pFor);
#   diff1 = norm(fp-f);
#   diff2 = norm(fp-f - h*dCdm*v);
#   println(h,"   ",diff1,"      ",diff2)
# end
#
#
Pj     = blkdiag(-P*Ne*G,[spzeros(size(P,1),nen) P*Ne],P*Ne,P*Ne,P*Ne)
z      = rand(Msh.nc)
JzMat  = -Pj*(dCdu\(dCdm*z))
println("Getting MaxwellTime sens mat vec product")
JzStep = getSensMatVec(z,sigma,pFor)
errInf = norm(JzMat-JzStep,Inf)/norm(JzMat,Inf)
errSL2 = norm(JzMat-JzStep)/norm(JzMat)
@test errSL2 < 1e-9
println("Rel. Inf and L2 norm errors in Jz are $errInf and $errSL2")

#Transpose
x       = rand(5*size(P,1))
tmp     = -dCdu'\(Pj'*x)
JtxMat  = dCdm'*tmp
println("Getting MaxwellTime sens transpose mat vec product")
JtxStep = getSensTMatVec(x,sigma,pFor)
#JtxStep,JtxDebug = getSensTMatVec(x,sigma,pFor)
errInfT = norm(JtxMat-JtxStep,Inf)/norm(JtxMat,Inf)
errSL2T = norm(JtxMat-JtxStep)/norm(JtxMat)
@test errSL2T < 1e-9
println("Rel. Inf and L2 norm errors in Jtx are $errInfT and $errSL2T")

#--------------------------------------------------------------------------

#Test multiple source survey
println("Testing multiple sources")

#setup e-dipole sources
ns = 2
nEx,nEy,nEz = getEdgeNumbering(Msh)
Tx = zeros(2,3,ns)
Sources = zeros(ne,ns)
Tx[:,:,1] = [1792.0 2048.0 1024.0;
             2048.0 2048.0 1024.0;]
Tx[:,:,2] = [2048.0 2304.0 1024.0;
             2304.0 2304.0 1024.0;]

for i=1:ns
  Sources[:,i] = getEdgeIntegralOfPolygonalChain(Msh,Tx[:,:,i])
end

#Get data at initial model
pFor = getMaxwellTimeParam(Msh,Sources,P',obsTimes,t0,dt,wave,sourceType,
                                 timeIntegrationMethod=timeIntegrationMethod)
pFor.cgTol = 1e-15
println("Getting data")
D,pFor = getData(sigma,pFor)

# D = d
# Jz = JzStep

println(" ")
println("==========  Derivative Test ======================")
println(" ")

function f(sigdum)
  d, = getData(sigdum,pFor)
  return d
end

df(zdum,sigdum) = getSensMatVec(zdum,sigdum,pFor)
pass,Error,Order = checkDerivativeMax(f,df,sigma;nSuccess=5)
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

end
