#Test sensitivities for BDF2 forward code with OcTree Mesh and grounded 
#source. Single source and multiple source configurations tested. Constant
#time-step size used.

using MaxwellTime
using JOcTree
using Base.Test
using jInv.LinearSolvers
using jInv.Utils
using MUMPS
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
Iobs = round(Integer,ceil(Msh.ne[1]*rand(16)))
P    = speye(ne)
P    = P[Iobs,:]
 
#Time stepping
dt      = 1e-4
nt      = 4
t       = [0.0;cumsum(dt*ones(nt))]
wave    = zeros(nt+1)
wave[1] = 1.0
# 
# #Random model
a     = rand(Msh.nc)
b     = 10*rand(Msh.nc) - 4
sigma = a.^b
 
#Get data at initial model
pFor   = getMaxwellTimeBDF2Param(Msh,Sources,P',dt,nt,wave)
println("Getting data")
d,pFor = getData(sigma,pFor)
ew     = pFor.fields
ehat   = pFor.ehat
#u0     = [phi0;vec(ew[:,1,2:end])];

#Form sensitivity matrix and related quantities
println("Forming sensitivity matrix")
mu    = 4*pi*1e-7  # magnetic permeability
Curl  = getCurlMatrix(Msh)
G     = getNodalGradientMatrix(Msh)
Msig  = getEdgeMassMatrix(Msh,vec(sigma))
Mmu   = getFaceMassMatrix(Msh,vec(zeros(size(sigma)).+1/mu))
Ne,Qe = getEdgeConstraints(Msh)
Nn,Qn = getNodalConstraints(Msh)
Nf,Qf = getFaceConstraints(Msh)

Curl = Qf*Curl*Ne
Msig = Ne'*Msig*Ne
Mmu  = Nf'*Mmu*Nf
G    = Qe*G*Nn
nen  = size(Ne,2)
nnn  = size(Nn,2)

dCdm   = [G'*Ne'*getdEdgeMassMatrix(Msh,-Ne*ew[:,1,1]);
          Ne'*getdEdgeMassMatrix(Msh,3/(2*dt)*Ne*(ehat[:,1]-ew[:,1,1]));
          Ne'*getdEdgeMassMatrix(Msh,1/dt*Ne*(1.5*ew[:,1,2]-0.75*ehat[:,1]-0.75*ew[:,1,1]));
          Ne'*getdEdgeMassMatrix(Msh,1/dt*Ne*(3/2*ew[:,1,3]-2*ew[:,1,2]+1/2*ew[:,1,1]));
          Ne'*getdEdgeMassMatrix(Msh,1/dt*Ne*(3/2*ew[:,1,4]-2*ew[:,1,3]+1/2*ew[:,1,2]));
          Ne'*getdEdgeMassMatrix(Msh,1/dt*Ne*(3/2*ew[:,1,5]-2*ew[:,1,4]+1/2*ew[:,1,3]));]

blnkn = spzeros(nnn,5*nen)
blnkee = spzeros(nen,nen)
blnken = spzeros(nen,nnn)
K     = Curl'*Mmu*Curl
A     = K + 3/(2*dt)*Msig
dCdu  = [G'*Msig*G blnkn;
         3/(2*dt)*Msig*G A blnkee blnkee blnkee blnkee;
         3/(4*dt)*Msig*G -3/(4*dt)*Msig A blnkee blnkee blnkee;
         -1/(2*dt)*Msig*G blnkee -2/dt*Msig A blnkee blnkee;
         blnken blnkee 1/(2*dt)*Msig -2/dt*Msig A blnkee;
         blnken blnkee blnkee 1/(2*dt)*Msig -2/dt*Msig A]
#          
# q = [G'*Ne'*Sources;zeros(4*size(ew,1))]

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
JzMat  = -Pj*solveMUMPS(dCdu,dCdm*z)
println("Getting MaxwellTime sens mat vec product")
JzStep = getSensMatVec(z,sigma,pFor)
errInf = norm(JzMat-JzStep,Inf)/norm(JzMat,Inf)
errSL2 = norm(JzMat-JzStep)/norm(JzMat)
@test errSL2 < 1e-9
println("Rel. Inf and L2 norm errors in Jz are $errInf and $errSL2")

#Transpose
x       = rand(5*size(P,1))
tmp     = -solveMUMPS(dCdu',Pj'*x)
JtxMat  = dCdm'*tmp
println("Getting MaxwellTime sens transpose mat vec product")
#JtxStep,JtxDebug,rhs = getSensTMatVec(x,sigma,pFor)
JtxStep = getSensTMatVec(x,sigma,pFor)
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
pFor   = pFor   = getMaxwellTimeBDF2Param(Msh,Sources,P',dt,nt,wave)
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