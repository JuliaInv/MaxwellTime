using MaxwellTime
using JOcTree
using Base.Test
using jInv.LinearSolvers
using jInv.Utils
using MUMPS

include("maxTime-derivative-test.jl")

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
Iobs = round(Integer,ceil(Msh.ne[1]*rand(16)))
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
chi   = 50*rand(Msh.nc)
mu    = mu0*(1+chi)
m0    = MaxwellTimeModel(sigma,mu,false,true)

#Get data at initial model
sourceType = :Galvanic
obsTimes   = cumsum(dt)
pFor       = getMaxwellTimeParam(Msh,Sources,P',obsTimes,dt,wave,sourceType)

# dummy functions for sensitivity test
function f(mudum)
  m = MaxwellTimeModel(sigma,mudum,false,true)
  d, = getData(m,pFor)
  return d
end

println(" ")
println("==========  Derivative Test ======================")
println(" ")

pass,Error,Order = checkDerivativeMaxMu(f,m0,pFor;nSuccess=4)
@test pass

# println(" ")
# println("==========  Adjoint Test ======================")
# println(" ")
# 
# v     = randn(size(Sources,2)*size(P,1)*length(dt));
# z     = rand(Msh.nc)
# Dsig  = Vector();
# push!(Dsig,UniformScaling(1.0));
# push!(Dsig,UniformScaling(1.0));
# 
# tic()
# Jz = getSensMatVec(m0,sigma,pFor)
# toc()
# I1 = dot(v,Jz)
# 
# tic()
# JTzLoc = getSensTMatVec(v,m0,pFor)
# JTz    = JTzLoc.mu
# toc()
# 
# I2 = dot(JTz,m0.mu)
# 
# println(I1,"      ",I2)
# println("Relative error:",abs(I1-I2)/abs(I1))
# @test abs(I1-I2)/abs(I1) < 1e-10

