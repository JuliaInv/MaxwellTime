

using MaxwellTime
using JOcTree
using PyPlot
include("analyticFields.jl")
include("/home/patrick/Dropbox/Sample/juliaSetup/importTEMobservationsUBC.jl")
include("/home/patrick/Dropbox/Sample/juliaSetup/pForSetupFuncs.jl")

L   = [10240.0;10240.0;10240.0]
x0  = -5120.0*ones(3)
h  = 10.0
n = [1024;1024;1024]
S = createOcTreeFromBox(
        x0[1], x0[2], x0[3],
        n[1], n[2], n[3],
           h,    h,    h,
        -200.0, 200.0, -200.0, 200.0, -200.0, 200.0,
        1, 2)
M = getOcTreeMeshFV(S,h*ones(3);x0=x0)
Xn = getNodalGrid(M)
Ex, Ey, Ez = getEdgeGrids(M)
E = [Ex;Ey;Ez]
Fx,Fy,Fz = getFaceGrids(M)
nEx,nEy,nEz = getEdgeNumbering(M)
l = 2.0
Tx = [  0.0  0.0 0.0;
        0.0    l 0.0;
          l    l 0.0;
          l  0.0 0.0;
        0.0  0.0 0.0]
MeS = getEdgeIntegralOfPolygonalChain(M,Tx)
a = sqrt((l^2)/pi)
Aloop = staticCurrentLoopVectorPotential(M,a,[l/2;l/2;0.0])
Ne,Qe, = getEdgeConstraints(M)
Nf,Qf, = getFaceConstraints(M)
C      = getCurlMatrix(M)
C      = Qf*C*Ne
Aloop  = Ne'*Aloop
b0     = C*Aloop
Div    = getDivergenceMatrixRec(M)
Div    = Div*Nf
# K = MaxwellTime.getMaxwellCurlCurlMatrix(M,fill(pi*4e-7,M.nc))
# sLoop = K*Aloop

Rx = copy(Tx)
offset = 10.0#150.0
Rx[:,1] = Rx[:,1] + offset
P = getEdgeIntegralOfPolygonalChain(M,Rx,normalize=true)

dt    = (5e-5)*ones(100)
t     = cumsum(dt)
wave  = zeros(length(dt)+1)
wave[1] = 1.0
obsTimes = cumsum(dt[1:end-1])

sourceType = :InductiveDiscreteWire
pForDisc   = getMaxwellTimeParam(M,MeS,P,obsTimes,dt,wave,sourceType)

sourceType = :InductiveAnalyticLoop
pForAnal   = getMaxwellTimeParam(M,Aloop,P,obsTimes,dt,wave,sourceType)

#Setup model
Xc = getCellCenteredGrid(M)
sigBg  = 1e-2#0.1
sigAir = 1e-8
sigma  = sigBg*ones(M.nc)
Iair   = find(Xc[:,3] .> 0.0)
sigma[Iair] = sigAir

d1,pForDisc = getData(sigma,pForDisc)
d2,pForAnal = getData(sigma,pForAnal)
#d2 = -d2

#Analytic fields
using SpecialFunctions
mu = pi*4e-7
moment = l^2
dhdz = moment*VMDdhdtz([offset+l/2 l/2 0.0],[l/2 l/2 0.0],a,sigBg,obsTimes)
dhdz = dhdz'
dbdz = mu*dhdz




p1 = find(d1 .> 0)
m1 = find(d1 .< 0)
p2 = find(d2 .> 0)
m2 = find(d2 .< 0)
p3 = find(dhdz .> 0)
m3 = find(dhdz .< 0)

figure()
loglog(obsTimes[p1],d1[p1],"b-")
loglog(obsTimes[m1],d1[m1],"b--")
loglog(obsTimes[p2],d2[p2],"r-")
loglog(obsTimes[m2],d2[m2],"r--")
# loglog(obsTimes[p3],dhdz[p3],"g-")
# loglog(obsTimes[m3],dhdz[m3],"g--")
loglog(obsTimes[p3],dbdz[p3],"m-")
loglog(obsTimes[m3],dbdz[m3],"m--")
