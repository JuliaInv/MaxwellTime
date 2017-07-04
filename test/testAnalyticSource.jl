function VMDdhdz(X,x0,rLoop,sigma,t)
  mu    = pi*4e-7
  theta = sqrt(mu*sigma./(4*t))
  rho   = sqrt((X[:,1]-x0[1]).^2 + (X[:,2] - x0[2]).^2)
  inSide = find(rho .<= rLoop)
  rho[inSide] = rLoop + 1.0
  dhdz  = zeros(length(rho),length(t))
  for i=1:length(t)
    dhdz[:,i] = -1./(2*pi*mu*sigma*rho.^5).*( 9*erf(theta[i]*rho) -
                     2*theta[i]*rho/sqrt(pi).*(9+6*theta[i]^2*rho.^2 +
                     4*theta[i]^4*rho.^4).*exp(-theta[i]^2.*rho.^2) )
  end
  return dhdz
end

using MaxwellTime
using JOcTree
using PyPlot
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
        -300.0, 300.0, -300.0, 300.0, -200.0, 200.0,
        1, 2)
M = getOcTreeMeshFV(S,h*ones(3);x0=x0)
Xn = getNodalGrid(M)
Ex, Ey, Ez = getEdgeGrids(M)
E = [Ex;Ey;Ez]
Fx,Fy,Fz = getFaceGrids(M)
nEx,nEy,nEz = getEdgeNumbering(M)
l = 5.0
Tx = [  0.0  0.0 0.0;
        0.0    l 0.0;
          l    l 0.0;
          l  0.0 0.0;
        0.0  0.0 0.0]
MeS = getEdgeIntegralOfPolygonalChain(M,Tx)
a = sqrt((l^2)/pi)
Aloop = staticCurrentLoopVectorPotential(M,a,[l/2;l/2;0.0])
# Ne,Qe, = getEdgeConstraints(M)
# Aloop = Ne'*Aloop
# K = MaxwellTime.getMaxwellCurlCurlMatrix(M,fill(pi*4e-7,M.nc))
# sLoop = K*Aloop

Rx = copy(Tx)
offset = 20.0
Rx[:,1] = Rx[:,1] + offset
P = getEdgeIntegralOfPolygonalChain(M,Rx)

dt    = (5e-4)*ones(40)
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
sigBg  = 0.1
sigAir = 1e-8
sigma  = sigBg*ones(M.nc)
Iair   = find(Xc[:,3] .> 0.0)
sigma[Iair] = sigAir

d1,pForDisc = getData(sigma,pForDisc)
d2,pForAnal = getData(sigma,pForAnal)
d2 = -d2

#Analytic fields
mu = pi*4e-7
dhdz = VMDdhdz([offset+l/2 l/2 0.0],[l/2 l/2 0.0],a,sigBg,obsTimes)
dhdz = dhdz'
dhdz = mu*dhdz

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
loglog(obsTimes[p3],dhdz[p3],"g-")
loglog(obsTimes[m3],dhdz[m3],"g--")
loglog(obsTimes[p3],(1/mu)*dhdz[p3],"m-")
loglog(obsTimes[m3],(1/mu)*dhdz[m3],"m--")
