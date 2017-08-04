

using MaxwellTime
using JOcTree
using PyPlot
using Base.Test
include("analyticFields.jl")
include("/home/patrick/Dropbox/Sample/juliaSetup/importTEMobservationsUBC.jl")
include("/home/patrick/Dropbox/Sample/juliaSetup/pForSetupFuncs.jl")

#L   = [10240.0;10240.0;10240.0] #For multiples of 10m cells
L = [30720;30720;30720]
#x0  = -5120.0*ones(3) # For multiples of 10m cells
x0 = -15360.0*ones(3)
h  = 30.0
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
Tx = [  h/2-l/2  h/2-l/2 0.0;
        h/2-l/2  h/2+l/2 0.0;
        h/2+l/2  h/2+l/2 0.0;
        h/2+l/2  h/2-l/2 0.0;
        h/2-l/2  h/2-l/2 0.0]



MeS = getEdgeIntegralOfPolygonalChain(M,Tx)
a = sqrt((l^2)/pi)
Aloop = staticCurrentLoopVectorPotential(M,a,[h;h;0.0])
Ne,Qe, = getEdgeConstraints(M)
Nf,Qf, = getFaceConstraints(M)
C      = getCurlMatrix(M)
C      = Qf*C*Ne
b0     = C*Ne'*Aloop
Div    = getDivergenceMatrixRec(M)
Div    = Div*Nf
@test norm(Div*b0) ≈ 0.0 atol = eps()
# K = MaxwellTime.getMaxwellCurlCurlMatrix(M,fill(pi*4e-7,M.nc))
# sLoop = K*Aloop
Rx = copy(Tx)
Rx[:,1] .-= 6
Rx[:,2] .+= 10
P = getEdgeIntegralOfPolygonalChain(M,Rx,normalize=true)

dt    = (5e-6)*ones(200)
#dt    = [(5e-6)*ones(25);(1e-5)*ones(25);(2.25e-5)*ones(25)]
t     = cumsum(dt)
wave  = zeros(length(dt)+1)
wave[1] = 1.0
obsTimes = cumsum(dt[1:end-1])

sourceType = :InductiveDiscreteWire
pForBE   = getMaxwellTimeParam(M,MeS,P,obsTimes,dt,wave,sourceType;
                                 timeIntegrationMethod=:BE)

pForBDF2   = getMaxwellTimeParam(M,MeS,P,obsTimes,dt,wave,sourceType;
                                 timeIntegrationMethod=:BDF2)
pForBDF2.cgTol = 1e-14
#Setup model
Xc = getCellCenteredGrid(M)
sigBg  = 1e-2
sigAir = 1e-8
sigma  = sigBg*ones(M.nc)
Iair   = find(Xc[:,3] .> 0.0)
sigma[Iair] = sigAir

d2,pForBDF2 = getData(sigma,pForBDF2)
d1,pForBE   = getData(sigma,pForBE)

#d2 = -d2

#Analytic fields
using SpecialFunctions
mu = pi*4e-7
moment = l^2
hz = hzAnalyticCentLoop(a,obsTimes,sigBg)
bzAnal = mu*hz
dbzdt  = mu*VMDdhdtzCentLoop(a,obsTimes,sigBg)



# Compute bz from approximate dbdtz
Fx,Fy,Fz = getFaceGrids(M)
F = [Fx;Fy;Fz]
itx      = find((0.<F[:,1].<h) .& (0.<F[:,2].<h) .& (F[:,3] .== 0))
function getb(b0,d,dt)
    b = zeros(length(dt)+1)
    b[1] = b0
    for i = 1:length(dt)
        b[i+1] = b[i] + dt[i]*d[i]
    end
    return b
end

b0tx = -(Nf*b0)[itx[1]]
bzDisc = getb(b0tx,d1,dt[1:end-1])
bzLoop = getb(b0tx,d2,dt[1:end-1])
bzDisc = bzDisc[2:end]
bzLoop = bzLoop[2:end]

# Plot b
p1 = find(bzDisc .> 0)
m1 = find(bzDisc .< 0)
p2 = find(bzLoop .> 0)
m2 = find(bzLoop .< 0)
p3 = find(bzAnal .> 0)
m3 = find(bzAnal .< 0)

# figure()
# loglog(obsTimes[p1],bzDisc[p1],"b-")
# loglog(obsTimes[m1],abs.(bzDisc[m1]),"b--")
# loglog(obsTimes[p2],bzLoop[p2],"r-")
# loglog(obsTimes[m2],abs.(bzLoop[m2]),"r--")
# loglog(obsTimes[p3],bzAnal[p3],"m-")
# loglog(obsTimes[m3],abs.(bzAnal[m3]),"m--")

# plot dbdt
p1 = find(d1 .> 0)
m1 = find(d1 .< 0)
p2 = find(d2 .> 0)
m2 = find(d2 .< 0)
p3 = find(dbzdt .> 0)
m3 = find(dbzdt .< 0)

figure()
loglog(obsTimes[p1],d1[p1],"bo-")
loglog(obsTimes[m1],abs.(d1[m1]),"b--")
loglog(obsTimes[p2],d2[p2],"ro-")
loglog(obsTimes[m2],abs.(d2[m2]),"r--")
loglog(obsTimes[p3],dbzdt[p3],"m-")
loglog(obsTimes[m3],abs.(dbzdt[m3]),"m--")

# plot dbdt err
err1 = abs.(d1-dbzdt)./dbzdt;
err2 = abs.(d2-dbzdt)./dbzdt;
figure()
loglog(err1)
loglog(err2)
