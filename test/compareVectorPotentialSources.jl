using MaxwellTime
using JOcTree
using PyPlot
using Base.Test
import CGI.MaxwellTime
using CGI.MaxwellUtils
include("analyticFields.jl")
include("/home/patrick/Dropbox/Sample/juliaSetup/importTEMobservationsUBC.jl")
include("/home/patrick/Dropbox/Sample/juliaSetup/pForSetupFuncs.jl")

L   = [10240.0;10240.0;10240.0] #For multiples of 10m cells
#L = [30720;30720;30720]
x0  = -5120.0*ones(3) # For multiples of 10m cells
#x0 = -15360.0*ones(3)
h  = 10.0
n = [2048; 2048; 2048]
# n = [4096;4096;4096]
S = createOcTreeFromBox(
        x0[1], x0[2], x0[3],
        n[1], n[2], n[3],
           h,    h,    h,
        -100.0, 100.0, -100.0, 100.0, -100.0, 100.0,
        1, 2)
M = getOcTreeMeshFV(S,h*ones(3);x0=x0)
Xn = getNodalGrid(M)
Ex, Ey, Ez = getEdgeGrids(M)
E = [Ex;Ey;Ez]
Fx,Fy,Fz = getFaceGrids(M)
nEx,nEy,nEz = getEdgeNumbering(M)
rLoop = 13.5; areaLoop = pi*rLoop^2
rModel = rLoop; lModel = sqrt(areaLoop); lLoop = sqrt(areaLoop)
# txLoc = [h/2; h/2; 0.0]
txLoc = [3.35; 8.2; 0.0]
Tx = [  txLoc[1]-lModel/2  txLoc[2]-lModel/2 txLoc[3];
        txLoc[1]-lModel/2  txLoc[2]+lModel/2 txLoc[3];
        txLoc[1]+lModel/2  txLoc[2]+lModel/2 txLoc[3];
        txLoc[1]+lModel/2  txLoc[2]-lModel/2 txLoc[3];
        txLoc[1]-lModel/2  txLoc[2]-lModel/2 txLoc[3]]


nv = 4
MeS = getEdgeIntegralOfPolygonalChain(M,Tx)

Aloop = staticCurrentLoopVectorPotential(M,rModel,txLoc)
#ADave = get_magnetic_vec_potential_loop(txx0, M, rLoop)
#Adiff = norm(Aloop-ADave)/norm(ADave)
# println("Error norm is $Adiff")
#
# plot(Aloop,"b")
# plot(ADave,"r")

# Ne,Qe, = getEdgeConstraints(M)
# Nf,Qf, = getFaceConstraints(M)
# C      = getCurlMatrix(M)
# C      = Qf*C*Ne
# b0     = C*Ne'*Aloop
# Div    = getDivergenceMatrixRec(M)
# Div    = Div*Nf
# @test norm(Div*b0) ≈ 0.0 atol = eps()
# K = MaxwellTime._getMaxwellCurlCurlMatrix(M,fill(pi*4e-7,M.nc))
mu     = fill(pi*4e-7,M.nc)
Curl   = getCurlMatrix(M)
Tf     = eltype(Curl)
Tn     = eltype(Curl.colptr)
Tn2    = eltype(M.S.SV.nzind)
Mmu    = getFaceMassMatrix(M,one(Tf)./mu)
Nf,Qf, = getFaceConstraints(M)
Ne,Qe, = getEdgeConstraints(M)
Curl   = Qf*Curl*Ne
Mmu    = Nf'*Mmu*Nf
K      = Curl'*Mmu*Curl
#sLoop = -K*Ne'*Aloop
sLoop = -K*Qe*Aloop
#sDave = -Qe*get_src_vec_potential_loop(txx0, M, rLoop)
Rx = [txLoc[1]-h/2  txLoc[2]-h/2 0.0;
      txLoc[1]-h/2  txLoc[2]+h/2 0.0;
      txLoc[1]+h/2  txLoc[2]+h/2 0.0;
      txLoc[1]+h/2  txLoc[2]-h/2 0.0;
      txLoc[1]-h/2  txLoc[2]-h/2 0.0]
Rx2 = copy(Tx)
P  = getEdgeIntegralOfPolygonalChain(M,Rx,normalize=true)
# P2 = getEdgeIntegralOfPolygonalChain(M,Rx2,normalize=true)

# new source
P2 = CGI.MaxwellTime.getZReceiverFaceInterpolation(M,txLoc)
nr = 1; ne = sum(M.ne)
Pnew = zeros(ne)

W, pts = CGI.MaxwellTime.getTrilinWeightsCC(txLoc, M)
for i = 1:8
    Wi = W[i]
    pti = pts[i, :]
    Pnew .-= Wi*CGI.MaxwellTime.get_dipole_src_A_potential(M, pti, [0.,0.,1.])
end

# dt    = (2.5e-6)*ones(120)
dt    = [(1.25e-6)*ones(20); (2.5e-6)*ones(20); (5e-6)*ones(60)]
t     = cumsum(dt)
t0    = 0.0
wave  = zeros(length(dt)+1)
wave[1] = 1.0
obsTimes = cumsum(dt[1:end-1])

sourceType = :InductiveLoopPotential
pForDisc   = getMaxwellTimeParam(M,sLoop,P2,obsTimes,t0,dt,wave,sourceType;
                                 timeIntegrationMethod=:BE,
                                 verbosity=2,
                                 storageLevel=:None)

sourceType = :InductiveLoopPotential
pForLoop   = getMaxwellTimeParam(M,sLoop,P,obsTimes,t0,dt,wave,sourceType;
                                 timeIntegrationMethod=:BE,
                                 verbosity=2,
                                 storageLevel=:None)

pForDave = getMaxwellTimeParam(M,sLoop,Pnew,obsTimes,t0,dt,wave,sourceType;
                                 timeIntegrationMethod=:BE,
                                 verbosity=2,
                                 storageLevel=:None)

#Setup model
Xc = getCellCenteredGrid(M)
sigBg  = 1e-2
sigAir = 1e-8
sigma  = sigBg*ones(M.nc)
Iair   = find(Xc[:,3] .> 0.0)
sigma[Iair] = sigAir

println("Getting rx=square tx loop data")
dDisc,pForDisc = getData(sigma,pForDisc)
println("Getting rx = 1 cell square loop data")
dLoop,pForLoop = getData(sigma,pForLoop)
println("Getting Dave's analytic loop rx data")
dDave,pForDave = getData(sigma,pForDave)
# dLoop = -dLoop

#Analytic fields
using SpecialFunctions
mu = pi*4e-7
# hz = hzAnalyticCentLoop(a,obsTimes,sigBg)
# bzAnal = mu*hz
dbzdt  = mu*VMDdhdtzCentLoop(rModel,obsTimes,sigBg)
#dbzdt  = (1/areaLoop)*mu*VMDdhdtzCentLoop(rLoop,obsTimes,sigBg)


# plot dbdt
p1 = find(dDisc .> 0)
m1 = find(dDisc .< 0)
p2 = find(dLoop .> 0)
m2 = find(dLoop .< 0)
p3 = find(dbzdt .> 0)
m3 = find(dbzdt .< 0)
p4 = find(dDave .> 0)
m4 = find(dDave .< 0)

figure()
loglog(obsTimes[p1],dDisc[p1],"b-o",label="big rx BDF2")
loglog(obsTimes[m1],abs.(dDisc[m1]),"b--")
loglog(obsTimes[p2],dLoop[p2],"r-",label="small loop rx BDF2")
loglog(obsTimes[m2],abs.(dLoop[m2]),"r--")
loglog(obsTimes[p3],dbzdt[p3],"m-",label="true solution")
loglog(obsTimes[m3],abs.(dbzdt[m3]),"m--")
loglog(obsTimes[p4],dDave[p4],"g-",label="trilinear interp rx BDF2")
loglog(obsTimes[m4],abs.(dDave[m4]),"g--")
xlabel("time (s)"); ylabel("db_z/dt")
legend()

# plot dbdt err
err1 = abs.(dDisc-dbzdt)./dbzdt;
err2 = abs.(dLoop-dbzdt)./dbzdt;
err3 = abs.(dDave-dbzdt)./dbzdt;
figure()
loglog(obsTimes,err1,"b-",label="big rx BDF2")
loglog(obsTimes,err2,"r-",label="small loop rx BDF2")
loglog(obsTimes,err3,"g-",label="trilinear interp rx BDF2")
legend()
