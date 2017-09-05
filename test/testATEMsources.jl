using MaxwellTime
using JOcTree

L = [30720;30720;30720]
#x0  = -5120.0*ones(3) # For multiples of 10m cells
x0 = -15360.0*ones(3)
h  = 15.0
n = 2*[1024;1024;1024]
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
rLoop = 13.5; areaLoop = pi*rLoop^2
lLoop = sqrt(areaLoop)
txHeight = 15.0
Tx = [  h/2-lLoop/2  h/2-lLoop/2 txHeight;
        h/2-lLoop/2  h/2+lLoop/2 txHeight;
        h/2+lLoop/2  h/2+lLoop/2 txHeight;
        h/2+lLoop/2  h/2-lLoop/2 txHeight;
        h/2-lLoop/2  h/2-lLoop/2 txHeight]



MeS = getEdgeIntegralOfPolygonalChain(M,Tx)
Aloop = staticCurrentLoopVectorPotential(M,rModel,[h/2;h/2;txHeight])
Ne,Qe, = getEdgeConstraints(M)
Nf,Qf, = getFaceConstraints(M)
Nn,Qn, = getNodalConstraints(M)
C      = getCurlMatrix(M)
C      = Qf*C*Ne
G      = getNodalGradientRec(M)
G      = Qe*G*Nn
b0     = C*Ne'*Aloop
Div    = getDivergenceMatrixRec(M)
Mmuinv = getFaceMassMatrix(M,fill(pi*4e-7,M.nc))
Mmuinv = Nf'*Mmuinv*Nf
Div    = Div*Nf
@test norm(Div*b0) ≈ 0.0 atol = eps()
s      = C'*Mmuinv*b0
@test norm(G'*s) ≈ 0.0 atol = eps()
# K = MaxwellTime.getMaxwellCurlCurlMatrix(M,fill(pi*4e-7,M.nc))
# sLoop = K*Aloop


Rx = copy(Tx)
P  = getEdgeIntegralOfPolygonalChain(M,Rx,normalize=true)

dt    = (1e-5)*ones(200)
t     = cumsum(dt)
t0    = 0.0
wave  = ones(length(dt)+1)
#wave[1] = 1.0
obsTimes = cumsum(dt[1:end-1])

sourceType = :InductiveDiscreteWire
pForDisc   = getMaxwellTimeParam(M,MeS,P,obsTimes,t0,dt,wave,sourceType;
                                 timeIntegrationMethod=:BE)

sourceType = :InductiveLoopPotential
pForAnal   = getMaxwellTimeParam(M,Aloop,P,obsTimes,t0,dt,wave,sourceType;
                                 timeIntegrationMethod=:BE)

#Setup model
Xc = getCellCenteredGrid(M)
sigBg  = 1e-2
sigAir = 1e-8
sigma  = sigBg*ones(M.nc)
Iair   = find(Xc[:,3] .> 0.0)
sigma[Iair] = sigAir

dDisc,pForDisc = getData(sigma,pForDisc)
dLoop,pForAnal = getData(sigma,pForAnal)
