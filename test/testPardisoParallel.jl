addprocs(2)
#@everywhere ENV["OMP_NUM_THREADS"] = 2

import MaxwellTime, JOcTree, jInv.Utils

@everywhere using MaxwellTime, JOcTree, jInv.Utils

L   = [10240.0;10240.0;10240.0] #For multiples of 10m cells
#L = [30720;30720;30720]
x0  = -5120.0*ones(3) # For multiples of 10m cells
#x0 = -15360.0*ones(3)
h  = 10.0
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
rModel = sqrt(1/pi); lModel = 1.0; lLoop = sqrt(areaLoop)
txHeight = 0.0
Tx = [  h/2-lModel/2  h/2-lModel/2 txHeight;
        h/2-lModel/2  h/2+lModel/2 txHeight;
        h/2+lModel/2  h/2+lModel/2 txHeight;
        h/2+lModel/2  h/2-lModel/2 txHeight;
        h/2-lModel/2  h/2-lModel/2 txHeight]



MeS = getEdgeIntegralOfPolygonalChain(M,Tx)
Rx = copy(Tx)
# Rx[:,1] .-= 6
# Rx[:,2] .+= 10
P = getEdgeIntegralOfPolygonalChain(M,Rx,normalize=true)

dt    = (5e-6)*ones(200)
t     = cumsum(dt)
t0    = 0.0
wave  = zeros(length(dt)+1)
wave[1] = 1.0
obsTimes = cumsum(dt[1:end-1])

sourceType = :InductiveDiscreteWire
pForRefs = Vector{RemoteChannel}(2)
for i in 1:2
    pForRefs[i] = initRemoteChannel(getMaxwellTimeParam,i+1,M,MeS,P,obsTimes,t0,dt,wave,sourceType;
                                 timeIntegrationMethod=:BE,
                                 EMsolverType=:Pardiso)
end

#Setup model
Xc = getCellCenteredGrid(M)
sigBg  = 1e-2
sigAir = 1e-8
sigma  = sigBg*ones(M.nc)
Iair   = find(Xc[:,3] .> 0.0)
sigma[Iair] = sigAir

@time dDisc,pForRefs = getData(sigma,pForRefs)
