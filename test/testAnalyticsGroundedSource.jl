using JOcTree
using MaxwellTime
using PyPlot
using SpecialFunctions
include("analyticFields.jl")

x0  = -20480.0*ones(3)
h  = 20.0
n = [2048;2048;2048]
S = createOcTreeFromBox(
        x0[1], x0[2], x0[3],
        n[1], n[2], n[3],
           h,    h,    h,
        -500.0, 500.0, -500.0, 500.0, -400.0, 400.0,
        1, 2)
M = getOcTreeMeshFV(S,h*ones(3);x0=x0)
#Xn = getNodalGrid(M)
Ex, Ey, Ez = getEdgeGrids(M)
Tx = [-20.0 0.0 0.0;
       0.0 0.0 0.0;]
MeS = getEdgeIntegralOfPolygonalChain(M,Tx)

#Setup model
Xc = getCellCenteredGrid(M)
sigBg  = 0.01
sigAir = 1e-8
sigma  = sigBg*ones(M.nc)
#Iair   = find(Xc[:,3] .> 0.0)
#sigma[Iair] = sigAir

#Get receivers
Eg = [Ex;Ey;Ez]
#Define Matrix P, where -P*C*e = dbdt_z at observation points
Iobs = find( (100.0 .< Eg[:,1] .< 125.0) .&
            ( Eg[:,2] .== 0.0) & (Eg[:,3] .== 0.0) )
Isurf = Isurf = find( (-1000.0 .< Ex[:,1] .< 1000.0) .& (-1000.0 .< Ex[:,2] .< 1000.0)
                     .& (Ex[:,3] .== 0.0) );
#Isurf  = find(Ex[:,3] .== 0.0)
P = speye(sum(M.ne))
P = P[Iobs,:]

#Time stepping
dt    = (5e-4)*ones(200)#;(1e-3)*ones(10);0.01*ones(10);0.1*ones(19)]
t     = cumsum(dt)
wave  = zeros(length(dt)+1)
wave[1] = 1.0

obsTimes = [0;t]
sourceType = :Galvanic
pFor = getMaxwellTimeParam(M,MeS,P',obsTimes,dt,wave,sourceType)
#pFor = MaxwellTimeParam(M,b0,P',dt,wave,zeros(0),"",solver)
println("Getting synthetic Data")
d,pFor = getData(sigma,pFor)
eComp = squeeze(pFor.fields[:,:,2:end-1],2)
iXplt = 1
expC   = find(d[2:end-1] .> 0.0)
exmC   = find(d[2:end-1] .< 0.0)


#Analytics
x0plt  = [-10.0;0.0;0.0]
ds = 20.0; I = 1.0
eDC    = EDCDipole(Ex,Ey,Ez,x0plt,sigBg,ds,I)
eAnal  = ETDDipole(Ex,Ey,Ez,x0plt,sigBg,ds,I,t)
Xplt   = Eg[Iobs,1]
eObs   = eAnal[Iobs,2:end-1]
expA    = find(eObs[iXplt,:] .> 0.0)
exmA    = find(eObs[iXplt,:] .< 0.0)

# it = 1;
# figure()
# plot(Xplt[find(dhdz[:,it].>0.0),1],log10(dhdz[find(dhdz[:,it] .> 0.0),it]),"x")
# plot(Xplt[find(dhdz[:,it].<0.0),1],log10(-dhdz[find(dhdz[:,it] .< 0.0),it]),"*")
# figure()hz
# plot(Xplt[find(d[:,it].>0.0),1],log10(d[find(d[:,it] .> 0.0),it]),"x")
# plot(Xplt[find(d[:,it].<0.0),1],log10(-d[find(d[:,it] .< 0.0),it]),"*")
figure()
loglog(t[expA],eObs[iXplt,expA],"b")
loglog(t[exmA],-eObs[iXplt,exmA],"b--")
loglog(t[expC],d[expC],"go-")
loglog(t[exmC],-d[exmC],"g--")
