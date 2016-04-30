using Mesh
using MaxwellTime
using LinearSolvers
using testFunctions

nm = zeros(4,4)
i0=1; imax = 4
nfine = [imax+1;imax+1;imax+1]
is=i0
Msh = []; S= []; h = []
#for is=i0:imax
n = [2^5;2^5;2^5]

h = 1./n
i,j,k = ndgrid(1:n[1],1:n[2],1:n[3])
bsz   = ones(Int,prod(n))
S     = sparse3(vec(i),vec(j),vec(k),bsz,n)
Msh = getOcTreeMeshFV(S,h)

# if (is==i0)
#   S,h = initCoarseOcTree([0.0 0.0 0.0;1.0 1.0 1.0],nfine,imax+1-i0)
#   Msh = getOcTreeMeshFV(S,h)
#   Xc  = getCellCenteredGrid(Msh)
#   Irf = find( (0.25 .< Xc[:,1] .< 0.75) & 
#            (0.25 .< Xc[:,2] .< 0.75) & (0.25 .< Xc[:,3] .< 0.75) )
#   tau = zeros(Msh.nc)
#   tau[Irf] = 1.0
#   S   = refineOcTree(S,tau,0.5)
#   Msh = getOcTreeMeshFV(S,h)
#   #tau = zeros(nnz(S.SV))
#   #Irf2 = find( (0.375 .< Xc[:,1] . < 0.625) & 
#   #         (0.375 .< Xc[:,2] . < 0.625) & (0.375 .< Xc[:,3] . < 0.625) )
# else
#   S = refineOcTree(S,ones(Msh.nc),0.5)
#   Msh = getOcTreeMeshFV(S,h)
# end

sigma = ones(Msh.nc)
Ex,Ey,Ez = getEdgeGrids(Msh)
Fx,Fy,Fz = getFaceGrids(Msh)
idxEdgeInt = find((0.0 .< [Ex[:,1];Ey[:,1];Ez[:,1]] .< 1.0) & 
           (0.0 .< [Ex[:,2];Ey[:,2];Ez[:,2]] .< 1.0) &
           (0.0 .< [Ex[:,3];Ey[:,3];Ez[:,3]] .< 1.0))
idxEdgeBdry = setdiff([1:sum(Msh.ne);],idxEdgeInt)
# idxFaceInt = find((0.0 .< [Fx[:,1];Fy[:,1];Fz[:,1]] .< 1.0) & 
#            (0.0 .< [Fx[:,2];Fy[:,2];Fz[:,2]] .< 1.0) &
#            (0.0 .< [Fx[:,3];Fy[:,3];Fz[:,3]] .< 1.0))
# idxFaceBdry = setdiff([1:sum(Msh.nf);],idxFaceInt)
Pbdry = speye(sum(Msh.ne))
Pbdry = Pbdry[idxEdgeInt,:]

A   = getA(Ex,Ey,Ez)
#CA  = getCA(Fx,Fy,Fz)
CCA = getCCAalt(Ex,Ey,Ez)
Me  = getEdgeMassMatrix(Msh,ones(Msh.nc))

#Define time dependence
#b   = exp(-t/tau)*CA
#e   = (1/tau)*exp(-t/tau)*A
dt0 = 1.0e-4; rfFac = 8
dt   = (dt0/rfFac)*ones(10*rfFac)
#dt   = dt/(2^(is-1))
t    = [0.0;cumsum(dt)]
tau  = 1.0e-4
wave = exp(-t/tau)
e    = (1/tau)*A*wave'
#b0   = CA
s    = Me*(CCA - (1/tau)*A)

P = speye(sum(Msh.ne)); P = P[:,idxEdgeBdry[1:10]]
solver = getMUMPSsolver([])
pFor   = getMaxwellTimeParam(Msh,s,P,dt,wave,solver)
d,pFor = getData(sigma,pFor,e[:,1])
eComp = squeeze(pFor.fields,2)
println("The norms are $(nm[is-i0+1,1]) and $(nm[is-i0+1,2])")
#end

println("$(nm[1:end-1,1]./nm[2:end,1]) and $(nm[1:end-1,2]./nm[2:end,2])")