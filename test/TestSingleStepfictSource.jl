using Mesh
using MaxwellTime
using LinearSolvers
using testFunctions

nm = zeros(4,4)
i0=1; imax = 4
nfine = [imax+1;imax+1;imax+1]
is=i0
#Msh = []; S= []; h = []
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
mu    = 1.0

Curl = getCurlMatrix(Msh)

Ex,Ey,Ez = getEdgeGrids(Msh)
Fx,Fy,Fz = getFaceGrids(Msh)
idxEdgeInt = find((0.0 .< [Ex[:,1];Ey[:,1];Ez[:,1]] .< 1.0) & 
           (0.0 .< [Ex[:,2];Ey[:,2];Ez[:,2]] .< 1.0) &
           (0.0 .< [Ex[:,3];Ey[:,3];Ez[:,3]] .< 1.0))
idxEdgeBdry = setdiff([1:sum(Msh.ne);],idxEdgeInt)
idxFaceInt = find((0.0 .< [Fx[:,1];Fy[:,1];Fz[:,1]] .< 1.0) & 
           (0.0 .< [Fx[:,2];Fy[:,2];Fz[:,2]] .< 1.0) &
           (0.0 .< [Fx[:,3];Fy[:,3];Fz[:,3]] .< 1.0))
idxFaceBdry = setdiff([1:sum(Msh.nf);],idxFaceInt)

#Get matrices
Curl = getCurlMatrix(Msh)
Msig = getEdgeMassMatrix(Msh,vec(sigma))
Mmu  = getFaceMassMatrix(Msh,vec(zeros(size(sigma)).+1/mu))
Ne   = getEdgeConstraints(Msh)
G    = getNodalGradientMatrix(Msh)
Nn   = getNodalConstraints(Msh)
Me   = getEdgeMassMatrix(Msh,vec(sigma))

A   = getA(Ex,Ey,Ez)
CA  = getCA(Fx,Fy,Fz)
CCA = getCCAalt(Ex,Ey,Ez)

#Define time dependence
dt  = (1.0e-4)*ones(5)
dt = dt/(2^(is-1))
t   = [0.0;cumsum(dt)]
tau = 1.0e-4
Tb  = exp(-t/tau)
Te  = (1.0/tau)*Tb
e   = zeros(length(A),length(t))
b   = zeros(length(CA),length(t))
js  = zeros(length(A),length(t))
dedt = copy(e)
dsdt = copy(js)
e   = A*Te'
b   = CA*Tb'
js  = CCA*Tb' - A*Te'
dedt = (-1.0/tau^2)*A*Tb'
dsdt = -CCA*Te' - dedt

dtinv = 1/dt[1]; Lm = Ne'*(Curl'*Mmu*Curl + dtinv*Msig)*Ne
rhs = dtinv*(Ne'*Msig*Ne*e[:,1]+(Tb[1]-Tb[2])*Ne'*Me*(CCA-(1/tau)*A))
e2Comp = solveMUMPS(Lm,rhs,1)
nm[is-i0+1,1] = norm(e[:,2]-e2Comp,Inf)
nm[is-i0+1,2] = norm(e[idxEdgeInt,2]-e2Comp[idxEdgeInt],Inf)

rhsSdComp = Ne'*Curl'*Mmu*Curl*Ne*e + Ne'*Msig*Ne*dedt
rhsSdAnal = -Me*dsdt 
nm[is-i0+1,3] = norm(rhsSdComp-rhsSdAnal,Inf)
nm[is-i0+1,4] = norm(rhsSdComp[idxEdgeInt,1]-rhsSdAnal[idxEdgeInt,1],Inf)
println("The norms are $(nm[is-i0+1,1]) and $(nm[is-i0+1,2])")
#end

println("$(nm[1:end-1,1]./nm[2:end,1]) and $(nm[1:end-1,2]./nm[2:end,2])")