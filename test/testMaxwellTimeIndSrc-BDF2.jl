using MaxwellTime
using Mesh
using Base.Test
using Utils
using MUMPS
using LinearSolvers
using KrylovMethods


L  = [4096, 4096, 2048.0]
x0 = [0.0,  0.0,  0.0]
n  = [8,  8,  8]
h  = L./n

xmin = 1548.0; xmax = 2548.0
ymin = 1548.0; ymax = 2548.0
zmin = 524.0;  zmax = 1524.0


S  = createOcTreeFromBox(x0[1],x0[2],x0[3],n[1],n[2],n[3],h[1],h[2],h[3],xmin,xmax,ymin,ymax,zmin,zmax,1,1)
M  = getOcTreeMeshFV(S,h;x0=x0)

Xe1, Xe2, Xe3 = getEdgeGrids(M)
Xc		      = getCellCenteredGrid(M)


# setup loop source
pc      = [x0[1] + n[1]*h[1]/2, x0[2] + n[2]*h[2]/2, x0[3] + n[3]*h[3]/2]

flightPath = -1000:500:1000
ns         = length(flightPath)
nr         = ns
Curl       = getCurlMatrix(M)
EX, EY, EZ = getEdgeNumbering(M)

Sources = zeros(size(Curl,2),ns)
P       = spzeros(size(Curl,2),ns)

for i=1:ns

   # Define sources
	p       = [pc[1]-10+flightPath[i]  pc[2]-10   pc[3]
				  pc[1]-10+flightPath[i]  pc[2]+10   pc[3]
			  	  pc[1]+10+flightPath[i]  pc[2]+10   pc[3]
			  	  pc[1]+10+flightPath[i]  pc[2]-10   pc[3]
			     pc[1]-10+flightPath[i]  pc[2]-10   pc[3]] 


# 	a1           = magnetostaticsCurrentPoly(Xe1, p)
# 	a2           = magnetostaticsCurrentPoly(Xe2, p)
# 	a3           = magnetostaticsCurrentPoly(Xe3, p)
# 	a            = [a1[:,1]; a2[:,2]; a3[:,3]] 
# 	Sources[:,i] = 1e8*Curl*a

	# Define receivers
	P[:,i]       = getEdgeIntegralOfPolygonalChain(M,p,EX,EY,EZ)
        Sources[:,i] = P[:,i]
end

t       = [1:6;]*1e-4 #[0; logspace(-6,-2,25)] #[0,1.3,2.7,4.5,6.4]*1e-8; #        
dt      = diff(t)
wave    = zeros(length(dt)+1); wave[1] = 1.0

sigma   = zeros(M.nc)+1e-2
sigma[Xc[:,3] .> 1024] = 1e-8

solver=getMUMPSsolver([])
#param = getMaxwellTimeParam(M,Sources,P,dt,wave,solver)
param = MaxwellTimeBDF2Param(M,Sources,P,dt[1],length(dt),wave,zeros(0),"",solver,zeros(1,1))

tic()
D,param = getData(sigma,param);
toc()

println(" ")
println("==========  Derivative Test ======================")
println(" ")

z = rand(size(sigma))*1e-2;z[sigma.<0.009] = 0;
Jz = getSensMatVec(z,sigma,param)
nrmold = [0.0 0.0]
nrm    = [0.0 0.0]
for i=1:10
	alpha = 2.0^(-i)
	D1, = getData(sigma+alpha*z,param)
	nrm = [norm(D1[:]-D[:]) norm(D1[:]-D[:]-alpha*Jz)]
	println(alpha,"      ",nrm[1],"       ",nrm[2],"       ",nrmold[2]/nrm[2])
	nrmold = copy(nrm)
end


println(" ")
println("==========  Adjoint Test ======================")
println(" ")

v  = randn(prod(size(D)))
tic()
Jz = getSensMatVec(z,sigma,param)
toc()
I1 = dot(v,Jz)

tic()
JTz = getSensTMatVec(v,sigma,param)
toc()

I2 = dot(JTz,z)

println(I1,"      ",I2)
println("Relative error:",abs(I1-I2)/abs(I1))
@test abs(I1-I2)/abs(I1) < 1e-5
