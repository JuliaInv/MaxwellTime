export getDataAlt

function getDataAlt(sigma::Array{Float64,1},param::MaxwellTimeParam)
#function getData(sigma,param)

# Solve the Backward Euler system
# Curl*e + b_t = 0
# Curl'*(Mmuinv*b) - Msig*e = s(t)
#
# By eliminating b and solving for e
# b_t =  -Curl*e
# (Curl'*Mmuinv*Curl + 1/dt*Msig)*e_{n+1} = - 1/dt*(Msig*e_{n} + s_n - s_{n+1})
#

	P        = param.Obs
	b0       = param.Sources
	dt       = param.dt
	wave     = param.wave
	Msh      = param.M
	mySolver = param.solver
	

	mu     = 4*pi*1e-7	
	Curl = getCurlMatrix(Msh)
	Msig = getEdgeMassMatrix(Msh,vec(sigma))
	Mmu  = getFaceMassMatrix(Msh,1/mu*ones(size(sigma)))
        Ne   = getEdgeConstraints(Msh)
        G    = getNodalGradientMatrix(Msh)
        Nn   = getNodalConstraints(Msh)
      
	# allocate space for fields
        numSrc  = size(b0,2)
	ew = zeros(size(Ne,2),numSrc,length(dt)+1);
	e0   = Msig\(Curl'*Mmu*b0)
	ew[:,:,1] = wave[1]*Ne'*e0

	
	# time step
	dtLast = 0.0
	A = []
	for i=1:length(dt)
            dtinv = 1.0/dt[i]
            rhs = dtinv*Ne'*Msig*Ne*ew[:,:,i]
            if (dt[i] != dtLast)
              mySolver.doClear = 1
              A = Ne'*(Curl'*Mmu*Curl + dtinv*Msig)*Ne
	      ew[:,:,i+1],mySolver = solveLinearSystem(A,rhs,mySolver,1,0)
	      mySolver.doClear = 0
	    else
	      ew[:,:,i+1],mySolver = solveLinearSystem(A,rhs,mySolver,1,0)
            end
	    dtLast = dt[i]
	end
	clear!(mySolver)
        mySolver.doClear = 1
	
	# compute the data
	D = zeros(size(P,2),numSrc,length(dt))
	println("size D is $(size(D)),size ew is $(size(ew))")
	for i=1:length(dt)
		D[:,:,i] = P'*(Ne*ew[:,:,i+1])
	end
	param.fields = ew
	
	return D, param
end