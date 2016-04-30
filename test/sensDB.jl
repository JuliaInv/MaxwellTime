function sensDB(z::Vector{Float64},sigma::Vector{Float64},param::MaxwellTimeParam)	

	# magnetic permeability
	mu = 4*pi*1e-7 
	# extract matrices
	P    = param.Obs
	b0   = param.Sources
	dt   = param.dt
	wave = param.wave
	
	
	Curl = getCurlMatrix(param.M)
	
	Msig = getEdgeMassMatrix(param.M,vec(sigma))
	Mmu  = getFaceMassMatrix(param.M,vec(zeros(size(sigma)).+1/mu))
   N    = getEdgeConstraints(param.M)
	
	Curl = Curl*N
	Msig = N'*Msig*N

	nt = length(dt)
	ns = size(param.Sources,2)
	ne = size(Curl,2)

	ew  = reshape(param.fields,ne,ns,nt+1);
	# Solve linear system
	lam = zeros(Float64,ne,ns,nt+1)
	Jv  = zeros(Float64,size(param.Obs,2),ns,nt)
	 
	for i=1:nt

	    # The linear system to be solved
	    Ke = Curl'*Mmu*Curl + 1/dt[i]*Msig
    	 for j = 1:ns
   	 	Gzi = (1/dt[i])*getdEdgeMassMatrix(param.M,N*(ew[:,j,i+1]-float(i>1)*ew[:,j,i]))
   	 	
	    	rhs = N'*Gzi*z + 1/dt[i]*Msig*lam[:,j,i]
	    	lam[:,j,i+1] = Ke\rhs
			#println(norm(Ke*lam[:,j,i+1]-rhs)/norm(rhs))
	    	# compute Jv
	   	Jvi     = P'*(N*lam[:,j,i+1])
	    	Jv[:,j,i] = -Jvi
		end        
	end

	Jv = Jv[:]

	return Jv
end