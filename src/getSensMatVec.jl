export getSensMatVec

function getSensMatVec(z::Vector{Float64},sigma::Vector{Float64},
                       param::MaxwellTimeParam)
	#getSensMatVec(z,sigma,param)
	#This function computes (dData/dsigma)*z for BE time-stepping
	#forward problem. It handles grounded and inductive sources, 
	#assuming DC data is integral of electric field and not a 
	#potential difference (i.e. electric field at t=0 is known) for  
	#grounded sources and e0=0 for inductive sources.
	#For grounded sources
	# dCdu = |G'*Msig*G                                    |
	#        |1/dt*Msig*G    K+1/dt*Msig                   |
	#        |               -1/dt*Msig    K+1/dt*Msig  ...|
	#        |                     .                       |
	#        |                     .                       |
	#        |                     .                       |
	dt       = param.dt
	wave     = param.wave
	Msh      = param.M
	EMsolver = param.EMsolver
	DCsolver = param.DCsolver
	storeDCf = param.storeDCfactors
	ew       = param.fields
	
	mu     = 4*pi*1e-7  # magnetic permeability
	Curl   = getCurlMatrix(Msh)
	Gin    = getNodalGradientMatrix(Msh)
	Msig   = getEdgeMassMatrix(Msh,vec(sigma))
	Mmu    = getFaceMassMatrix(Msh,vec(zeros(size(sigma)).+1/mu))
        Ne,Qe, = getEdgeConstraints(Msh)
        Nn,Qn  = getNodalConstraints(Msh)
        Nf,Qf  = getFaceConstraints(Msh)
	
	Curl = Qf*Curl*Ne
	Msig = Ne'*Msig*Ne
	Mmu   = Nf'*Mmu*Nf
	G    = Qe*Gin*Nn
	s    = Ne'*copy(param.Sources)
	P    = Ne'*param.Obs

	ns = size(s,2)
	ne = size(Curl,2)
        nt = length(dt)
        
	#Check source type
	groundedSource = any( abs(G'*s) .> 1e-12) ? true : false

	#Initialize intermediate and output arrays
	lam = zeros(Float64,ne,ns,2)
	Jv  = zeros(Float64,size(P,2),ns,nt)
	
	#Do the DC part of dCdm if source is grounded.
	#Note that this code assumes (for grounded sources) that DC data
	#Are computed as the integral of electric field over a receiver
	#dipole and not as a potential difference
	if groundedSource
	  A = G'*Msig*G
	  Jvdc = zeros(Float64,size(P,2),ns)
	  for j = 1:ns
	    Gzi = G'*Ne'*getdEdgeMassMatrix(Msh,-Ne*ew[:,j,1])
	    rhs = Gzi*z
	    lam0,DCsolver    = solveDC!(A,rhs,DCsolver)
	    lam[:,j,1]       = -G*lam0 #Taking gradient of lam0
	                               #Prepares for data projection and
	                               #preps lam for use as rhs in first 
	                               #time-step
	    DCsolver.doClear = 0
	    Jvdc[:,j] = -P'*lam[:,j,1]
	  end
	  if ~storeDCf
	    clear!(DCsolver)
	    DCsolver.doClear = 1
	  end
	end
	
	nFacs = 0
	for i=1:nt
	  A = Curl'*Mmu*Curl + 1/dt[i]*Msig
	  if ( (i==1) || (dt[i] != dt[i-1]) )
	    nFacs += 1
	  end
    	  for j = 1:ns
   	 	Gzi = (1/dt[i])*Ne'*getdEdgeMassMatrix(Msh,Ne*(ew[:,j,i+1]-
   	 	       ew[:,j,i]))
   	 	rhs = Gzi*z + 1/dt[i]*Msig*lam[:,j,1]
   	 	lam[:,j,2],EMsolver = solveMaxTime!(A,rhs,Msig,Msh,dt,i,
   	 	                                    nFacs,EMsolver)
	    	
	    	# compute Jv
	    	Jv[:,j,i]  = -P'*(lam[:,j,2])
	    	lam[:,j,1] = lam[:,j,2]
	  end
	end
        Jv = groundedSource ? [vec(Jvdc);vec(Jv)] : vec(Jv)
	return Jv
end

#-------------------------------------------------------

function getSensMatVec(z::Vector{Float64},sigma::Vector{Float64},
                       param::MaxwellTimeBDF2Param)	
	#getSensMatVec(z,sigma,param)
	#This function computes (dData/dsigma)*z for BDF2 time-stepping
	#forward problem with constant step-size. It handles grounded 
	#and inductive sources, assuming DC data is integral of electric  
	#field and not a potential difference (i.e. electric field at t=0  
	#is known) for grounded sources and e0=0 for inductive sources.
	#e_1 is computed by interpolation between two be steps of size
	#2*dt/3. For grounded sources:
	#dCdu= |G'*Msig*G                                                            |
	#      |(-dt*K + Msig)*G Msig                                                |
	#      |1/(2*dt)*Msig*G  -2/dt*Msig      K+3/(2*dt)*Msig                     |
	#      |                 1/(2*dt)*Msig  -2/dt*Msig        K+3/(2*dt)*Msig ...|
	#      |                      .                                              |
	#      |                      .                                              |
	#      |                      .                                              |
	dt       = param.dt
	wave     = param.wave
	Msh      = param.M
	EMsolver = param.EMsolver
	storeEMf = param.storeEMfactors
	DCsolver = param.DCsolver
	storeDCf = param.storeDCfactors
	
	ew       = param.fields
	ehat     = param.ehat
	
	mu     = 4*pi*1e-7  # magnetic permeability
	Curl   = getCurlMatrix(Msh)
	Gin    = getNodalGradientMatrix(Msh)
	Msig   = getEdgeMassMatrix(Msh,vec(sigma))
	Mmu    = getFaceMassMatrix(Msh,vec(zeros(size(sigma)).+1/mu))
        Ne,Qe, = getEdgeConstraints(Msh)
        Nn,Qn  = getNodalConstraints(Msh)
        Nf,Qf  = getFaceConstraints(Msh)
	
	Curl = Qf*Curl*Ne
	Msig = Ne'*Msig*Ne
	Mmu   = Nf'*Mmu*Nf
	G    = Qe*Gin*Nn
	s    = Ne'*copy(param.Sources)
	P    = Ne'*param.Obs

	ns = size(s,2)
	ne = size(Curl,2)
        nt = param.nt
        
	#Check source type
	groundedSource = any( abs(G'*s) .> 1e-12) ? true : false

	#Initialize intermediate and output arrays
	lam = zeros(Float64,ne,ns,3)
	Jv  = zeros(Float64,size(P,2),ns,nt)
	
	#Do the DC part if source is grounded.
	#Note that this code assumes (for grounded sources) that DC data
	#Are computed as the integral of electric field over a receiver
	#dipole and not as a potential difference
	if groundedSource
	  A = G'*Msig*G
	  Jvdc = zeros(Float64,size(P,2),ns)
	  for j = 1:ns
	    Gzi = G'*Ne'*getdEdgeMassMatrix(Msh,-Ne*ew[:,j,1])
	    rhs = Gzi*z
	    lam0,DCsolver    = solveDC!(A,rhs,DCsolver)
	    lam[:,j,1]       = -G*lam0
	    DCsolver.doClear = 0
	    Jvdc[:,j] = -P'*lam[:,j,1]
	  end
	  if ~storeDCf
	    clear!(DCsolver)
	    DCsolver.doClear = 1
	  end
	end
	
	#Do BE and interpolation for first time-step
	A = Curl'*Mmu*Curl + 3/(2*dt)*Msig
	for j = 1:ns
	  Gzi                 = 3/(2*dt)*Ne'*getdEdgeMassMatrix(Msh,Ne*(ehat[:,j]-ew[:,j,1])) 
	  rhs                 = Gzi*z + 3/(2*dt)*Msig*lam[:,j,1]
	  lmTmp,EMsolver      = solveMaxTime!(A,rhs,Msig,Msh,dt,EMsolver)
	  EMsolver.doClear    = 0
	  Gzi                 = Ne'*getdEdgeMassMatrix(Msh,1/dt*Ne*(1.5*ew[:,j,2]-0.75*ehat[:,j]-0.75*ew[:,j,1]))
	  rhs                 = Gzi*z + 3/(4*dt)*Msig*lam[:,j,1] + 3/(4*dt)*Msig*lmTmp
	  lam[:,j,2],EMsolver = solveMaxTime!(A,rhs,Msig,Msh,dt,EMsolver)
	  Jv[:,j,1]           = -P'*lam[:,j,2]
	end
	
	#Do the rest of the time-steps
	for i=2:nt
    	  for j = 1:ns
   	 	Gzi = (1/dt)*Ne'*getdEdgeMassMatrix(Msh,Ne*(1.5*ew[:,j,i+1]-
   	 	       2*ew[:,j,i]+0.5*ew[:,j,i-1]))
   	 	rhs = Gzi*z + 1/dt*Msig*(2*lam[:,j,2]-0.5*lam[:,j,1])
   	 	lam[:,j,3],EMsolver = solveMaxTime!(A,rhs,Msig,Msh,dt,EMsolver)
	    	# compute Jv
	    	Jv[:,j,i]  = -P'*(lam[:,j,3])
	    	lam[:,j,1] = lam[:,j,2]
	    	lam[:,j,2] = lam[:,j,3]
	  end
	end
	if ~storeEMf
          clear!(EMsolver)
          EMsolver.doClear = 1
        end
        Jv = groundedSource ? [vec(Jvdc);vec(Jv)] : vec(Jv)
	return Jv
end

#-------------------------------------------------------

function getSensMatVec(z::Vector{Float64},sigma::Vector{Float64},
                       param::MaxwellTimeSEParam)	
	if isempty(param.Sens)
		tempParam = getMaxwellTimeParam(param.M,param.Sources,param.Obs, param.dt, param.wave,
	                                param.EMsolver,param.DCsolver,param.storeDCfactors)
		# compute sensitivity matrix
		nt = length(param.dt)
		ns = size(param.Sources,2)
		nr = size(param.Obs,2)
		
		J = zeros(nt*ns*nr, length(sigma))
		for k=1:size(J,1)
			v      = zeros(nt*ns*nr)
			v[k]   = 1.0
			JTv    = getSensTMatVec(v,sigma,tempParam)
			J[k,:] = vec(JTv)
		end
		param.Sens = J
		param.fields=[]
	end	
	return param.Sens*z
end

