export getSensTMatVec

function getSensTMatVec(z::Array{Float64,1},sigma::Array{Float64,1},
                        param::MaxwellTimeParam)
	#getSensTMatVec(z,sigma,param)
	#This function computes (dData/dsigma)^T*z for BE time-stepping
	#forward problem. It handles grounded and inductive sources, 
	#assuming DC data is integral of electric field and not a 
	#potential difference (i.e. electric field at t=0 is known) for  
	#grounded sources and e0=0 for inductive sources.
	# magnetic permeability
	mu = 4*pi*1e-7 
	# extract param quantities
	dt       = param.dt
	wave     = param.wave
	Msh      = param.M
	EMsolver = param.EMsolver
	DCsolver = param.DCsolver
	storeDCf = param.storeDCfactors
	ew       = param.fields
	
	
	Curl  = getCurlMatrix(Msh)
	Gin   = getNodalGradientMatrix(Msh)
	Msig  = getEdgeMassMatrix(Msh,vec(sigma))
	Mmu   = getFaceMassMatrix(Msh,vec(zeros(size(sigma)).+1/mu))
        Ne,Qe = getEdgeConstraints(Msh)
        Nn,Qn = getNodalConstraints(Msh)
        Nf,Qf = getFaceConstraints(Msh)
	
	Curl = Qf*Curl*Ne
	Msig = Ne'*Msig*Ne
	Mmu  = Nf'*Mmu*Nf
	G    = Qe*Gin*Nn
	s     = Ne'*param.Sources
        P     = Ne'*param.Obs

	nt = length(dt)
	ns = size(s,2)
	ne = size(Curl,2)
	nr = size(P,2)

	lam = zeros(Float64,ne,ns,2)
	
	#Check source type
	groundedSource = any( abs(G'*s) .> 1e-12) ? true : false
	
	if groundedSource
	  z  = reshape(z,nr,ns,nt+1)
	  pz = zeros(Float64,ne,ns,nt)
	  for i=1:ns
	    for j=1:nt
	      pz[:,i,j] = P*z[:,i,j+1]
	    end
	  end
	else
          z  = reshape(z,nr,ns,nt)
	  pz = zeros(Float64,ne,ns,nt)
	  for i=1:ns
	    for j=1:nt
	      pz[:,i,j] = P*z[:,i,j]
	    end
	  end
	end
	
	JTv = 0
	dt = [dt[:];dt[end]]
	dtLast = 0.0 #Size of last time step, used to check if 
	             #factorization of Forward mod. matrix is needed.
	nFacs = 0
	for i=length(dt)-1:-1:1
          A = Curl'*Mmu*Curl + 1/dt[i]*Msig
          if dt[i] == dtLast
            nFacs += 1
          end
   	  for j = 1:ns
	    rhs = pz[:,j,i] + 1/dt[i+1]*Msig*lam[:,j,2]
	    if norm(rhs) > 1e-20
	      lam[:,j,1],Emsolver = solveMaxTime!(A,rhs,Msig,Msh,dt,i,
	                                          nFacs,EMsolver)
	    end	
	    Gzi = (1/dt[i])*getdEdgeMassMatrix(Msh,Ne*(ew[:,j,i+1]-ew[:,j,i]))
 	    JTv = JTv - Gzi'*(Ne*lam[:,j,1])
 	    lam[:,j,2] = lam[:,j,1]
	  end
	  dtLast = dt[i]
	end

	#Do the DC part if source is grounded.
	#Note that this code assumes (for grounded sources) that DC data
	#Are computed as the integral of electric field over a receiver
	#dipole and not as a potential difference
	if groundedSource
	  A = G'*Msig*G
	  for j = 1:ns
	    pz0 = -G'*P*z[:,j,1]
	    rhs = pz0 - 1/dt[1]*G'*Msig*lam[:,j,2]
	    if norm(rhs) > 1e-20
	      lam0,DCsolver = solveDC!(A,rhs,DCsolver)
	    end
	    Gzi = G'*Ne'*getdEdgeMassMatrix(Msh,-Ne*ew[:,j,1])
	    JTv = JTv -Gzi'*lam0
	  end
	  if ~storeDCf
	    clear!(DCsolver)
	    DCsolver.doClear = 1
	  end
	end
	return JTv
end

#----------------------------------------------------------------------

function getSensTMatVec(z::Array{Float64,1},sigma::Array{Float64,1},
                        param::MaxwellTimeBDF2Param)
	#getSensTMatVec(z,sigma,param)
	#This function computes (dData/dsigma)^T*z for BDF2 time-stepping
	#forward problem. It handles grounded and inductive sources, 
	#assuming DC data is integral of electric field and not a 
	#potential difference (i.e. electric field at t=0 is known) for  
	#grounded sources and e0=0 for inductive sources. FE is used
	#for first time step.
	
	# magnetic permeability
	mu = 4*pi*1e-7 
	# extract param quantities
	dt       = param.dt
	wave     = param.wave
	Msh      = param.M
	EMsolver = param.EMsolver
	storeEMf = param.storeEMfactors
	DCsolver = param.DCsolver
	storeDCf = param.storeDCfactors
	
	ew       = param.fields
	ehat     = param.ehat
	
	Curl  = getCurlMatrix(Msh)
	Gin   = getNodalGradientMatrix(Msh)
	Msig  = getEdgeMassMatrix(Msh,vec(sigma))
	Mmu   = getFaceMassMatrix(Msh,vec(zeros(size(sigma)).+1/mu))
        Ne,Qe = getEdgeConstraints(Msh)
        Nn,Qn = getNodalConstraints(Msh)
        Nf,Qf = getFaceConstraints(Msh)
	
	Curl = Qf*Curl*Ne
	Msig = Ne'*Msig*Ne
	Mmu  = Nf'*Mmu*Nf
	G    = Qe*Gin*Nn
	s     = Ne'*param.Sources
        P     = Ne'*param.Obs

	nt = param.nt
	ns = size(s,2)
	ne = size(Curl,2)
	nr = size(P,2)

	lam = zeros(Float64,ne,ns,3)
	
	#Check source type
	groundedSource = any( abs(G'*s) .> 1e-12) ? true : false
	
	if groundedSource
	  z  = reshape(z,nr,ns,nt+1)
	  pz = zeros(Float64,ne,ns,nt)
	  for i=1:ns
	    for j=1:nt
	      pz[:,i,j] = P*z[:,i,j+1]
	    end
	  end
	else
	  z  = reshape(z,nr,ns,nt)
	  pz = zeros(Float64,ne,ns,nt)
	  for i=1:ns
	    for j=1:nt
	      pz[:,i,j] = P*z[:,i,j]
	    end
	  end
	end
	
	JTv = zeros(length(sigma))
	A = Curl'*Mmu*Curl + 3/(2*dt)*Msig
	for i=nt:-1:2
   	  for j = 1:ns
	    rhs = pz[:,j,i] + 1/dt*Msig*(2*lam[:,j,2]-0.5*lam[:,j,3])
	    if norm(rhs) > 1e-20
	      lam[:,j,1],EMsolver = solveMaxTime!(A,rhs,Msig,Msh,dt,EMsolver)
	      EMsolver.doClear = 0
	    end	
	    Gzi = (1/dt)*getdEdgeMassMatrix(Msh,Ne*(1.5*ew[:,j,i+1]-2*ew[:,j,i]+
	                                    0.5*ew[:,j,i-1])) 
 	    JTv = JTv - Gzi'*(Ne*lam[:,j,1])
 	    lam[:,j,3] = lam[:,j,2]
 	    lam[:,j,2] = lam[:,j,1]
	  end
	end
# 	clear!(mySolver)
	
	#Do first time step stuff.
	for j = 1:ns
	  rhs                 = pz[:,j,1] + 1/dt*Msig*(2*lam[:,j,2]-0.5*lam[:,j,3])
	  lam[:,j,1],EMsolver = solveMaxTime!(A,rhs,Msig,Msh,dt,EMsolver)
	  Gzi                 = getdEdgeMassMatrix(Msh,1/dt*Ne*(1.5*ew[:,j,2]-0.75*ehat[:,j]-0.75*ew[:,j,1]))
	  JTv                 = JTv - Gzi'*(Ne*lam[:,j,1])
	  lam[:,j,3]          = lam[:,j,2]
	  lam[:,j,2]          = lam[:,j,1]
	end
	
	#Do ehat stuff
	lmTmp = zeros(size(lam,1),ns)
	for j = 1:ns
	  rhs = 3/(4*dt)*Msig*lam[:,j,2]
	  lmTmp[:,j],EMsolver = solveMaxTime!(A,rhs,Msig,Msh,dt,EMsolver)
	  Gzi                 = getdEdgeMassMatrix(Msh,3/(2*dt)*Ne*(ehat[:,j]-ew[:,j,1]))
	  JTv                 = JTv - Gzi'*Ne*lmTmp[:,j]
	end
	if ~storeEMf
	  clear!(EMsolver)
	  EMsolver.doClear = 1
	end
	
	#Do the DC part if source is grounded. 
	#Note that this code assumes (for grounded sources) that DC data
	#Are computed as the integral of electric field over a receiver
	#dipole and not as a potential difference
	if groundedSource
	  A0 = G'*Msig*G
	  for j = 1:ns
	    rhs              = -G'*P*z[:,j,1] +
	                       1/(2*dt)*G'*Msig*lam[:,j,3] -
	                       3/(4*dt)*G'*Msig*lam[:,j,2] - 
	                       3/(2*dt)*G'*Msig*lmTmp[:,j]
	    if norm(rhs) > 1e-20
	      lam0,DCsolver = solveDC!(A0,rhs,DCsolver)
	    end
	    Gzi = getdEdgeMassMatrix(Msh,-Ne*ew[:,j,1])
	    JTv = JTv - Gzi'*(Ne*G*lam0)
	  end
	  if ~storeDCf
	    clear!(DCsolver)
	    DCsolver.doClear = 1
	  end
	end
	return JTv
end

#----------------------------------------------------------------------

function getSensTMatVec(z::Vector{Float64},sigma::Vector{Float64},
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
	
	return param.Sens'*z
end