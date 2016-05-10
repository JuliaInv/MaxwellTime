export getData

function getData(sigma::Array{Float64,1},param::MaxwellTimeParam)
#function getData(sigma,param)

# Solve the Backward Euler system
# Curl*e + b_t = 0
# Curl'*(Mmuinv*b) - Msig*e = s(t)
#
# By eliminating b and solving for e
# b_t =  -Curl*e
# (Curl'*Mmuinv*Curl + 1/dt*Msig)*e_{n+1} = - 1/dt*(Msig*e_{n} + s_n - s_{n+1})
#

	dt       = param.dt
	wave     = param.wave
	Msh      = param.M
	mySolver = param.solver
	

	mu    = 4*pi*1e-7	
	Curl  = getCurlMatrix(Msh)
	Msig  = getEdgeMassMatrix(Msh,vec(sigma))
	Mmu   = getFaceMassMatrix(Msh,1/mu*ones(size(sigma)))
        Ne,Qe = getEdgeConstraints(Msh)
        G     = getNodalGradientMatrix(Msh)
        Nn,Qn = getNodalConstraints(Msh)
        Nf,Qf = getFaceConstraints(Msh)
        
        G     = Qe*G*Nn
        Curl  = Qf*Curl*Ne
        Mmu   = Nf'*Mmu*Nf
        Msig  = Ne'*Msig*Ne
        s     = Ne'*copy(param.Sources)
        P     = Ne'*param.Obs
#         G     = Ne'*G*Nn
#         Curl  = Curl*Ne
#         Msig  = Ne'*Msig*Ne
#         s     = Ne'*copy(param.Sources)
#         P     = Ne'*param.Obs
      
	# allocate space for fields
	numSrc = size(s,2)
        #numSrc  = size(b0,2)
	ew = zeros(size(Ne,2),numSrc,length(dt)+1);
	#e0   = Msig\(Curl'*Mmu*b0)
	#ew[:,:,1] = wave[1]*Ne'*e0

	#Compute e0. Check divergence first to determine source type
	if any( abs(G'*s) .> 1e-12)
	  groundedSource = true
	  Adc = G'*Msig*G
	  phi0,mySolver = solveLinearSystem(Adc,G'*s,mySolver,1,0)
	  ew[:,:,1] = -G*phi0
	  clear!(mySolver)
	else
	  groundedSource = false
	end
	
	# time step
	dtLast = 0.0
	A = []
	for i=1:length(dt)
            dtinv = 1.0/dt[i]
            rhs = dtinv*(Msig*ew[:,:,i]+(wave[i]-wave[i+1])*s)
            if (dt[i] != dtLast)
              mySolver.doClear = 1
              A = Curl'*Mmu*Curl + dtinv*Msig
            end
            ew[:,:,i+1],mySolver = solveMaxTime(A,rhs,Msig,Msh,dt[i],mySolver)
	    mySolver.doClear = 0
	    dtLast = dt[i]
	end
	clear!(mySolver)
        mySolver.doClear = 1
	
	# compute the data
	#Note (for grounded sources) that DC data
	#Are computed as the integral of electric field over a receiver
	#dipole and not as a potential difference
	if groundedSource
	  D = zeros(size(P,2),numSrc,length(dt)+1)
	  D[:,:,1] = P'*ew[:,:,1]
	  offset = 1
	else
	  D = zeros(size(P,2),numSrc,length(dt))
	  offset = 0
	end
	for i=1:length(dt)
		D[:,:,i+offset] = P'*ew[:,:,i+1]
	end
	param.fields = ew
	
	return D, param
end

#--------------------------------------------------------------------

function getData(sigma::Array{Float64,1},param::MaxwellTimeBDF2Param)
  #function getData(sigma,param)
  #Use BDF2 backward differentiation formula for second order time-stepping
	dt       = param.dt
	nt       = param.nt
	wave     = param.wave
	Msh      = param.M
	EMsolver = param.EMsolver
	storeEMf = param.storeEMfactors
	DCsolver = param.DCsolver
	storeDCf = param.storeDCfactors
	

	mu    = 4*pi*1e-7	
	Curl  = getCurlMatrix(Msh)
	Msig  = getEdgeMassMatrix(Msh,vec(sigma))
	Mmu   = getFaceMassMatrix(Msh,1/mu*ones(size(sigma)))
        Ne,Qe = getEdgeConstraints(Msh)
        G     = getNodalGradientMatrix(Msh)
        Nn,Qn = getNodalConstraints(Msh)
        Nf,Qf = getFaceConstraints(Msh)
        
        G     = Qe*G*Nn
        Curl  = Qf*Curl*Ne
        Mmu   = Nf'*Mmu*Nf
        Msig  = Ne'*Msig*Ne
        s     = Ne'*copy(param.Sources)
        P     = Ne'*param.Obs
        
        # allocate space for fields
	numSrc = size(s,2)
	ew = zeros(size(Ne,2),numSrc,nt+1);

	#Compute e_0. Check divergence first to determine source type
	if any( abs(G'*s) .> 1e-12)
	  groundedSource   = true
	  DCsolver.doClear = 1
	  Adc              = G'*Msig*G
	  #phi0,mySolver    = solveLinearSystem(Adc,full(G'*s),DCsolver,1,0)
	  phi0,DCSolver    = solveDC!(Adc,G'*s,DCsolver)
	  ew[:,:,1]        = -G*phi0
	  DCsolver.doClear = 0
	  if ~storeDCf
	    clear!(DCsolver)
	    DCsolver.doClear = 1
	  end
	else
	  groundedSource = false
	end

        #Element-wise linear interp.
        EMsolver.doClear = 1
	A    = Curl'*Mmu*Curl + 3/(2*dt)*Msig
	rhs  = 3/(2*dt)*( Msig*ew[:,:,1] + (wave[1]-wave[2])*s )
	ehat,EMsolver = solveMaxTime!(A,rhs,Msig,Msh,2/(3*dt),EMsolver)
	EMsolver.doClear = 0
	rhs  = 3/(2*dt)*( Msig*ehat )
	ehat2,EMsolver = solveMaxTime!(A,rhs,Msig,Msh,2/(3*dt),EMsolver)
	ew[:,:,2] = 0.5*(ehat+ehat2)
	param.ehat = ehat
	
	#Continue time-stepping using BDF2
	for i=2:nt
	  rhs = -3/(2*dt)*( (wave[i+1]-(4/3)*wave[i]+wave[i-1]/3)*s + 
	        Msig*(-4/3*ew[:,:,i] + ew[:,:,i-1]/3 ) )
	  ew[:,:,i+1],EMsolver = solveMaxTime!(A,rhs,Msig,Msh,dt,EMsolver)
	end
	if ~storeEMf
	  clear!(EMsolver)
	  EMsolver.doClear = 1
	end
	
	# compute the data
	#Note (for grounded sources) that DC data
	#Are computed as the integral of electric field over a receiver
	#dipole and not as a potential difference
	if groundedSource
	  D = zeros(size(P,2),numSrc,nt+1)
	  D[:,:,1] = P'*ew[:,:,1]
	  offset = 1
	else
	  D = zeros(size(P,2),numSrc,nt)
	  offset = 0
	end
	for i=1:nt
		D[:,:,i+offset] = P'*ew[:,:,i+1]
	end
	param.fields = ew
	
	return D, param
end

#--------------------------------------------------------------------

function getData(sigma::Array{Float64,1},param::MaxwellTimeTRBDF2Param)
  #function getData(sigma,param)
  #Use BDF2 backward differentiation formula for second order time-stepping
	dt       = param.dt
	nt       = param.nt
	wave     = param.wave
	Msh      = param.M
	mySolver = param.solver
	

	mu    = 4*pi*1e-7	
	Curl  = getCurlMatrix(Msh)
	Msig  = getEdgeMassMatrix(Msh,vec(sigma))
	Mmu   = getFaceMassMatrix(Msh,1/mu*ones(size(sigma)))
        Ne,Qe = getEdgeConstraints(Msh)
        G     = getNodalGradientMatrix(Msh)
        Nn,Qn = getNodalConstraints(Msh)
        Nf,Qf = getFaceConstraints(Msh)
        
        G     = Qe*G*Nn
        Curl  = Qf*Curl*Ne
        Mmu   = Nf'*Mmu*Nf
        Msig  = Ne'*Msig*Ne
        s     = Ne'*copy(param.Sources)
        P     = Ne'*param.Obs
        
        # allocate space for fields
	numSrc = size(s,2)
        #numSrc  = size(b0,2)
	ew = zeros(size(Ne,2),numSrc,nt+1);

	#Compute e_0. Check divergence first to determine source type
	mySolver.doClear = 1
	if any( abs(G'*s) .> 1e-12)
	  groundedSource = true
	  Adc = G'*Msig*G
	  phi0,mySolver = solveLinearSystem(Adc,G'*s,mySolver,1,0)
	  ew[:,:,1] = -G*phi0
	  clear!(mySolver)
	else
	  groundedSource = false
	end


        #Setup gamma related stuff
        gma     = 2-sqrt(2)
        gmaFctr = 2/(gma*dt)
        f1      = 1/(gma*(2-gma))
        f2      = ((1-gma)^2)/(gma*(2-gma))

        #Factor the matrix once and for all
        mySolver.doClear = 1
        K = Curl'*Mmu*Curl
        A = K + gmaFctr*Msig
	
	#Time-stepping
	for i=1:nt
	  
	  #Trapezoidal rule to get e_(n+gma)
	  if i > 1
	    rhs = -K*ew[:,:,i] + gmaFctr*( Msig*ew[:,:,i] + (wave[i]-wave[i+1])*s )
	  else
	    rhs = gmaFctr*( Msig*ew[:,:,i] + (wave[i]-wave[i+1])*s )
	  end
	  eGma,mySolver = solveMaxTime(A,rhs,Msig,Msh,dt,mySolver)
	  mySolver.doClear = 0
	  
	  #BDF-2 to get e_(n+1)
	  wvGma = 0.0
	  rhs = gmaFctr*( f1*Msig*eGma - f2*Msig*ew[:,:,i] + (wvGma-wave[i]-wave[i+1])*s )
	  ew[:,:,i+1],mySolver = solveMaxTime(A,rhs,Msig,Msh,dt,mySolver)
	end
	mySolver.doClear = 1
	
	# compute the data
	#Note (for grounded sources) that DC data
	#Are computed as the integral of electric field over a receiver
	#dipole and not as a potential difference
	if groundedSource
	  D = zeros(size(P,2),numSrc,nt+1)
	  D[:,:,1] = P'*ew[:,:,1]
	  offset = 1
	else
	  D = zeros(size(P,2),numSrc,nt)
	  offset = 0
	end
	for i=1:nt
		D[:,:,i+offset] = P'*ew[:,:,i+1]
	end
	param.fields = ew
	
	return D, param #,ehat
end

function getData(sigma::Array{Float64,1},param::MaxwellTimeSEParam)

#function getTDdata(sigma::Array{Float64,1},M::AbstractMesh,P::SparseMatrixCSC,b0,dt::Vector,wave::Vector)
# d,e = getTDdata(sigma::Vector,M::Mesh,P::SparseMatrixCSC,b0::Vector,dt::Vector,wave::Vector)
# Solve the Backward Euler system
# Curl*e + b_t = 0
# Curl'*(Mmuinv*b) - Msig*e = s(t)
#
# By eliminating b and solving for e
# b_t =  -Curl*e
# (Curl'*Mmuinv*Curl + 1/dt*Msig)*en = - 1/dt*Msig*e_{n-1} + s(n+1)
#
	tempParam = getMaxwellTimeParam(param.M,param.Sources,param.Obs, param.dt, param.wave)
	param.Sens=[]
	D,tempParam = getData(sigma,tempParam)
	param.fields = tempParam.fields
	return D, param
end