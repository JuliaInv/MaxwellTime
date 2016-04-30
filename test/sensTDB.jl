function sensTDB(z,m,param)

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
	nr = size(param.Obs,2)

	ew  = reshape(param.fields,ne,ns,nt+1);

	lam = zeros(Float64,ne,ns,nt+1)
	
	s = zeros(Float64,ne,ns,nt)
	z = reshape(z,nr,ns,nt)
	for i=1:ns
		for j=1:nt
		  s[:,i,j] = N'*param.Obs*z[:,i,j]
	  end
	end
	
	JTv = 0
	dt = [dt[:];dt[end]]
	for i=length(dt)-1:-1:1
       Ke = Curl'*Mmu*Curl + 1/dt[i]*Msig
   	 for j = 1:ns

			 rhs = s[:,j,i] + 1/dt[i+1]*Msig*lam[:,j,i+1]
			 lam[:,j,i] = Ke\rhs
			 	 
		    Gzi = (1/dt[i])*getdEdgeMassMatrix(param.M,N*(ew[:,j,i+1]-float(i>1)*ew[:,j,i]))
		     
 		    JTv = JTv - Gzi'*(N*lam[:,j,i])
		 end
	 end
	 return JTv
end
