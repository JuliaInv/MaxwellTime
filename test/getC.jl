function C(u,sigma::Array{Float64,1},param::MaxwellTimeParam)
	dt       = param.dt[1]
	nt       = length(param.dt)
	wave     = param.wave
	Msh      = param.M

	mu     = 4*pi*1e-7	
	Curl   = getCurlMatrix(Msh)
	Msig   = getEdgeMassMatrix(Msh,vec(sigma))
	Mmu    = getFaceMassMatrix(Msh,1/mu*ones(size(sigma)))
        Ne,Qe, = getEdgeConstraints(Msh)
        G      = getNodalGradientMatrix(Msh)
        Nn,Qn, = getNodalConstraints(Msh)
        Nf,Qf, = getFaceConstraints(Msh)
        
        G     = Qe*G*Nn
        Curl  = Qf*Curl*Ne
        Mmu   = Nf'*Mmu*Nf
        Msig  = Ne'*Msig*Ne
        s     = Ne'*copy(param.Sources)
        
        K     = Curl'*Mmu*Curl
        A     = K + 3/(2*dt)*Msig
        #Ainv  = inv(A)
        Ainv   = A\Msig
#         println(size(Ainv))
#         println(typeof(Ainv))
#         println(size(dt))
#         println(typeof(dt))
        
        nen  = size(Ne,2)
        nnn  = size(Nn,2)
        blnkn = spzeros(nnn,4*nen)
        blnkee = spzeros(nen,nen)
        blnken = spzeros(nen,nnn)
        dCdu  = [G'*Msig*G blnkn;
                (3/(4*dt)*Msig+9/(8*dt^2)*Msig*Ainv)*G A blnkee blnkee blnkee;
                -1/(2*dt)*Msig*G -2/dt*Msig A blnkee blnkee;
                blnken 1/(2*dt)*Msig -2/dt*Msig A blnkee;
                blnken blnkee 1/(2*dt)*Msig -2/dt*Msig A]
         
        q = [G'*vec(s);zeros(4*size(Curl,2))]
        
        return dCdu*u-q
end