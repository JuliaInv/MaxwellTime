export getSensTMatVec,interpLocalToGlobal

# For backward compatibility
function getSensTMatVec{T<:Real}(z::Vector{T},sigma::Vector{T},
                                 param::MaxwellTimeParam)
    m   = MaxwellTimeModel(sigma,fill(convert(eltype(sigma),4*pi*1e-7),
                           param.M.nc))
    JTz = getSensTMatVec(z,m,param)
    return JTz.sigma
end

function getSensTMatVec{T<:Real}(z::Vector{T},model::MaxwellTimeModel,
                                 param::MaxwellTimeParam)

    if param.sensitivityMethod == :Explicit
        if model.invertSigma && model.invertMu
            error("Using explicit sensitivities while inverting sigma and mu simultaneously is not supported")
        end
        mod = model.invertSigma ? model.sigma : model.mu
        if isempty(param.Sens)
            nt = length(param.dt)
            ns = size(param.Sources,2)
            nr = size(param.Obs,2)
            J  = zeros(nt*ns*nr, length(mod))
            for k=1:size(J,1)
            	v        = zeros(nt*ns*nr)
            	v[k]     = 1.0
            	JTvStruc = getSensTMatVec(v,model,param)
            	JTv      = model.invertSigma ? JTvStruc.sigma : JTvStruc.mu
            	J[k,:] = vec(JTv)
            end
            param.Sens = J
            param.fields=[]
	end	
	JTz = param.Sens'*z
	return MaxwellTimeModel(JTz,JTz,model.invertSigma,model.invertMu)
    else # Implicit sensitivities
        # Dispatch based on forward modelling integration method,
        # e.g. BE vs BDF2
        sensTFunc = sensitivityTFunctions[param.timeIntegrationMethod]
        return sensTFunc(z,model,param)
    end
end
    
function getSensTMatVecBE{T<:Real}(z::Vector{T},model::MaxwellTimeModel,
                                   param::MaxwellTimeParam)
    #getSensTMatVec(z,sigma,param)
    #This function computes (dData/dsigma)^T*z for BE time-stepping
    #forward problem. It handles grounded and inductive sources, 
    #assuming DC data is integral of electric field and not a 
    #potential difference (i.e. electric field at t=0 is known) for  
    #grounded sources and e0=0 for inductive sources.
    # magnetic permeability
	
    # Unpack model into conductivity and magnetic permeability
    sigma       = model.sigma
    mu          = model.mu
    invertSigma = model.invertSigma
    invertMu    = model.invertMu	

    #Unpack param
    M             = param.M
    storageLevel  = param.storageLevel
    EMsolvers     = param.EMsolvers   
    dt            = param.dt
    ew            = param.fields
    s             = param.Sources

    # Get matrices
    Ne,Qe, = getEdgeConstraints(M)
    Msig   = getEdgeMassMatrix(M,vec(sigma))
    Msig   = Ne'*Msig*Ne
    P      = Ne'*param.Obs
    if invertMu || (param.storageLevel == :None)
        Nf,Qf, = getFaceConstraints(M)
        Curl   = getCurlMatrix(M)
        Curl   = Qf*Curl*Ne
        Mmu    = getFaceMassMatrix(M,1./mu)
        Mmu    = Nf'*Mmu*Nf
    end
    if param.storageLevel == :None
        K = Curl'*Mmu*Curl
    else
        K = spzeros(T,0,0)
    end

    #Initialize intermediate and output arrays
    ns  = size(s,2)
    ne  = size(ew,1)
    nt  = length(dt)
    nr  = size(P,2)
    lam = zeros(T,ne,ns,2)
    
    # Multiply by transpose of time interpolation matrix
    # To map from data space to data at all times space
    ptz = param.ObsTimes'*z
    
    # Check source type
    groundedSource = param.sourceType == :Galvanic ? true : false
    if groundedSource
        ptz  = reshape(ptz,nr,ns,nt+1)
        pz   = zeros(T,ne,ns,nt)
        for j=1:nt
            pz[:,:,j] = P*ptz[:,:,j+1]
        end
    else
        ptz  = reshape(ptz,nr,ns,nt)
        pz   = zeros(T,ne,ns,nt)
        for j=1:nt
            pz[:,:,j] = P*ptz[:,:,j]
        end
    end
    
    JTvSigma    = zeros(T,length(sigma))
    JTvMu       = zeros(T,length(mu))
    dt          = [dt[:];dt[end]]
    uniqueSteps = unique(dt)
    A           = spzeros(T,0,0)
    dtLast      = 0.0 # Size of last time step, used to check if 
                      # factorization of Forward mod. matrix is needed.
    for i=length(dt)-1:-1:1
        if dt[i] != dtLast
            A,iSolver = getBEMatrix!(dt[i],A,K,Msig,param,uniqueSteps)
        end
        
        rhs = pz[:,:,i] + 1/dt[i+1]*Msig*lam[:,:,2]
        lam[:,:,1],EMsolvers[iSolver] = solveMaxTimeBE!(A,rhs,Msig,M,dt,i,storageLevel,
                                                         EMsolvers[iSolver])
        for j = 1:ns
            if invertSigma
                Gzi        = (1/dt[i])*getdEdgeMassMatrix(M,Ne*(ew[:,j,i+1]-ew[:,j,i]))
                JTvSigma   = JTvSigma - Gzi'*(Ne*lam[:,j,1])  
            end
            if invertMu
                Gzi        = (1/dt[i])*getdFaceMassMatrix(M,mu,Nf*Curl*ew[:,j,i+1])
                JTvMu      = JTvMu - Gzi'*(Curl*Nf*lam[:,j,1])  
            end
        end
        lam[:,:,2] = lam[:,:,1]
        dtLast = dt[i]
    end
    if storageLevel != :Factors
        for solver in EMsolvers
          clear!(EMsolver)
          solver.doClear = 1
        end
    end 
    
    #Do the DC part if source is grounded.
    #Note that this code assumes (for grounded sources) that DC data
    #Are computed as the integral of electric field over a receiver
    #dipole and not as a potential difference
    if groundedSource && invertSigma
        DCsolver = param.DCsolver
        Nn,Qn,   = getNodalConstraints(M)
        Gin      = getNodalGradientMatrix(M)
        G        = Qe*Gin*Nn
        A        = getDCmatrix(Msig,G,param) # Forms matrix only if needed
        pz0      = -G'*P*ptz[:,:,1]
        rhs      = pz0 - 1/dt[1]*G'*Msig*lam[:,:,2]
        lam0,DCsolver = solveDC!(A,rhs,DCsolver)
        for j = 1:ns
            Gzi      = G'*Ne'*getdEdgeMassMatrix(M,-Ne*ew[:,j,1])
            JTvSigma = JTvSigma -Gzi'*lam0[:,j]
        end
        if param.storageLevel != :Factors
            clear!(DCsolver)
            DCsolver.doClear = 1
        end
    end
    return MaxwellTimeModel(JTvSigma,JTvMu,invertSigma,invertMu)
end

#----------------------------------------------------------------------

function getSensTMatVecBDF2ConstDT{T<:Real}(z::Vector{T},model::MaxwellTimeModel,
                                   param::MaxwellTimeParam)
	#getSensTMatVec(z,sigma,param)
	#This function computes (dData/dsigma)^T*z for BDF2 time-stepping
	#forward problem. It handles grounded and inductive sources, 
	#assuming DC data is integral of electric field and not a 
	#potential difference (i.e. electric field at t=0 is known) for  
	#grounded sources and e0=0 for inductive sources. FE is used
	#for first time step.

    # Unpack model into conductivity and magnetic permeability
    sigma       = model.sigma
    mu          = model.mu
    invertSigma = model.invertSigma
    invertMu    = model.invertMu	
    
    if invertMu
        error("Inverting for mu with bdf2 time-stepping not yet supported")
    end

    #Unpack param
    M             = param.M
    storageLevel  = param.storageLevel
    EMsolver      = param.EMsolvers[1]
    dt            = param.dt[1]
    ew            = param.fields
    s             = param.Sources
    ehat          = param.AuxFields

    # Get matrices
    Ne,Qe, = getEdgeConstraints(M)
    Msig   = getEdgeMassMatrix(M,vec(sigma))
    Msig   = Ne'*Msig*Ne
    P      = Ne'*param.Obs
    if invertMu || (param.storageLevel == :None)
        Nf,Qf, = getFaceConstraints(M)
        Curl   = getCurlMatrix(M)
        Curl   = Qf*Curl*Ne
        Mmu    = getFaceMassMatrix(M,1./mu)
        Mmu    = Nf'*Mmu*Nf
    else
        Curl = spzeros(T,0,0)
        Mmu  = spzeros(T,0,0)
    end

    #Initialize intermediate and output arrays
    ns  = size(s,2)
    ne  = size(ew,1)
    nt  = length(param.dt)
    nr  = size(P,2)
    lam = zeros(T,ne,ns,3)
	
    # Check source type
    groundedSource = param.sourceType == :Galvanic ? true : false
	
    if groundedSource
      z  = reshape(z,nr,ns,nt+1)
      pz = zeros(T,ne,ns,nt)
      for i=1:ns
        for j=1:nt
          pz[:,i,j] = P*z[:,i,j+1]
        end
      end
    else
      z  = reshape(z,nr,ns,nt)
      pz = zeros(T,ne,ns,nt)
      for i=1:ns
        for j=1:nt
          pz[:,i,j] = P*z[:,i,j]
        end
      end
    end
	
#     JTvSigma = 0
#     JTvMu    = 0
    JTv = 0
    A        = getBDF2ConstDTmatrix(Msig,Mmu,Curl,dt,param)
    for i=nt:-1:2
      for j = 1:ns
        rhs = pz[:,j,i] + 1/dt*Msig*(2*lam[:,j,2]-0.5*lam[:,j,3])
        if norm(rhs) > 1e-20
          lam[:,j,1],EMsolver = solveMaxTimeBDF2ConstDT!(A,rhs,Msig,M,dt,EMsolver)
          EMsolver.doClear = 0
        end	
        Gzi = (1/dt)*getdEdgeMassMatrix(M,Ne*(1.5*ew[:,j,i+1]-2*ew[:,j,i]+
                                        0.5*ew[:,j,i-1])) 
        JTv = JTv - Gzi'*(Ne*lam[:,j,1])
        lam[:,j,3] = lam[:,j,2]
        lam[:,j,2] = lam[:,j,1]
      end
    end
    
    #Do first time step stuff.
    for j = 1:ns
      rhs                 = pz[:,j,1] + 1/dt*Msig*(2*lam[:,j,2]-0.5*lam[:,j,3])
      lam[:,j,1],EMsolver = solveMaxTimeBDF2ConstDT!(A,rhs,Msig,M,dt,EMsolver)
      Gzi                 = getdEdgeMassMatrix(M,1/dt*Ne*(1.5*ew[:,j,2]-0.75*ehat[:,j]-0.75*ew[:,j,1]))
      JTv                 = JTv - Gzi'*(Ne*lam[:,j,1])
      lam[:,j,3]          = lam[:,j,2]
      lam[:,j,2]          = lam[:,j,1]
    end
    
    #Do ehat stuff
    lmTmp = zeros(size(lam,1),ns)
    for j = 1:ns
      rhs = 3/(4*dt)*Msig*lam[:,j,2]
      lmTmp[:,j],EMsolver = solveMaxTimeBDF2ConstDT!(A,rhs,Msig,M,dt,EMsolver)
      Gzi                 = getdEdgeMassMatrix(M,3/(2*dt)*Ne*(ehat[:,j]-ew[:,j,1]))
      JTv                 = JTv - Gzi'*Ne*lmTmp[:,j]
    end
    if storageLevel != :Factors
      clear!(EMsolver)
      EMsolver.doClear = 1
    end
    
    #Do the DC part if source is grounded. 
    #Note that this code assumes (for grounded sources) that DC data
    #Are computed as the integral of electric field over a receiver
    #dipole and not as a potential difference
    if groundedSource
        DCsolver = param.DCsolver
        Nn,Qn,   = getNodalConstraints(M)
        Gin      = getNodalGradientMatrix(M)
        G        = Qe*Gin*Nn
        A        = getDCmatrix(Msig,G,param) # Forms matrix only if needed
        for j = 1:ns
            rhs           = -G'*P*z[:,j,1] +
                             1/(2*dt)*G'*Msig*lam[:,j,3] -
                             3/(4*dt)*G'*Msig*lam[:,j,2] - 
                             3/(2*dt)*G'*Msig*lmTmp[:,j]
            lam0,DCsolver = solveDC!(A,rhs,DCsolver)
            Gzi           = getdEdgeMassMatrix(M,-Ne*ew[:,j,1])
            JTv           = JTv - Gzi'*(Ne*G*lam0)
        end
        if param.storageLevel != :Factors
            clear!(DCsolver)
            DCsolver.doClear = 1
        end
    end
    return MaxwellTimeModel(JTv,[zero(T)],invertSigma,invertMu)
end

#-------------------------------------------------------

function getSensTMatVecBDF2{T<:Real}(z::Vector{T},model::MaxwellTimeModel,
                                   param::MaxwellTimeParam)
    error("Sensitivities for variable step-size bdf2 not implemented")
end
#-------------------------------------------------------


function getSensTMatVecTRBDF2{T<:Real}(z::Vector{T},model::MaxwellTimeModel,
                                   param::MaxwellTimeParam)
    error("Sensitivities for TRBDF2 not implemented")
end

#----------------------------------------------------------------------

sensitivityTFunctions = Dict(zip(supportedIntegrationMethods,
                                [getSensTMatVecBE;getSensTMatVecBDF2;
                                 getSensTMatVecBDF2ConstDT;
                                 getSensTMatVecTRBDF2]))