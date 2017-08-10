export getSensTMatVec,interpLocalToGlobal

# For backward compatibility
function getSensTMatVec{T<:Real}(z::Vector{T},sigma::Vector{T},
                                 param::MaxwellTimeParam)
    m   = MaxwellTimeModel(sigma,fill(convert(eltype(sigma),4*pi*1e-7),
                           param.Mesh.nc))
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
            J  = zeros(size(param.ObsTimes,1), length(mod))
            sensTFunc = sensitivityTFunctions[param.timeIntegrationMethod]
            for k=1:size(J,1)
            	v        = zeros(size(param.ObsTimes,1))
            	v[k]     = 1.0
            	JTvStruc = sensTFunc(v,model,param)
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
        #jtz,tmp = sensTFunc(z,model,param)
        return sensTFunc(z,model,param) #,tmp
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
    Mesh          = param.Mesh
    storageLevel  = param.storageLevel
    EMsolvers     = param.EMsolvers
    dt            = param.dt
    ew            = param.fields
    s             = param.Sources

    # Get matrices
    Ne,Qe, = getEdgeConstraints(Mesh)
    Msig   = getEdgeMassMatrix(Mesh,vec(sigma))
    Msig   = Ne'*Msig*Ne
    P      = Ne'*param.Obs
    if invertMu || (param.storageLevel == :None)
        Nf,Qf,    = getFaceConstraints(Mesh)
        Curl      = getCurlMatrix(Mesh)
        Curl      = Qf*Curl*Ne
        Mmu       = getFaceMassMatrix(Mesh,1./mu)
        Mmu       = Nf'*Mmu*Nf
        DmuinvDmu = spdiagm(-1./(mu.^2))
        curllam     = zeros(T,size(Mesh.Qf,1))
        nfcurllam   = zeros(T,size(Mesh.Nf,1))
        Gzitnfcurllam = zeros(T,size(DmuinvDmu,2))
    end
    if param.storageLevel == :None
        K = getMaxwellCurlCurlMatrix!(param,model)
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
    iSolver     = 0
    dtLast      = -1.0 # Size of last time step, used to check if
                      # factorization of Forward mod. matrix is needed.
    Nelam       = zeros(T,size(Ne,1))
    GzitNelam   = zeros(T,Mesh.nc)
    for i=length(dt)-1:-1:1
        if dt[i] != dtLast
            A,iSolver = getBEMatrix!(dt[i],A,K,Msig,param,uniqueSteps)
        end

        rhs = pz[:,:,i] + 1/dt[i+1]*Msig*lam[:,:,2]
        lam[:,:,1],EMsolvers[iSolver] = solveMaxTimeBE!(A,rhs,Msig,Mesh,dt,i,storageLevel,
                                                         EMsolvers[iSolver])

        for j = 1:ns
            if invertSigma
                A_mul_B!(Nelam,Ne,lam[:,j,1]) # Nelam = Ne*lam[:,j,1]
                JTvSigma   .-= (1/dt[i])*dEdgeMassMatrixTrTimesVector(Mesh,sigma,Ne*(ew[:,j,i+1]-ew[:,j,i]),Nelam)
            end
            if invertMu
                curle      = Curl*ew[:,j,i+1]
                Gzi        = getdFaceMassMatrix(Mesh,mu,Nf*curle)*DmuinvDmu
                A_mul_B!(curllam,Curl,lam[:,j,1]) #curllam = Curl*lam[:,j,1]
                A_mul_B!(nfcurllam,Nf,curllam) #nfcurllam  = Nf*curllam
                At_mul_B!(Gzitnfcurllam,Gzi,nfcurllam)
                JTvMu      .-= Gzitnfcurllam
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
        Nn,Qn,   = getNodalConstraints(Mesh)
        Gin      = getNodalGradientMatrix(Mesh)
        G        = Qe*Gin*Nn
        A        = getDCmatrix(Msig,G,param) # Forms matrix only if needed
        pz0      = -G'*P*ptz[:,:,1]
        rhs      = pz0 - 1/dt[1]*G'*Msig*lam[:,:,2]
        lam0,DCsolver = solveDC!(A,rhs,DCsolver)
        for j = 1:ns
            Gzi      = G'*Ne'*getdEdgeMassMatrix(Mesh,-Ne*ew[:,j,1])
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

function getSensTMatVecBDF2{T<:Real}(z::Vector{T},model::MaxwellTimeModel,
    param::MaxwellTimeParam)
    # Unpack model into conductivity and magnetic permeability
    sigma       = model.sigma
    mu          = model.mu
    invertSigma = model.invertSigma
    invertMu    = model.invertMu

    if invertMu
        error("Inverting for mu with bdf2 time-stepping not yet supported")
    end

    #Unpack param
    Mesh          = param.Mesh
    storageLevel  = param.storageLevel
    EMsolvers     = param.EMsolvers
    dt            = param.dt
    ew            = param.fields
    s             = param.Sources
    ehat          = param.AuxFields

    # Get matrices
    Ne,Qe, = getEdgeConstraints(Mesh)
    Msig   = getEdgeMassMatrix(Mesh,vec(sigma))
    Msig   = Ne'*Msig*Ne
    P      = Ne'*param.Obs
    if invertMu || (param.storageLevel == :None)
        Nf,Qf, = getFaceConstraints(Mesh)
        Curl   = getCurlMatrix(Mesh)
        Curl   = Qf*Curl*Ne
        Mmu    = getFaceMassMatrix(Mesh,1./mu)
        Mmu    = Nf'*Mmu*Nf
    else
        Curl = spzeros(T,0,0)
        Mmu  = spzeros(T,0,0)
    end
    K = getMaxwellCurlCurlMatrix!(param,model)

    #Initialize intermediate and output arrays
    ns  = size(s,2)
    ne  = size(ew,1)
    nt  = length(dt)
    nr  = size(P,2)
    lam = zeros(T,ne,ns,3)

    #For debugging
    Nn,Qn,   = getNodalConstraints(Mesh)
    Gin      = getNodalGradientMatrix(Mesh)
    G        = Qe*Gin*Nn
    nn = size(G,2)
    #tmp = zeros(T,nn+5*ne)

    # Multiply by transpose of time interpolation matrix
    # To map from data space to data at all times space
    ptz = param.ObsTimes'*z

    # Check source type
    groundedSource = param.sourceType == :Galvanic ? true : false

    if groundedSource
        ptz  = reshape(ptz,nr,ns,nt+1)
        pz   = zeros(T,ne,ns,nt)
        for i=1:ns
            for j=1:nt
                pz[:,i,j] = P*ptz[:,i,j+1]
            end
        end
    else
        ptz  = reshape(ptz,nr,ns,nt)
        pz   = zeros(T,ne,ns,nt)
        for i=1:ns
            for j=1:nt
                pz[:,i,j] = P*ptz[:,i,j]
            end
        end
    end

#      JTvSigma = 0
#      JTvMu    = 0
    JTv = zeros(T,Mesh.nc)
    uniqueSteps = unique(dt)
    A           = spzeros(T,0,0)
    iSolver     = 0
    #A,iSolver   = getBDF2ConstDTmatrix!(dt[end],A,K,Msig,param,uniqueSteps)
    dt          = [dt;dt[end];dt[end]]
    nelam       = zeros(T,size(Ne,1))
    Gzitnelam   = zeros(T,Mesh.nc)
    for i=nt:-1:1
        tau1 = i > 1 ? dt[i]/dt[i-1] : one(T)
        tau2 = dt[i+1]/dt[i]
        tau3 = dt[i+2]/dt[i+1]
        g1   = (1+2*tau1)/(1+tau1)
        g2   = 1 + tau1
        g3   = (tau1^2)/(1+tau1)
        g2p  = 1 + tau2
        g3p  = (tau3^2)/(1+tau3)
        if (dt[i] != dt[i+1]) || (i==nt)
            A,iSolver = getBDF2ConstDTmatrix!(dt[i],A,K,Msig,param,uniqueSteps)
        end
        if (i>1) && (dt[i] != dt[i-1])
            Atr  = K + (g1/dt[i])*Msig
            if EMsolvers[iSolver].doClear == 1
                clear!(EMsolvers[iSolver])
                EMsolvers[iSolver].Ainv = hasMUMPS ? factorMUMPS(A,1) : cholfact(A)
                EMsolvers[iSolver].doClear = 0
            end
            M = Y -> begin
                         Y = hasMUMPS ? applyMUMPS(EMsolvers[iSolver].Ainv,Y) : EMsolvers[iSolver].Ainv\Y
                         return Y
                     end
        end
        for j = 1:ns
            rhs  = pz[:,j,i] + Msig*(g2p*lam[:,j,2]/dt[i+1]-g3p*lam[:,j,3]/dt[i+2])
            #println("At step $i, norms are $(norm(pz[:,j,i])), $(norm(rhs))")
            if norm(rhs) > 1e-20
                #println("hit norm if at step $i")
                if (i>1) && (dt[i] != dt[i-1])
                    #println("hit cg if at step $i, seeing rhs norm $(norm(rhs))")
                    lam[:,j,1],cgFlag,err,iterTmp, = cg(Atr,rhs,
                       x=vec(lam[:,j,1]),M=M,maxIter=20,tol=param.cgTol)
                    #println("On iter $i $cgFlag $err $(norm(lam[:,j,1])) $(norm(solveMUMPS(Atr,rhs,1)))")
                    if cgFlag != 0
                        warn("getSensTMatVec: cg failed to converge at time step $i. Exited with flag $cgFlag. Reached residual $err with tolerance $(param.cgTol)")
                    end
                else
                    lam[:,j,1],EMsolvers[iSolver] = solveMaxTimeBDF2ConstDT!(A,rhs,Msig,Mesh,dt[i],EMsolvers[iSolver])
                    EMsolvers[iSolver].doClear = 0
                end
            end
            if i>1
                Gzi = getdEdgeMassMatrix(Mesh,1/dt[i]*Ne*(g1*ew[:,j,i+1]-
                              g2*ew[:,j,i] + g3*ew[:,j,i-1]))
            else
                Gzi = getdEdgeMassMatrix(Mesh,1/dt[1]*Ne*(1.5*ew[:,j,2]-
                                              0.75*ehat[:,j]-0.75*ew[:,j,1]))
            end
            A_mul_B!(nelam,Ne,lam[:,j,1]) # nelam = Ne*lam[:,j,1]
            At_mul_B!(Gzitnelam,Gzi,nelam) # Gzi'*nelam
            JTv .-= Gzitnelam
            lam[:,j,3] = lam[:,j,2]
            lam[:,j,2] = lam[:,j,1]
            #tmp[nn+i*ne+1:nn+(i+1)*ne] = lam[:,j,1]
        end
    end

    #Do ehat stuff
    lmTmp = zeros(size(lam,1),ns)
    for j = 1:ns
        rhs = 3/(4*dt[1])*Msig*lam[:,j,2]
        lmTmp[:,j],EMsolvers[iSolver] = solveMaxTimeBDF2ConstDT!(A,rhs,Msig,
                                                  Mesh,dt[1],EMsolvers[iSolver])
        Gzi = getdEdgeMassMatrix(Mesh,3/(2*dt[1])*Ne*(ehat[:,j]-ew[:,j,1]))
        A_mul_B!(nelam,Ne,lmTmp[:,j])
        At_mul_B!(Gzitnelam,Gzi,nelam) # Gzi'*nelam
        JTv .-= Gzitnelam
        #tmp[nn+1:nn+ne] = lam[:,j,1]
    end
    if storageLevel != :Factors
        for solver in EMsolvers
            clear!(solver)
            solver.doClear = 1
        end
    end

    #Do the DC part if source is grounded.
    #Note that this code assumes (for grounded sources) that DC data
    #Are computed as the integral of electric field over a receiver
    #dipole and not as a potential difference
    if groundedSource
        DCsolver = param.DCsolver
        Nn,Qn,   = getNodalConstraints(Mesh)
        Gin      = getNodalGradientMatrix(Mesh)
        G        = Qe*Gin*Nn
        A        = getDCmatrix(Msig,G,param) # Forms matrix only if needed
        tau3     = dt[2]/dt[1]
        g3       = (tau3^2)/(1+tau3)
        for j = 1:ns
            rhs           = -G'*P*ptz[:,j,1] +
                             g3/dt[2]*G'*Msig*lam[:,j,3] -
                             3/(4*dt[1])*G'*Msig*lam[:,j,2] -
                             3/(2*dt[1])*G'*Msig*lmTmp[:,j]
            lam0,DCsolver = solveDC!(A,rhs,DCsolver)
            Gzi           = getdEdgeMassMatrix(Mesh,-Ne*ew[:,j,1])
            JTv           .-= Gzi'*(Ne*G*lam0)
            #tmp[1:nn] = lam0
        end
        if param.storageLevel != :Factors
            clear!(DCsolver)
            DCsolver.doClear = 1
        end
    end
    return MaxwellTimeModel(JTv,[zero(T)],invertSigma,invertMu) #,tmp
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
    M             = param.Mesh
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
    K = getMaxwellCurlCurlMatrix!(param,model)

    #Initialize intermediate and output arrays
    ns  = size(s,2)
    ne  = size(ew,1)
    nt  = length(param.dt)
    nr  = size(P,2)
    lam = zeros(T,ne,ns,3)

    # Multiply by transpose of time interpolation matrix
    # To map from data space to data at all times space
    ptz = param.ObsTimes'*z

    # Check source type
    groundedSource = param.sourceType == :Galvanic ? true : false

    if groundedSource
      ptz  = reshape(ptz,nr,ns,nt+1)
      pz   = zeros(T,ne,ns,nt)
      for i=1:ns
        for j=1:nt
          pz[:,i,j] = P*ptz[:,i,j+1]
        end
      end
    else
      ptz  = reshape(ptz,nr,ns,nt)
      pz   = zeros(T,ne,ns,nt)
      for i=1:ns
        for j=1:nt
          pz[:,i,j] = P*ptz[:,i,j]
        end
      end
    end

#     JTvSigma = 0
#     JTvMu    = 0
    JTv = 0
    uniqueSteps = [dt]
    A           = spzeros(T,0,0)
    A,iSolver   = getBDF2ConstDTmatrix!(dt,A,K,Msig,param,uniqueSteps)
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
            rhs           = -G'*P*ptz[:,j,1] +
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

function getSensTMatVecTRBDF2{T<:Real}(z::Vector{T},model::MaxwellTimeModel,
                                   param::MaxwellTimeParam)
    error("Sensitivities for TRBDF2 not implemented")
end

#----------------------------------------------------------------------

sensitivityTFunctions = Dict(zip(supportedIntegrationMethods,
                                [getSensTMatVecBE;getSensTMatVecBDF2;
                                 getSensTMatVecBDF2ConstDT;
                                 getSensTMatVecTRBDF2]))
