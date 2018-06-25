export getSensTMatVec,interpLocalToGlobal

# For backward compatibility
function getSensTMatVec{Tf<:Real}(z::Vector{Tf},sigma::Vector{Tf},
                                 param::MaxwellTimeParam)
    m   = MaxwellTimeModel(sigma,fill(convert(eltype(sigma),4*pi*1e-7),
                           param.Mesh.nc))
    JTz = getSensTMatVec(z,m,param)
    return JTz.values["sigmaCell"]
end

function getSensTMatVec{Tf<:Real}(z::Vector{Tf},model::MaxwellTimeModel,
                                 param::MaxwellTimeParam)

    if param.sensitivityMethod == :Explicit
        invertSigma = in("sigmaCell", model.activeInversionProperties)
        invertMu    = in(   "muCell", model.activeInversionProperties)
        if invertSigma && invertMu
            error("Using explicit sensitivities while inverting sigma and mu simultaneously is not supported")
        end
        mod = invertSigma ? model.values["sigmaCell"] : model.values["muCell"]
        if isempty(param.Sens)
            ndata     = size(param.ObsTimes,1)
            J         = Array{Float64}(ndata, length(mod))
            sensTFunc = sensitivityTFunctions[param.timeIntegrationMethod]
            v = Array{Float64}(ndata)
            for k=1:ndata
            	fill!(v, 0.0)
            	v[k]     = 1.0
            	JTvStruc = sensTFunc(v,model,param)
            	JTv      = invertSigma ? JTvStruc.values["sigmaCell"] : JTvStruc.values["muCell"]
            	J[k,:] = vec(JTv)
            end
            param.Sens = J
            param.fields=[]

            # Clear all factorizations.
            EMsolvers = param.EMsolvers
            DCsolver = param.DCsolver
            for solver in EMsolvers
              clear!(solver)
              solver.doClear = 1
            end
            clear!(DCsolver)
            DCsolver.doClear = 1
            gc()

        end
	     JTz = param.Sens'*z
         key = invertSigma ? "sigmaCell" : "muCell"
	     return MaxwellTimeModel(Dict(key=>JTz),model.activeInversionProperties)
    else # Implicit sensitivities
        # Dispatch based on forward modelling integration method,
        # e.g. BE vs BDF2
        sensTFunc = sensitivityTFunctions[param.timeIntegrationMethod]
        #jtz,tmp = sensTFunc(z,model,param)
        return sensTFunc(z,model,param) #,tmp
    end
end

function getSensTMatVecBE{Tf<:Real}(z::Vector{Tf},model::MaxwellTimeModel,
                                   param::MaxwellTimeParam)
    #getSensTMatVec(z,sigma,param)
    #This function computes (dData/dsigma)^T*z for BE time-stepping
    #forward problem. It handles grounded and inductive sources,
    #assuming DC data is integral of electric field and not a
    #potential difference (i.e. electric field at t=0 is known) for
    #grounded sources and e0=0 for inductive sources.
    # magnetic permeability

    # Unpack model into conductivity and magnetic permeability
    sigma = model.values["sigmaCell"]
    sigma = param.modUnits == :res ? 1./sigma : sigma
    mu          = model.values["muCell"]
    invertSigma = in("sigmaCell", model.activeInversionProperties)
    invertMu    = in(   "muCell", model.activeInversionProperties)

    #Unpack param
    Mesh          = param.Mesh
    storageLevel  = param.storageLevel
    EMsolvers     = param.EMsolvers
    dt            = param.dt
    ew            = param.fields
    s             = param.Sources

    # Get matrices
    Tn     = typeof(param.Mesh.nc)
    Ne,Qe, = getEdgeConstraints(Mesh)
    Msig   = getEdgeMassMatrix(Mesh,vec(sigma))
    Msig   = Ne'*Msig*Ne
    P      = param.Obs
    DrhoDsig = param.modUnits == :res ? spdiagm(-(sigma.^2)) : UniformScaling(1.0)
    if invertMu || (param.storageLevel == :None)
        Nf,Qf,    = getFaceConstraints(Mesh)
        Curl      = getCurlMatrix(Mesh)
        Curl      = Qf*Curl*Ne
        Mmu       = getFaceMassMatrix(Mesh,1./mu)
        Mmu       = Nf'*Mmu*Nf
        DmuinvDmu = spdiagm(-1./(mu.^2))
        curllam     = zeros(Tf,size(Mesh.Qf,1))
        nfcurllam   = zeros(Tf,size(Mesh.Nf,1))
        Gzitnfcurllam = zeros(Tf,size(DmuinvDmu,2))
    end

    #Initialize intermediate and output arrays
    ns  = size(s,2)
    ne  = size(ew,1)
    nt  = length(dt)
    nr  = size(P,2)
    lam = zeros(Tf,ne,ns) # only one lam is needed

    # Multiply by transpose of time interpolation matrix
    # To map from data space to data at all times space
    ptz = param.ObsTimes'*z

    # Check source type
    groundedSource = param.sourceType == :Galvanic ? true : false
    if groundedSource
        ptz  = reshape(ptz,nr,ns,nt+1)
#        pz   = zeros(T,ne,ns,nt)
#        for j=1:nt
#            pz[:,:,j] = P*ptz[:,:,j+1]
#        end
    else
        ptz  = reshape(ptz,nr,ns,nt)
#        pz   = zeros(T,ne,ns,nt)
#        for j=1:nt
#            pz[:,:,j] = P*ptz[:,:,j]
#        end
    end

    JTv   = Dict{String,Vector{Float64}}()
    if invertSigma
        JTv["sigmaCell"] = zeros(Tf,length(sigma))
    end
    if invertMu
        JTv["muCell"] = zeros(Tf,length(mu))
    end
    dt          = [dt[:];dt[end]]
    uniqueSteps = unique(dt)
    A           = speye(Tf,Tn,size(Ne,2)) #spzeros(T,0,0)
    iSolver     = 0
    dtLast      = -1.0 # Size of last time step, used to check if
                      # factorization of Forward mod. matrix is needed.
    Nelam       = Array{Tf}(size(Ne,1))
   # GzitNelam   = zeros(T,Mesh.nc)
    dfields = Array{Float64}(size(Ne,1))

    for i=length(dt)-1:-1:1
        if dt[i] != dtLast
            A,iSolver = getBEMatrix!(dt[i],model,Msig,param,uniqueSteps)
        end

        ij = groundedSource ? i+1 : i
       # rhs = pz[:,:,i] + 1/dt[i+1]*Msig*lam[:,:,2]
        rhs = P*ptz[:,:,ij] + (1/dt[i+1])*Msig*lam


        lam, EMsolvers[iSolver] = solveMaxTimeBE!(A,rhs,Msig,Mesh,dt,i,storageLevel,
                                                  EMsolvers[iSolver])

        for j = 1:ns
            if invertSigma
                A_mul_B!(Nelam,Ne,lam[:,j]) # Nelam = Ne*lam[:,j]
                A_mul_B!(dfields, Ne, ew[:,j,i+1]-ew[:,j,i])  # dfields = Ne*(ew[:,j,i+1]-ew[:,j,i])
                JTv["sigmaCell"] .-= (1/dt[i])*(DrhoDsig'*dEdgeMassMatrixTrTimesVector(Mesh,sigma, dfields, Nelam))
            end
            if invertMu
                curle      = Curl*ew[:,j,i+1]
                Gzi        = getdFaceMassMatrix(Mesh,mu,Nf*curle)*DmuinvDmu
                A_mul_B!(curllam,Curl,lam[:,j]) #curllam = Curl*lam[:,j]
                A_mul_B!(nfcurllam,Nf,curllam) #nfcurllam  = Nf*curllam
                At_mul_B!(Gzitnfcurllam,Gzi,nfcurllam)
                JTv["muCell"]      .-= Gzitnfcurllam
            end
        end
       # lam[:,:,2] = lam[:,:,1]
        dtLast = dt[i]
    end  # i
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
    if groundedSource && invertSigma
        DCsolver = param.DCsolver
        Nn,Qn,   = getNodalConstraints(Mesh)
        Gin      = getNodalGradientMatrix(Mesh)
        G        = Qe*Gin*Nn
        A        = getDCmatrix(Msig,G,param) # Forms matrix only if needed
        pz0      = -G'*P*ptz[:,:,1]
        rhs      = pz0 - (1/dt[1])*G'*(Msig*lam)
        lam0,DCsolver = solveDC!(A,rhs,DCsolver)
        for j = 1:ns
            Gzi      = G'*DrhoDsig*Ne'*getdEdgeMassMatrix(Mesh,sigma,-Ne*ew[:,j,1])
            JTv["sigmaCell"] .-= Gzi'*lam0[:,j]
        end
        if param.storageLevel != :Factors
            clear!(DCsolver)
            DCsolver.doClear = 1
        end
    end
    return MaxwellTimeModel(JTv, model.activeInversionProperties)
end

#----------------------------------------------------------------------

function getSensTMatVecBDF2{Tf<:Real}(z::Vector{Tf},model::MaxwellTimeModel,
    param::MaxwellTimeParam)
    # Unpack model into conductivity and magnetic permeability
    sigma       = model.values["sigmaCell"]
    mu          = model.values["muCell"]
    invertSigma = in("sigmaCell", model.activeInversionProperties)
    invertMu    = in(   "muCell", model.activeInversionProperties)

    #Unpack param
    Mesh          = param.Mesh
    storageLevel  = param.storageLevel
    EMsolvers     = param.EMsolvers
    dt            = param.dt
    ew            = param.fields
    s             = param.Sources
    ehat          = param.AuxFields

    # Get matrices
    Tn     = typeof(param.Mesh.nc)
    Ne,Qe, = getEdgeConstraints(Mesh)
    Msig   = getEdgeMassMatrix(Mesh,vec(sigma))
    Msig   = Ne'*Msig*Ne
    P      = param.Obs
    if invertMu || (param.storageLevel == :None)
        Nf,Qf, = getFaceConstraints(Mesh)
        Curl   = getCurlMatrix(Mesh)
        Curl   = Qf*Curl*Ne
        Mmu    = getFaceMassMatrix(Mesh,1./mu)
        Mmu    = Nf'*Mmu*Nf
        DmuinvDmu = spdiagm(-1./(mu.^2))
    else
        Curl = spzeros(Tf,Tn,0,0)
        Mmu  = spzeros(Tf,Tn,0,0)
    end
    K = getMaxwellCurlCurlMatrix!(param,model)

    #Initialize intermediate and output arrays
    ns  = size(s,2)
    ne  = size(ew,1)
    nt  = length(dt)
    nr  = size(P,2)
    lam = zeros(Tf,ne,ns,3)

    # Multiply by transpose of time interpolation matrix
    # To map from data space to data at all times space
    ptz = param.ObsTimes'*z

    # Check source type
    groundedSource = param.sourceType == :Galvanic ? true : false

    if groundedSource
        ptz  = reshape(ptz,nr,ns,nt+1)
        pz   = zeros(Tf,ne,ns,nt)
        for i=1:ns
            for j=1:nt
                pz[:,i,j] = P*ptz[:,i,j+1]
            end
        end
    else
        ptz  = reshape(ptz,nr,ns,nt)
        pz   = zeros(Tf,ne,ns,nt)
        for i=1:ns
            for j=1:nt
                pz[:,i,j] = P*ptz[:,i,j]
            end
        end
    end

    JTv   = Dict{String,Vector{Float64}}()
    if invertSigma
        JTv["sigmaCell"] = zeros(Tf,length(sigma))
    end
    if invertMu
        JTv["muCell"] = zeros(Tf,length(mu))
    end
    uniqueSteps = unique(dt)
    A           = speye(Tf,Tn,size(Ne,2)) #spzeros(T,0,0)
    iSolver     = 0
    dt          = [dt;dt[end];dt[end]]
    for i=nt:-1:1
        tau1 = i > 1 ? dt[i]/dt[i-1] : one(Tf)
        tau2 = dt[i+1]/dt[i]
        tau3 = dt[i+2]/dt[i+1]
        g1   = (1+2*tau1)/(1+tau1)
        g2   = 1 + tau1
        g3   = (tau1^2)/(1+tau1)
        g2p  = 1 + tau2
        g3p  = (tau3^2)/(1+tau3)
        if (dt[i] != dt[i+1]) || (i==nt)
            A,iSolver = getBDF2ConstDTmatrix!(dt[i],model,Msig,param,uniqueSteps)
        end
        if (i>1) && (dt[i] != dt[i-1])
            Atr  = K + (g1/dt[i])*Msig
            if EMsolvers[iSolver].doClear == 1
                clear!(EMsolvers[iSolver])
                factorLinearSystem!(A,EMsolvers[iSolver])
                EMsolvers[iSolver].doClear = 0
            end
            M(X) = solveLinearSystem!(A,X,similar(X),EMsolvers[iSolver])[1]
        end

        rhs  = @views pz[:,:,i] + Msig*(g2p*lam[:,:,2]/dt[i+1]-g3p*lam[:,:,3]/dt[i+2])
        if (i>1) && (dt[i] != dt[i-1])
            for j in 1:ns
                if norm(rhs) > 1e-20
                    lam[:,j,1],cgFlag,err,iterTmp, = cg(Atr,rhs[:,j],
                        x=vec(lam[:,j,1]),M=M,maxIter=20,tol=param.cgTol)
                    if cgFlag != 0
                        warn("getSensTMatVec: cg failed to converge at time step $i. Exited with flag $cgFlag. Reached residual $err with tolerance $(param.cgTol)")
                    end
                end
            end
        else
            lam[:,:,1],EMsolvers[iSolver] = solveMaxTimeBDF2ConstDT!(A,rhs,Msig,Mesh,dt[i],EMsolvers[iSolver])
            EMsolvers[iSolver].doClear = 0
        end

        for j = 1:ns
            if invertSigma
                if i>1
                    JTv["sigmaCell"] .-= @views dEdgeMassMatrixTrTimesVector(Mesh,sigma,
                      Ne*((g1.*ew[:,j,i+1].-g2.*ew[:,j,i] .+ g3.*ew[:,j,i-1])./dt[i]),
                      Ne*lam[:,j,1])
                else
                    JTv["sigmaCell"] .-= @views dEdgeMassMatrixTrTimesVector(Mesh,sigma,
                      Ne*((1.5.*ew[:,j,2].-0.75.*ehat[:,j].-0.75.*ew[:,j,1])./dt[i]),
                      Ne*lam[:,j,1])
                end
            end
            if invertMu
                JTv["muCell"] .-= @views DmuinvDmu*dFaceMassMatrixTrTimesVector(Mesh,mu,
                  Nf*(Curl*ew[:,j,i+1]),
                  Nf*(Curl*lam[:,j,1]))
            end
        end
        lam[:,:,3] = lam[:,:,2]
        lam[:,:,2] = lam[:,:,1]
    end

    #Do ehat stuff
    lmTmp = zeros(size(lam,1),ns)
    for j = 1:ns
        rhs = 3/(4*dt[1])*Msig*lam[:,j,2]
        lmTmp[:,j],EMsolvers[iSolver] = solveMaxTimeBDF2ConstDT!(A,rhs,Msig,
                                                  Mesh,dt[1],EMsolvers[iSolver])
        if invertSigma
            JTv["sigmaCell"] .-= @views dEdgeMassMatrixTrTimesVector(Mesh,sigma,
              Ne*(1.5.*(ehat[:,j].-ew[:,j,1])./dt[1]),
              Ne*lmTmp[:,j])
        end
        if invertMu
            JTv["muCell"] .-= @views DmuinvDmu*dFaceMassMatrixTrTimesVector(Mesh,mu,
              Nf*(Curl*ehat[:,j]),
              Nf*(Curl*lmTmp[:,j]))
        end
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
    if groundedSource && invertSigma
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
            JTv["sigmaCell"] .-= @views dEdgeMassMatrixTrTimesVector(Mesh,sigma,
                                   -Ne*ew[:,j,1],Ne*(G*lam0))
        end
        if param.storageLevel != :Factors
            clear!(DCsolver)
            DCsolver.doClear = 1
        end
    end
    return MaxwellTimeModel(JTv, model.activeInversionProperties)
end

#----------------------------------------------------------------------

function getSensTMatVecBDF2ConstDT{Tf<:Real}(z::Vector{Tf},model::MaxwellTimeModel,
                                   param::MaxwellTimeParam)
	#getSensTMatVec(z,sigma,param)
	#This function computes (dData/dsigma)^T*z for BDF2 time-stepping
	#forward problem. It handles grounded and inductive sources,
	#assuming DC data is integral of electric field and not a
	#potential difference (i.e. electric field at t=0 is known) for
	#grounded sources and e0=0 for inductive sources. FE is used
	#for first time step.

    # Unpack model into conductivity and magnetic permeability
    sigma       = model.values["sigmaCell"]
    mu          = model.values["muCell"]
    invertSigma = in("sigmaCell", model.activeInversionProperties)
    invertMu    = in(   "muCell", model.activeInversionProperties)

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
    Tn     = typeof(param.Mesh.nc)
    Ne,Qe, = getEdgeConstraints(M)
    Msig   = getEdgeMassMatrix(M,vec(sigma))
    Msig   = Ne'*Msig*Ne
    P      = param.Obs
    if invertMu || (param.storageLevel == :None)
        Nf,Qf, = getFaceConstraints(M)
        Curl   = getCurlMatrix(M)
        Curl   = Qf*Curl*Ne
        Mmu    = getFaceMassMatrix(M,1./mu)
        Mmu    = Nf'*Mmu*Nf
    else
        Curl = spzeros(Tf,Tn,0,0)
        Mmu  = spzeros(Tf,Tn,0,0)
    end

    #Initialize intermediate and output arrays
    ns  = size(s,2)
    ne  = size(ew,1)
    nt  = length(param.dt)
    nr  = size(P,2)
    lam = zeros(Tf,ne,ns,3)

    # Multiply by transpose of time interpolation matrix
    # To map from data space to data at all times space
    ptz = param.ObsTimes'*z

    # Check source type
    groundedSource = param.sourceType == :Galvanic ? true : false

    if groundedSource
      ptz  = reshape(ptz,nr,ns,nt+1)
      pz   = zeros(Tf,ne,ns,nt)
      for i=1:ns
        for j=1:nt
          pz[:,i,j] = P*ptz[:,i,j+1]
        end
      end
    else
      ptz  = reshape(ptz,nr,ns,nt)
      pz   = zeros(Tf,ne,ns,nt)
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
    A           = speye(Tf,Tn,size(Ne,2)) #spzeros(T,0,0)
    A,iSolver   = getBDF2ConstDTmatrix!(dt,model,Msig,param,uniqueSteps)
    for i=nt:-1:2
      for j = 1:ns
        rhs = pz[:,j,i] + 1/dt*Msig*(2*lam[:,j,2]-0.5*lam[:,j,3])
        if norm(rhs) > 1e-20
          lam[:,j,1],EMsolver = solveMaxTimeBDF2ConstDT!(A,rhs,Msig,M,dt,EMsolver)
          EMsolver.doClear = 0
        end
        Gzi = (1/dt)*getdEdgeMassMatrix(M,sigma,Ne*(1.5*ew[:,j,i+1]-2*ew[:,j,i]+
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
      Gzi                 = getdEdgeMassMatrix(M,sigma,1/dt*Ne*(1.5*ew[:,j,2]-0.75*ehat[:,j]-0.75*ew[:,j,1]))
      JTv                 = JTv - Gzi'*(Ne*lam[:,j,1])
      lam[:,j,3]          = lam[:,j,2]
      lam[:,j,2]          = lam[:,j,1]
    end

    #Do ehat stuff
    lmTmp = zeros(size(lam,1),ns)
    for j = 1:ns
      rhs = 3/(4*dt)*Msig*lam[:,j,2]
      lmTmp[:,j],EMsolver = solveMaxTimeBDF2ConstDT!(A,rhs,Msig,M,dt,EMsolver)
      Gzi                 = getdEdgeMassMatrix(M,sigma,3/(2*dt)*Ne*(ehat[:,j]-ew[:,j,1]))
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
            Gzi           = getdEdgeMassMatrix(M,sigma,-Ne*ew[:,j,1])
            JTv           = JTv - Gzi'*(Ne*G*lam0)
        end
        if param.storageLevel != :Factors
            clear!(DCsolver)
            DCsolver.doClear = 1
        end
    end
    MaxwellTimeModel(Dict("sigmaCell"=>JTv),["sigmaCell"])
end

#-------------------------------------------------------

function getSensTMatVecTRBDF2{Tf<:Real}(z::Vector{Tf},model::MaxwellTimeModel,
                                   param::MaxwellTimeParam)
    error("Sensitivities for TRBDF2 not implemented")
end

#----------------------------------------------------------------------

sensitivityTFunctions = Dict(zip(supportedIntegrationMethods,
                                [getSensTMatVecBE;getSensTMatVecBDF2;
                                 getSensTMatVecBDF2ConstDT;
                                 getSensTMatVecTRBDF2]))
