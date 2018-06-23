export getSensMatVec,interpGlobalToLocal


# For backward compatibility
function getSensMatVec{Tf<:Real}(z::Vector{Tf},sigma::Vector{Tf},
                       param::MaxwellTimeParam)
    zm = MaxwellTimeModel(z,[one(eltype(z))])
    m  = MaxwellTimeModel(sigma,fill(convert(eltype(sigma),4*pi*1e-7),
                          param.Mesh.nc))
    return getSensMatVec(zm,m,param)
end

"""
function getSensMatVec(DsigDmz::MaxwellTimeModel,model::MaxwellTimeModel,
                       param::MaxwellTimeParam)

This routine computes local forward problem inversion parametrization
independent contribution to sensitivity matrix vector products J*z. The full
contribution to J from a single forward problem is

J = (d D/d sigma)\*Mesh2Mesh\*(d sigma/d m) + (d D/d mu)\*Mesh2Mesh\*(d mu/d m)

where m is the vector of inversion model parameters. In a jInv inversion
Mesh2Mesh*(d sigma/d m) and Mesh2Mesh*(d mu/d m) are computed externally
using the modfun and interpGlobalToLocal functions.

The MaxwellTime specific version of interpGlobalToLocal returns a composite
object DsigDmz of type MaxwellTimeModel, which contains Mesh2Mesh*(d sigma/d m)
and Mesh2Mesh*(d mu/d m). DsigDmz is input to this function.

The input to this function DsigDmz is a composite object with two fields,
the vectors (d sigma/d m)\*z and (d mu/d m)\*z for some vector z in model space.

It returns J\*z. We save computing by including two additional fields in the
MaxwellTimeModel object, invertSigma and invertMu.
Recall that

(d D/dsigma) = (A^-1)(dC/dsigma)
(d D/dmu)    = (A^-1)(dC/dmu)

No need to compute dC/dmu if we're not inverting for mu, for example.

Note that we don't use a different forward problem type for explicit sensitivities.
This function checks value of the symbol param.sensitivityMethod. If it is set to
:Explicit then explicit sensitivities are computed, otherwise implicit sensitivities
are used.
"""
function getSensMatVec(x::MaxwellTimeModel,model::MaxwellTimeModel,
                       param::MaxwellTimeParam)

    activeProps = model.activeInversionProperties

    if param.sensitivityMethod == :Explicit
        if in("sigmaCell",activeProps) && in("muCell",activeProps)
            error("Using explicit sensitivities while inverting sigma and mu simultaneously is not supported")
        end
        mod = in("sigmaCell",activeProps) ? model.values["sigmaCell"] : model.values["muCell"]
        if isempty(param.Sens)
           #  error("Why are we here?")
            ndata     = size(param.ObsTimes,1)
            J         = zeros(ndata, length(mod))
            sensTFunc = sensitivityTFunctions[param.timeIntegrationMethod]
            v = Array{Float64}(ndata)
            for k=1:ndata
            	fill!(v, 0.0)
            	v[k]     = 1.0
            	JTvStruc = sensTFunc(v,model,param)
            	JTv      = in("sigmaCell",activeProps) ? JTvStruc.values["sigmaCell"] : JTvStruc.values["muCell"]
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
        z = in("sigmaCell",activeProps) ? x.values["sigmaCell"] : x.values["muCell"]
        return param.Sens*z
    else # Implicit sensitivities
        # Dispatch based on forward modelling integration method,
        # e.g. BE vs BDF2
        sensFunc = sensitivityFunctions[param.timeIntegrationMethod]
        Jx        = sensFunc(x,model,param)
        return Jx
    end

end


function getSensMatVecBE{Tf<:Real}(x::MaxwellTimeModel{Tf},
                                   model::MaxwellTimeModel,
                                   param::MaxwellTimeParam)

    # getSensMatVecBE(z,sigma,param)
    # This function computes (dData/dsigma)*DsigDmz + (dData/dmu)*DmuDmz
    # for BE time-stepping It handles grounded and inductive sources,
    # assuming DC data is integral of electric field and not a
    # potential difference (i.e. electric field at t=0 is known) for
    # grounded sources and e0=0 for inductive sources.
    # For grounded sources
    #  dCdu = |G'*Msig*G                                    |
    #         |1/dt*Msig*G    K+1/dt*Msig                   |
    #         |               -1/dt*Msig    K+1/dt*Msig  ...|
    #         |                     .                       |
    #         |                     .                       |
    #         |                     .                       |

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
    Tn       = typeof(param.Mesh.nc)
    Ne,Qe,   = getEdgeConstraints(Mesh)
    Msig     = getEdgeMassMatrix(Mesh,vec(sigma))
    Msig     = Ne'*Msig*Ne
    P        = param.Obs
    DrhoDsig = param.modUnits == :res ? spdiagm(-(sigma.^2)) : UniformScaling(one(Tf))
    if invertMu || (param.storageLevel == :None)
        Nf,Qf,    = getFaceConstraints(Mesh)
        Curl      = getCurlMatrix(Mesh)
        Curl      = Qf*Curl*Ne
        Mmu       = getFaceMassMatrix(Mesh,1./mu)
        Mmu       = Nf'*Mmu*Nf
        DmuinvDmu = spdiagm(-1./(mu.^2))
    end

    if invertMu & (length(mu)>3*Mesh.nc)
        error("MaxwellTime:getSensMatVec: Inverting fully anisotropic mu not supported")
    end

    #Initialize intermediate and output arrays
    ns = size(s,2)
    ne = size(ew,1)
    nt = length(dt)
    lam = zeros(Tf,ne,ns)  # only one lam is needed

    Jv  = zeros(Tf,size(P,2),ns,nt)

    #Do the DC part of dCdm if source is grounded and we're inverting for sigma.
    #Note that this code assumes (for grounded sources) that DC data
    #Are computed as the integral of electric field over a receiver
    #dipole and not as a potential difference
    groundedSource = param.sourceType == :Galvanic ? true : false
    if groundedSource && invertSigma
        DCsolver = param.DCsolver
        Nn,Qn, = getNodalConstraints(Mesh)
        Gin    = getNodalGradientMatrix(Mesh)
        G      = Qe*Gin*Nn
        A      = getDCmatrix(Msig,G,param) # Forms matrix only if needed
        Jvdc   = zeros(Tf,size(P,2),ns)
        rhs    = zeros(Tf,size(G,2),ns)
        for j = 1:ns
            Gzi      = G'*Ne'*getdEdgeMassMatrix(Mesh,sigma,-Ne*ew[:,j,1])
            rhs[:,j] = Gzi*DrhoDsig*x.values["sigmaCell"]
        end
        lam0,DCsolver = solveDC!(A,rhs,DCsolver)
        lam    = -G*lam0 #Taking gradient of lam0
                         #Prepares for data projection and
                         #preps lam for use as rhs in first
                         #time-step
        Jvdc          = -P'*lam
        DCsolver.doClear = 0
        if param.storageLevel != :Factors
            clear!(DCsolver)
            DCsolver.doClear = 1
        end
    elseif groundedSource && ~invertSigma
       Jvdc   = zeros(Tf,size(P,2),ns)
    end

    if invertMu

    end

    iSolver     = 0
    uniqueSteps = unique(dt)
    Tn          = typeof(param.Mesh.nc)
    A           = speye(Tf,Tn,size(Ne,2)) #spzeros(T,0,0)
    for i=1:nt
        # Form A when needed. It is not stored or formed if
        # param.storageLevel == :Factors
        if ( (i==1) || (dt[i] != dt[i-1]) )
            A,iSolver = getBEMatrix!(dt[i],model,Msig,param,uniqueSteps)
        end
        rhs = (1/dt[i])*(Msig*lam)
        if invertSigma
            addDCDsigmaBE!(Mesh,sigma,Ne,ew[:,:,i:i+1],dt[i],DrhoDsig*x.values["sigmaCell"],ns,rhs)
        end
        if invertMu
            addDCDmuBE!(Mesh,mu,Nf,Curl,ew[:,:,i+1],DmuinvDmu*x.values["muCell"],ns,rhs)
        end

        lam,EMsolvers[iSolver] = solveMaxTimeBE!(A,rhs,Msig,Mesh,dt,i,
                                         storageLevel,EMsolvers[iSolver])

        # compute Jv
        Jv[:,:,i]  = -P'*lam
       # lam[:,:,1] = lam[:,:,2]
    end
    if storageLevel != :Factors
        for solver in EMsolvers
          clear!(solver)
          solver.doClear = 1
        end
    end
    #Jv = (groundedSource && invertSigma) ? [vec(Jvdc);vec(Jv)] : vec(Jv)
    Jv = groundedSource ? [vec(Jvdc);vec(Jv)] : vec(Jv)
    return param.ObsTimes*Jv
end

# Prepare sigma and mu contributions to rhs efficiently. See commented out code
# below for less efficient but more mathematically transparent and numerically
# equivalent code.
function addDCDsigmaBE!(Mesh,sigma,Ne,ew,dt,DsigDmz,ns,rhs)
    Gzi  = zeros(size(Ne,1))
    tmpv = zeros(size(Ne,2))
    for j = 1:ns
        Gzi = getdEdgeMassMatrix(Mesh,sigma,Ne*(ew[:,j,2]-ew[:,j,1]))*DsigDmz
        rhs[:,j] += (1/dt)*At_mul_B!(tmpv,Ne,Gzi)
    end
end

function addDCDmuBE!(Mesh,mu,Nf,Curl,ew,dz,ns,rhs)
    Curle       = zeros(size(Curl,1))
    NfCurle     = zeros(size(Nf,1))
    Gzi         = zeros(size(Nf,1))
    NftGzi      = zeros(size(Nf,2))
    CurltNftGzi = zeros(size(Curl,2))
    for j = 1:ns
        A_mul_B!(Curle,Curl,ew[:,j])
        A_mul_B!(NfCurle,Nf,Curle)
        G = getdFaceMassMatrix(Mesh,mu,NfCurle)
        A_mul_B!(Gzi,G,dz)
        At_mul_B!(NftGzi,Nf,Gzi)
        rhs[:,j] += At_mul_B!(CurltNftGzi,Curl,NftGzi)
    end
end

# for j = 1:ns
#     if invertSigma
#       Gzi           = (1/dt[i])*Ne'*getdEdgeMassMatrix(M,sigma,Ne*(ew[:,j,i+1]-
#                                                     ew[:,j,i]))
#       rhsSigma[:,j] = Gzi*DsigDmz
#     end
#     if invertMu
#         Gzi        = Curl'*Nf'*getdFaceMassMatrix(M,1./mu,Nf*Curl*ew[:,j,i+1])*
#                      DmuinvDmu
#         rhsMu[:,j] = Gzi*DmuDmz
#     end
# end

#-------------------------------------------------------

function getSensMatVecBDF2{Tf<:Real}(x::MaxwellTimeModel{Tf},
                                    model::MaxwellTimeModel,
                                    param::MaxwellTimeParam)
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
    # Unpack model into conductivity and magnetic permeability
    sigma       = model.values["sigmaCell"]
    mu          = model.values["muCell"]
    invertSigma = in("sigmaCell", model.activeInversionProperties)
    invertMu    = in("muCell", model.activeInversionProperties)

    #Unpack param
    Mesh         = param.Mesh
    storageLevel = param.storageLevel
    EMsolvers    = param.EMsolvers
    dt           = param.dt
    ew           = param.fields
    ehat         = param.AuxFields
    s            = param.Sources

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
        DmuinvDmz = -x.values["muCell"]./(mu.^2)
        # Equivalent to
        # DmuinvDmu = spdiagm(-1./(mu.^2))
        # DmuinvDmz = DmuinvDmu*x.values["muCell"]
    else
        Curl = spzeros(Tf,Tn,0,0)
        Mmu  = spzeros(Tf,Tn,0,0)
    end
    K = getMaxwellCurlCurlMatrix!(param,model)

    #Initialize intermediate and output arrays
    ns = size(s,2)
    ne = size(ew,1)
    nt = length(param.dt)
    lam = zeros(Tf,ne,ns,3)
    Jv  = zeros(Tf,size(P,2),ns,nt)

    #Do the DC part of dCdm if source is grounded and we're inverting for sigma.
    #Note that this code assumes (for grounded sources) that DC data
    #Are computed as the integral of electric field over a receiver
    #dipole and not as a potential difference
    groundedSource = param.sourceType == :Galvanic ? true : false
    if groundedSource && invertSigma
        DCsolver = param.DCsolver
        Nn,Qn, = getNodalConstraints(Mesh)
        Gin    = getNodalGradientMatrix(Mesh)
        G      = Qe*Gin*Nn
        A      = getDCmatrix(Msig,G,param) # Forms matrix only if needed
        Jvdc   = zeros(Tf,size(P,2),ns)
        for j = 1:ns
            Gzi = G'*Ne'*getdEdgeMassMatrix(Mesh,sigma,-Ne*ew[:,j,1])
            rhs = Gzi*x.values["sigmaCell"]
            lam0,DCsolver    = solveDC!(A,rhs,DCsolver)
            lam[:,j,1]       = -G*lam0 #Taking gradient of lam0
                                       #Prepares for data projection and
                                       #preps lam for use as rhs in first
                                       #time-step
            Jvdc[:,j]        = -P'*view(lam,:,j,1)
            DCsolver.doClear = 0
        end
        if param.storageLevel != :Factors
            clear!(DCsolver)
            DCsolver.doClear = 1
        end
    end

    #Do BE and interpolation for first time-step
    uniqueSteps = unique(dt)
    A           = speye(Tf,Tn,size(Ne,2)) # spzeros(T,0,0)
    A,iSolver   = getBDF2ConstDTmatrix!(dt[1],model,Msig,param,uniqueSteps)
    rhs = 3/(2*dt[1])*Msig*lam[:,:,1]
    for j = 1:ns
        if invertSigma
            @views begin rhs[:,j] .+= 3/(2*dt[1]).*(Ne'*dEdgeMassMatrixTimesVector(Mesh,sigma,
                                  Ne*(ehat[:,j]-ew[:,j,1]),
                                  x.values["sigmaCell"])) end
        end
        if invertMu
            @views begin rhs[:,j] .+= Curl'*(Nf'*dFaceMassMatrixTimesVector(Mesh,1./mu,
                           Nf*(Curl*ehat[:,j]),
                           DmuinvDmz)) end
        end
    end
    lmTmp,EMsolvers[1]      = solveMaxTimeBDF2ConstDT!(A,rhs,Msig,Mesh,dt[1],EMsolvers[1])
    EMsolvers[1].doClear    = 0

    rhs = 3/(4*dt[1])*Msig*lam[:,:,1] + 3/(4*dt[1])*Msig*lmTmp
    for j = 1:ns
        if invertSigma
            @views begin rhs[:,j] .+= Ne'*dEdgeMassMatrixTimesVector(Mesh,sigma,
                         Ne*((1.5.*ew[:,j,2] .- 0.75.*ehat[:,j] .- 0.75.*ew[:,j,1])./dt[1]),
                         x.values["sigmaCell"]) end
        end
        if invertMu
            @views begin rhs[:,j] .+= Curl'*(Nf'*dFaceMassMatrixTimesVector(Mesh,1./mu,
                           Nf*(Curl*ew[:,j,2]),
                           DmuinvDmz)) end
        end
    end
    lam[:,:,2],EMsolvers[1] = solveMaxTimeBDF2ConstDT!(A,rhs,Msig,Mesh,dt[1],EMsolvers[1])
    Jv[:,:,1]               = -P'*view(lam,:,:,2)

    #Do the rest of the time-steps
    for i=2:nt
        if dt[i] != dt[i-1]
            A,iSolver = getBDF2ConstDTmatrix!(dt[i],model,Msig,param,uniqueSteps)
            if EMsolvers[iSolver].doClear == 1
                clear!(EMsolvers[iSolver])
                factorLinearSystem!(A,EMsolvers[iSolver])
                EMsolvers[iSolver].doClear = 0
            end
            M(X) = solveLinearSystem!(A,X,similar(X),EMsolvers[iSolver])[1]
            tau = dt[i]/dt[i-1]
            g1  = (1+2*tau)/(1+tau)
            g2  = 1 + tau
            g3  = (tau^2)/(1+tau)
            Atr = K + (g1/dt[i])*Msig
            rhs = Msig*((g2.*lam[:,:,2].-g3.*lam[:,:,1])./dt[i])
            for j = 1:ns
                if invertSigma
                    @views begin rhs[:,j] .+= Ne'*dEdgeMassMatrixTimesVector(Mesh,sigma,
                                Ne*((g1.*ew[:,j,i+1].-g2.*ew[:,j,i].+g3.*ew[:,j,i-1])./dt[i]),
                                x.values["sigmaCell"]) end
                end
                if invertMu
                    @views begin rhs[:,j] .+= Curl'*(Nf'*dFaceMassMatrixTimesVector(Mesh,1./mu,
                                   Nf*(Curl*ew[:,j,i+1]),
                                   DmuinvDmz)) end
                end
                lam[:,j,3],cgFlag,err,iterTmp, = cg(Atr,rhs[:,j],
                   x=vec(lam[:,j,3]),M=M,maxIter=20,tol=param.cgTol)
                if cgFlag != 0
                    warn("getSensMatVec: cg failed to converge at time step $i. Reached residual $err with tolerance $(param.cgTol)")
                end
                # compute Jv
                @views Jv[:,j,i]  = -P'*(lam[:,j,3])
                @views lam[:,j,1] = lam[:,j,2]
                @views lam[:,j,2] = lam[:,j,3]
            end
        else
            @views rhs = Msig*((2.0.*lam[:,:,2].-0.5.*lam[:,:,1])./dt[i])
            for j = 1:ns
                if invertSigma
                    @views begin rhs[:,j] .+= Ne'*dEdgeMassMatrixTimesVector(Mesh,sigma,
                                Ne*((1.5.*ew[:,j,i+1].-2.0.*ew[:,j,i].+0.5.*ew[:,j,i-1])./dt[i]),
                                x.values["sigmaCell"]) end
                end
                if invertMu
                    @views begin rhs[:,j] .+= Curl'*(Nf'*dFaceMassMatrixTimesVector(Mesh,1./mu,
                                   Nf*(Curl*ew[:,j,i+1]),
                                   DmuinvDmz)) end
                end
            end
       	    lam[:,:,3],EMsolvers[iSolver] = solveMaxTimeBDF2ConstDT!(A,rhs,Msig,Mesh,dt[i],EMsolvers[iSolver])
            # compute Jv
            @views Jv[:,:,i]  = -P'*(lam[:,:,3])
            @views lam[:,:,1] = lam[:,:,2]
            @views lam[:,:,2] = lam[:,:,3]
        end
    end
    if storageLevel != :Factors
        for solver in EMsolvers
          clear!(solver)
          solver.doClear = 1
        end
    end
    Jv = groundedSource ? [vec(Jvdc);vec(Jv)] : vec(Jv)
    return param.ObsTimes*Jv
end
#-------------------------------------------------------


function getSensMatVecTRBDF2{Tf<:Real}(DsigDmz::Vector{Tf},DmuDmz::Vector{Tf},
                                    model::MaxwellTimeModel,param::MaxwellTimeParam)
    error("Sensitivities for TRBDF2 not implemented")
end
#-------------------------------------------------------

function getSensMatVecBDF2ConstDT{Tf<:Real}(x::MaxwellTimeModel{Tf},
                                  model::MaxwellTimeModel,param::MaxwellTimeParam)
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
    # Unpack model into conductivity and magnetic permeability
    sigma       = model.values["sigmaCell"]
    mu          = model.values["muCell"]
    invertSigma = in("sigmaCell", model.activeInversionProperties)
    invertMu    = in(   "muCell", model.activeInversionProperties)

    #Unpack param
    M            = param.Mesh
    storageLevel = param.storageLevel
    EMsolver     = param.EMsolvers[1]
    dt           = param.dt[1]
    ew           = param.fields
    ehat         = param.AuxFields
    s            = param.Sources

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
    ns = size(s,2)
    ne = size(ew,1)
    nt = length(param.dt)
    lam = zeros(Tf,ne,ns,3)
    Jv  = zeros(Tf,size(P,2),ns,nt)

    #Do the DC part of dCdm if source is grounded and we're inverting for sigma.
    #Note that this code assumes (for grounded sources) that DC data
    #Are computed as the integral of electric field over a receiver
    #dipole and not as a potential difference
    groundedSource = param.sourceType == :Galvanic ? true : false
    if groundedSource && invertSigma
        DCsolver = param.DCsolver
        Nn,Qn, = getNodalConstraints(M)
        Gin    = getNodalGradientMatrix(M)
        G      = Qe*Gin*Nn
        A      = getDCmatrix(Msig,G,param) # Forms matrix only if needed
        Jvdc   = zeros(Tf,size(P,2),ns)
        for j = 1:ns
            Gzi = G'*Ne'*getdEdgeMassMatrix(M,sigma,-Ne*ew[:,j,1])
            rhs = Gzi*x.values["sigmaCell"]
            lam0,DCsolver    = solveDC!(A,rhs,DCsolver)
            lam[:,j,1]       = -G*lam0 #Taking gradient of lam0
                                       #Prepares for data projection and
                                       #preps lam for use as rhs in first
                                       #time-step
            Jvdc[:,j]        = -P'*lam[:,j,1]
            DCsolver.doClear = 0
        end
        if param.storageLevel != :Factors
            clear!(DCsolver)
            DCsolver.doClear = 1
        end
    end

    #Do BE and interpolation for first time-step
    if invertMu
        error("Inverting for mu with bdf2 time-stepping not yet supported")
    end
    uniqueSteps = [dt]
    A           = speye(Tf,Tn,size(Ne,2)) #spzeros(T,0,0)
    A,iSolver   = getBDF2ConstDTmatrix!(dt,model,Msig,param,uniqueSteps)
    for j = 1:ns
        Gzi                 = 3/(2*dt)*Ne'*getdEdgeMassMatrix(M,sigma,Ne*(ehat[:,j]-ew[:,j,1]))
        rhs                 = Gzi*x.values["sigmaCell"] + 3/(2*dt)*Msig*lam[:,j,1]
        lmTmp,EMsolver      = solveMaxTimeBDF2ConstDT!(A,rhs,Msig,M,dt,EMsolver)
        EMsolver.doClear    = 0
        Gzi                 = Ne'*getdEdgeMassMatrix(M,sigma,1/dt*Ne*(1.5*ew[:,j,2]-0.75*ehat[:,j]-0.75*ew[:,j,1]))
        rhs                 = Gzi*x.values["sigmaCell"] + 3/(4*dt)*Msig*lam[:,j,1] + 3/(4*dt)*Msig*lmTmp
        lam[:,j,2],EMsolver = solveMaxTimeBDF2ConstDT!(A,rhs,Msig,M,dt,EMsolver)
        Jv[:,j,1]           = -P'*lam[:,j,2]
    end

    #Do the rest of the time-steps
    for i=2:nt
        for j = 1:ns
       	       Gzi = (1/dt)*Ne'*getdEdgeMassMatrix(M,sigma,Ne*(1.5*ew[:,j,i+1]-
       	              2*ew[:,j,i]+0.5*ew[:,j,i-1]))
       	       rhs = Gzi*x.values["sigmaCell"] + 1/dt*Msig*(2*lam[:,j,2]-0.5*lam[:,j,1])
       	       lam[:,j,3],EMsolver = solveMaxTimeBDF2ConstDT!(A,rhs,Msig,M,dt,EMsolver)
            # compute Jv
            Jv[:,j,i]  = -P'*(lam[:,j,3])
            lam[:,j,1] = lam[:,j,2]
            lam[:,j,2] = lam[:,j,3]
        end
    end
    if storageLevel != :Factors
        clear!(EMsolver)
        EMsolver.doClear = 1
    end
    Jv = groundedSource ? [vec(Jvdc);vec(Jv)] : vec(Jv)
    return param.ObsTimes*Jv
end

#-------------------------------------------------------

# function getBEMatrix{T<:Real,N}(Msig::SparseMatrixCSC{T,N},Mmu::SparseMatrixCSC{T,N},
#                        Curl::SparseMatrixCSC{T,N},
#                        dt::T,iStepSize::N,param::MaxwellTimeParam)
#
#     storageLevel = param.storageLevel
#     if storageLevel == :Matrices
#         matNum = param.sourceType == :Galvanic ? iStepSize+1 : iStepSize
#         A = param.Matrices[matNum]
#     elseif storageLevel == :None
#         A = Curl'*Mmu*Curl + (1/dt)*Msig
#     else
#         A = spzeros(T,0,0) # Empty sparse matrix placeholder argument
#     end
#     return A
# end

# function getBDF2ConstDTmatrix{T<:Real,N}(Msig::SparseMatrixCSC{T,N},
#                                          Mmu::SparseMatrixCSC{T,N},
#                                          Curl::SparseMatrixCSC{T,N},
#                                          dt::T,param::MaxwellTimeParam)
#
#     storageLevel = param.storageLevel
#     if storageLevel == :Matrices
#         matNum = param.sourceType == :Galvanic ? 2 : 1
#         A = param.Matrices[matNum]
#     elseif storageLevel == :None
#         A = A = Curl'*Mmu*Curl + 3/(2*dt)*Msig
#     else
#         A = spzeros(T,0,0) # Empty sparse matrix placeholder argument
#     end
#     return A
# end

sensitivityFunctions = Dict(zip(supportedIntegrationMethods,
                                [getSensMatVecBE;getSensMatVecBDF2;
                                 getSensMatVecBDF2ConstDT;
                                 getSensMatVecTRBDF2]))
