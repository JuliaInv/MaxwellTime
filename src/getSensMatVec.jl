export getSensMatVec,interpGlobalToLocal


# For backward compatibility
function getSensMatVec{T<:Real}(z::Vector{T},sigma::Vector{T},
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

J = (d D/d sigma)*Mesh2Mesh*(d sigma/d m) + (d D/d mu)*Mesh2Mesh*(d mu/d m)

where m is the vector of inversion model parameters. In a jInv inversion 
Mesh2Mesh*(d sigma/d m) and Mesh2Mesh*(d mu/d m) are computed externally
using the modfun and interpGlobalToLocal functions.

The MaxwellTime specific version of interpGlobalToLocal returns a composite
object DsigDmz of type MaxwellTimeModel, which contains Mesh2Mesh*(d sigma/d m)
and Mesh2Mesh*(d mu/d m). DsigDmz is input to this function.

The input to this function DsigDmz is a composite object with two fields,
the vectors (d sigma/d m)*z and (d mu/d m)*z for some vector z in model space.

It returns J*z. We save computing by including two additional fields in the
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
function getSensMatVec(DsigDmz::MaxwellTimeModel,model::MaxwellTimeModel,
                       param::MaxwellTimeParam)

    zsig = DsigDmz.sigma
    zmu  = DsigDmz.mu
                       
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
	z = invertSigma ? zsig : zmu
	return param.Sens*z
    else # Implicit sensitivities
        # Dispatch based on forward modelling integration method,
        # e.g. BE vs BDF2
        sensFunc = sensitivityFunctions[param.timeIntegrationMethod]
        Jz       = sensFunc(zsig,zmu,model,param)
        return Jz
    end
    
end


function getSensMatVecBE{T<:Real}(DsigDmz::Vector{T},DmuDmz::Vector{T},
                                  model::MaxwellTimeModel,param::MaxwellTimeParam)
                       
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
    sigma       = model.sigma
    mu          = model.mu
    invertSigma = model.invertSigma
    invertMu    = model.invertMu
    
    #Unpack param
    M             = param.Mesh
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
        Nf,Qf,    = getFaceConstraints(M)
        Curl      = getCurlMatrix(M)
        Curl      = Qf*Curl*Ne
        Mmu       = getFaceMassMatrix(M,1./mu)
        Mmu       = Nf'*Mmu*Nf
        DmuinvDmu = spdiagm(-1./(mu.^2))
    end
    if param.storageLevel == :None
        K = Curl'*Mmu*Curl
    else
        K = spzeros(T,0,0)
    end
    
    if invertMu & (length(mu)>3*M.nc)
        error("MaxwellTime:getSensMatVec: Inverting fully anisotropic mu not supported")
    end
    
    #Initialize intermediate and output arrays
    ns = size(s,2)
    ne = size(ew,1)
    nt = length(dt)
    lam = zeros(T,ne,ns,2)
    Jv  = zeros(T,size(P,2),ns,nt)
    
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
        Jvdc   = zeros(T,size(P,2),ns)
        rhs    = zeros(T,size(G,2),ns)
        for j = 1:ns
            Gzi      = G'*Ne'*getdEdgeMassMatrix(M,sigma,-Ne*ew[:,j,1])
            rhs[:,j] = Gzi*DsigDmz
        end
        lam0,DCsolver = solveDC!(A,rhs,DCsolver)
        lam[:,:,1]    = -G*lam0 #Taking gradient of lam0
                                #Prepares for data projection and
                                #preps lam for use as rhs in first 
                                #time-step
        Jvdc          = -P'*lam[:,:,1]
        DCsolver.doClear = 0
        if param.storageLevel != :Factors
            clear!(DCsolver)
            DCsolver.doClear = 1
        end
    end
    
    if invertMu
        
    end
    
    iSolver     = 0
    uniqueSteps = unique(dt)
    A           = spzeros(T,0,0)
    rhsSigma    = invertSigma ? zeros(T,ne,ns) : 0
    rhsMu       = invertMu ? zeros(T,ne,ns) : 0
    for i=1:nt
        # Form A when needed. It is not stored or formed if
        # param.storageLevel == :Factors
        if ( (i==1) || (dt[i] != dt[i-1]) )        
            A,iSolver = getBEMatrix!(dt[i],A,K,Msig,param,uniqueSteps)
        end
        for j = 1:ns
            if invertSigma
     	        Gzi           = (1/dt[i])*Ne'*getdEdgeMassMatrix(M,sigma,Ne*(ew[:,j,i+1]-
                                                            ew[:,j,i]))
     	        rhsSigma[:,j] = Gzi*DsigDmz + 1/dt[i]*Msig*lam[:,j,1]
            end
            if invertMu
                Gzi        = Curl'*Nf'*getdFaceMassMatrix(M,1./mu,Nf*Curl*ew[:,j,i+1])*
                             DmuinvDmu
                rhsMu[:,j] = Gzi*DmuDmz
            end
            rhs = rhsSigma + rhsMu
        end
     	lam[:,:,2],EMsolvers[iSolver] = 
     	    solveMaxTimeBE!(A,rhs,Msig,M,dt,i,storageLevel,
     	                    EMsolvers[iSolver])
        
        # compute Jv
        Jv[:,:,i]  = -P'*(lam[:,:,2])
        lam[:,:,1] = lam[:,:,2]
    end
    if storageLevel != :Factors
        for solver in EMsolvers
          clear!(EMsolvers)
          solver.doClear = 1
        end
    end  
    Jv = (groundedSource && invertSigma) ? [vec(Jvdc);vec(Jv)] : vec(Jv)
    return param.ObsTimes*Jv
end

#-------------------------------------------------------

function getSensMatVecBDF2{T<:Real}(DsigDmz::Vector{T},DmuDmz::Vector{T},
                                    model::MaxwellTimeModel,param::MaxwellTimeParam)
    error("Sensitivities for variable step-size bdf2 not implemented")
end
#-------------------------------------------------------


function getSensMatVecTRBDF2{T<:Real}(DsigDmz::Vector{T},DmuDmz::Vector{T},
                                    model::MaxwellTimeModel,param::MaxwellTimeParam)
    error("Sensitivities for TRBDF2 not implemented")
end
#-------------------------------------------------------

function getSensMatVecBDF2ConstDT{T<:Real}(DsigDmz::Vector{T},DmuDmz::Vector{T},
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
    sigma       = model.sigma
    mu          = model.mu
    invertSigma = model.invertSigma
    invertMu    = model.invertMu
    
    #Unpack param
    M            = param.Mesh
    storageLevel = param.storageLevel
    EMsolver     = param.EMsolvers[1]
    dt           = param.dt[1]
    ew           = param.fields
    ehat         = param.AuxFields
    s            = param.Sources

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
    ns = size(s,2)
    ne = size(ew,1)
    nt = length(param.dt)
    lam = zeros(T,ne,ns,3)
    Jv  = zeros(T,size(P,2),ns,nt)
    
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
        Jvdc   = zeros(T,size(P,2),ns)
        for j = 1:ns
            Gzi = G'*Ne'*getdEdgeMassMatrix(M,sigma,-Ne*ew[:,j,1])
            rhs = Gzi*DsigDmz
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
    A = getBDF2ConstDTmatrix(Msig,Mmu,Curl,dt,param)
    for j = 1:ns
      Gzi                 = 3/(2*dt)*Ne'*getdEdgeMassMatrix(M,Ne*(ehat[:,j]-ew[:,j,1])) 
      rhs                 = Gzi*DsigDmz + 3/(2*dt)*Msig*lam[:,j,1]
      lmTmp,EMsolver      = solveMaxTimeBDF2ConstDT!(A,rhs,Msig,M,dt,EMsolver)
      EMsolver.doClear    = 0
      Gzi                 = Ne'*getdEdgeMassMatrix(M,1/dt*Ne*(1.5*ew[:,j,2]-0.75*ehat[:,j]-0.75*ew[:,j,1]))
      rhs                 = Gzi*DsigDmz + 3/(4*dt)*Msig*lam[:,j,1] + 3/(4*dt)*Msig*lmTmp
      lam[:,j,2],EMsolver = solveMaxTimeBDF2ConstDT!(A,rhs,Msig,M,dt,EMsolver)
      Jv[:,j,1]           = -P'*lam[:,j,2]
    end
    
    #Do the rest of the time-steps
    for i=2:nt
      for j = 1:ns
     	Gzi = (1/dt)*Ne'*getdEdgeMassMatrix(M,Ne*(1.5*ew[:,j,i+1]-
     	       2*ew[:,j,i]+0.5*ew[:,j,i-1]))
     	rhs = Gzi*DsigDmz + 1/dt*Msig*(2*lam[:,j,2]-0.5*lam[:,j,1])
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
    return Jv
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

function getBDF2ConstDTmatrix{T<:Real,N}(Msig::SparseMatrixCSC{T,N},
                                         Mmu::SparseMatrixCSC{T,N},
                                         Curl::SparseMatrixCSC{T,N},
                                         dt::T,param::MaxwellTimeParam)

    storageLevel = param.storageLevel
    if storageLevel == :Matrices
        matNum = param.sourceType == :Galvanic ? 2 : 1
        A = param.Matrices[matNum]
    elseif storageLevel == :None
        A = A = Curl'*Mmu*Curl + 3/(2*dt)*Msig
    else
        A = spzeros(T,0,0) # Empty sparse matrix placeholder argument
    end
    return A
end

sensitivityFunctions = Dict(zip(supportedIntegrationMethods,
                                [getSensMatVecBE;getSensMatVecBDF2;
                                 getSensMatVecBDF2ConstDT;
                                 getSensMatVecTRBDF2]))
