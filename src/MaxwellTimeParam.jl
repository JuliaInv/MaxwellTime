export MaxwellTimeParam, getMaxwellTimeParam

# MaxwellTime forward problem type. Symbols are used to select categorical options.
# Inner and outer constructors are required to allow instantiating the type with
# last three fields left uninitialized. Use getMaxwellTimeParam function
# for instantiation. It validates input and sets things to default values.
mutable struct MaxwellTimeParam{Tf<:Real,Tn<:Integer,TSem<:AbstractSolver,TSdc<:AbstractSolver} <: ForwardProbType
    Mesh::AbstractMesh
    Sources::AbstractArray{Tf}
    Obs::AbstractArray{Tf}
    ObsTimes::AbstractArray{Tf}
    t0::Tf
    dt::Vector{Tf}
    wave::Vector{Tf}
    sourceType::Symbol
    storageLevel::Symbol
    sensitivityMethod::Symbol
    timeIntegrationMethod::Symbol
    EMsolvers::Vector{TSem}
    DCsolver::TSdc
    modUnits::Symbol
    K::SparseMatrixCSC{Tf,Tn}
    Matrices::Vector{SparseMatrixCSC}
    fields::Array{Float64}
    Sens::Array{Tf}
    cgTol::Tf  #Only used by BDF2. Tolerance for cg solves used when stepsize changes
    AuxFields::Array{Tf} # Stores extra electric fields outside
                        # of regular time-stepping. Currently used
                        # to initialize BDF2 time stepping without additional BE
                        # factorization. Needs to be stored for use in sensitivity
                        # computations

    # Use custom constructors to allow creation of objects with
    # Matrices, fields, explicit sensitivity matrix, and possibly K
    # left uninitialized
    # MaxwellTimeParam{Tf,TSem,TSdc}(Mesh::AbstractMesh,Sources::AbstractArray{Tf},Obs::AbstractArray{Tf},
    #                  ObsTimes::AbstractArray{Tf},t0::Tf,dt::Vector{Tf},wave::Vector{Tf},sourceType::Symbol,
    #                  storageLevel::Symbol,sensitivityMethod::Symbol,
    #                  timeIntegrationMethod::Symbol,EMsolvers::Vector{TSem},
    #                  DCsolver::TSdc,
    #                  modUnits::Symbol) where {Tf<:Real,TSem<:AbstractSolver,TSdc<:AbstractSolver} = new(
    #                                     Mesh,Sources,Obs,ObsTimes,t0,dt,wave,sourceType,
    #                                     storageLevel,sensitivityMethod,
    #                                     timeIntegrationMethod,
    #                                     EMsolvers,DCsolver,modUnits)

    MaxwellTimeParam{Tf,Tn,TSem,TSdc}(Mesh::AbstractMesh,Sources::AbstractArray{Tf},Obs::AbstractArray{Tf},
                     ObsTimes::AbstractArray{Tf},t0::Tf,dt::Vector{Tf},wave::Vector{Tf},sourceType::Symbol,
                     storageLevel::Symbol,sensitivityMethod::Symbol,
                     timeIntegrationMethod::Symbol,EMsolvers::Vector{TSem},
                     DCsolver::TSdc,modUnits::Symbol,K::SparseMatrixCSC{Tf,Tn}) where
                     {Tf<:Real,Tn<:Integer,TSem<:AbstractSolver,TSdc<:AbstractSolver} = new(
                                     Mesh,Sources,Obs,ObsTimes,t0,dt,wave,
                                     sourceType,storageLevel,sensitivityMethod,
                                     timeIntegrationMethod,
                                     EMsolvers,DCsolver,modUnits,K)
end

# Unfortunately parametric types need these matching outer and inner constructors
# MaxwellTimeParam{Tf,TSem,TSdc}(Mesh::AbstractMesh,Sources::AbstractArray{Tf},Obs::AbstractArray{Tf},
# 	                 ObsTimes::AbstractArray{Tf},t0::Tf,dt::Vector{Tf},wave::Vector{Tf},sourceType::Symbol,
# 	                 storageLevel::Symbol,sensitivityMethod::Symbol,
# 	                 timeIntegrationMethod::Symbol,EMsolvers::Vector{TSem},
# 	                 DCsolver::TSdc,modUnits::Symbol) = MaxwellTimeParam{Tf,TSem,TSdc}(
# 	                                                 Mesh,Sources,Obs,ObsTimes,
#                                                      t0,dt,wave,
# 	                                                 sourceType,storageLevel,
# 	                                                 sensitivityMethod,
#                                                      timeIntegrationMethod,
#                                                      EMsolvers,DCsolver,modUnits)

MaxwellTimeParam{Tf,Tn,TSem,TSdc}(Mesh::AbstractMesh,Sources::AbstractArray{Tf},Obs::AbstractArray{Tf},
	                 ObsTimes::AbstractArray{Tf},t0::Tf,dt::Vector{Tf},wave::Vector{Tf},sourceType::Symbol,
	                 storageLevel::Symbol,sensitivityMethod::Symbol,
	                 timeIntegrationMethod::Symbol,EMsolvers::Vector{TSem},
	                 DCsolver::TSdc,modUnits::Symbol,
                     K::SparseMatrixCSC{Tf,Tn}) = MaxwellTimeParam{Tf,Tn,TSem,TSdc}(
	                                                 Mesh,Sources,Obs,ObsTimes,t0,dt,wave,
	                                                 sourceType,storageLevel,
	                                                 sensitivityMethod,
                                                     timeIntegrationMethod,
                                                     EMsolvers,DCsolver,
                                                     modUnits,K)

# Supported options for categorical settings. See getMaxwellTimeParam
# docstring below for documentation.
supportedIntegrationMethods = [:BE; :BDF2; :BDF2Const; :TRBDF2]
supportedSourceTypes        = [:InductiveDiscreteWire;
                               :InductiveLoopPotential; :Galvanic]
supportedStorageLevels      = [:Factors; :Matrices; :None]
supportedSensitivityMethods = [:Implicit; :Explicit]
supportedModUnits           = [:con; :res]

# Conditionally set default solver, depending on whether or not user
# has MUMPS installed
defaultSolver = hasMUMPS ? :MUMPS : :juliaSolver
#defaultSolver = :Pardiso
"""
function param = getMaxwellTimeParam(Mesh,Sources,Obs,ObsTimes,dt,
                                     wave,sourceType;kwargs)

Input:  Mandatory arguments:

        Mesh::AbstractMesh
        Sources - Array of size ne X ns where ne is number of mesh edges
                  and ns number of sources. Each column contains integral
                  of source wire path approximated onto the mesh. Normally
                  computed using getEdgeIntegralofPolygonalChain.
        Obs -     Linear measurement matrix. Computes data at times
                  at which electric field is measured as Dtmp = Obs'\*efield(tn)
        ObsTimes - Vector of observation times. Data is first computed at
                   times with efields then linearly interpolated to observation
                   times.
        t0  - Real number. Initial time. Needed to link ObsTimes and dt.
        dt  - vector holding size of each time-step. For direct solvers,
              the number of factorizations = length(unique(dt))
        wave - vector of length 1+length(dt) giving the magnitude of the
               current at each time in [0;cumsum(dt)]
        sourceType::Symbol - Options are :Galvanic, :InductiveDiscreteWire,
                             and :InductiveLoopPotential.
                             Inductive sources assume zero electric field
                             at t=0 (initial time). Solve DC problem to
                             compute initial electric fields. DC electric
                             field data are computed for galvanic sources

        Optional keyword arguments:

        storageLevel::Symbol - Options are :Factors, :Matrices, :None.
                               Default value :Factors. If a direct solver
                               is being used and :Factors is chosen,
                               factorizations computed during forward
                               modelling will be stored for future use in
                               sensitivity computations. A new factorization
                               is computed for each unique step-size. the
                               :Matrices stores matrices but not
                               factorizations. It is mostly for iterative
                               solvers but may also be used to gain small
                               savings in computation time when using
                               on problems too large to store multiple
                               factorizations simultaneously.

        sensitivityMethod::Symbol - Options are :Implicit, :Explicit. Default
                               value is :Implicit.

        timeIntegrationMethod::Symbol - Options are :BE for backward Euler,
                               :BDF2 for variable step-size second order
                               backward differentiation formula, :BDF2Const
                               for BDF2 with constant step-size, :TRBDF2
                               for unstable trapezoidal BDF2 method. It's
                               only in here for experimentation. Perhaps
                               we'll add exponential integrators at some
                               point too.

        EMsolverType::Symbol - Linear solver for EM time-stepping solves.
                               Currently supported options
                               are :MUMPS, :Pardiso and :JuliaSolver (to use
                               Julia's built in sparse Cholesky routine,
                               which I believe is cholmod, from SuiteSparse).
                               Default is :MUMPS

        DCsolverType::Symbol - Solver for DC problem for galvanic sources.
                               Same options as EMSolverType.



"""
function getMaxwellTimeParam{S<:Real}(Mesh::AbstractMesh,
                                                 Sources::AbstractArray{S},
                                                 Obs::AbstractArray{S},
                                                 ObsTimes::Vector{S},
                                                 t0::S,
                                                 dt::Vector{S},
                                                 wave::Vector{S},
                                                 sourceType::Symbol;
			                         storageLevel::Symbol=:Factors,
			                         sensitivityMethod::Symbol=:Implicit,
			                         timeIntegrationMethod::Symbol=:BE,
			                         EMsolverType::Symbol=defaultSolver,
			                         DCsolverType::Symbol=defaultSolver,
                                     modUnits::Symbol=:con)

    # Check that user has chosen valid settings for categorical options
    in(timeIntegrationMethod,supportedIntegrationMethods) || error("Unsupported integration method")
    in(sourceType,supportedSourceTypes) || error("Unsupported source type")
    in(storageLevel,supportedStorageLevels) || error("Unknown storageLevel selection")
    in(sensitivityMethod,supportedSensitivityMethods) || error("Invalid sensitivity method")
    in(modUnits,supportedModUnits) || error("Invalid model unit selected")

    # Check consistency of dt and wave
    if length(wave) != (length(dt)+1)
        error("length(wave) != length(dt)+1")
    end

    # Setup solvers based on input options
    if DCsolverType == :MUMPS
        DCsolver = hasMUMPS ? getMUMPSsolver([],1,0,1) : error("Unable to load MUMPS")
    elseif DCsolverType ==:Pardiso
        DCsolver = hasPardiso ? getjInvPardisoSolver([],1,0,2,1) : error("Unable to load Pardiso")
        warn("DC problem can be unstable w Pardiso for large problems")
    elseif DCsolverType == :juliaSolver
        DCsolver = getJuliaSolver(sym=1)
    end

    if EMsolverType == :MUMPS
        baseEMSolver = hasMUMPS ? getMUMPSsolver([],1,0,1) : error("Unable to load MUMPS")
    elseif EMsolverType ==:Pardiso
        baseEMSolver = hasPardiso ? getjInvPardisoSolver([],1,0,2,1) : error("Unable to load Pardiso")
    elseif EMsolverType == :juliaSolver
        baseEMSolver = getJuliaSolver(sym=1)
    else
        error("Solver $(EMsolverType) not supported by MaxwellTime")
    end

    # For direct solvers we use one factorization per unique time-step size.
    # When storing the factors we need one jInv.LinearSolvers solver object
    # per factorization
    nFacs = length(unique(dt))
    if (storageLevel == :Factors) && (timeIntegrationMethod != :BDF2Const)
        EMsolvers = [copySolver(baseEMSolver) for i=1:nFacs]
    else
        EMsolvers = [baseEMSolver]
    end

    # The following is a no-op for conformal meshes.
    # For non-conformal meshes, restrict source to active edges.
    # This is done here to save having to do it every iteration
    # in an inversion.
    Ne,Qe, = getEdgeConstraints(Mesh)
    s      = Ne'*Sources

    # Get rid of late times that are greater than observed.
    maxobstime = maximum(ObsTimes)
    timevalues = t0 + cumsum(dt)
    maxidx = searchsortedfirst(timevalues, maxobstime) + 1
    deleteat!(wave, maxidx+1:length(wave))
    deleteat!(dt, maxidx:length(dt))


    # ObsTimeMat interpolates observations from step times
    # to observation times
    ObsTimeMat = getObsTimeMatrix(ObsTimes,t0,dt,size(Obs,2),size(s,2),sourceType)

    K = getMaxwellCurlCurlMatrix(Mesh,fill(pi*4e-7,Mesh.nc))
    #println("Is returned to getParam as type $(typeof(K))")

    if sourceType != :Galvanic
        clear!(Mesh.FX) ; clear!(Mesh.FY) ; clear!(Mesh.FZ)
        clear!(Mesh.EX) ; clear!(Mesh.EY) ; clear!(Mesh.EZ)
        clear!(Mesh.NFX); clear!(Mesh.NFY); clear!(Mesh.NFZ)
        clear!(Mesh.NEX); clear!(Mesh.NEY); clear!(Mesh.NEZ)
        clear!(Mesh.NN)
    end
    clear!(Mesh.NC)

    if sourceType == :InductiveLoopPotential
        s = -0.5*K*s
    end

    return MaxwellTimeParam(Mesh,s,Obs,ObsTimeMat,t0,dt,wave,sourceType,storageLevel,
                            sensitivityMethod,timeIntegrationMethod,
                            EMsolvers,DCsolver,modUnits,K)
end

function getObsTimeMatrix{S<:Real}(ObsTimes::Vector{S},t0::S,dt::Vector{S},
                                   nr::Integer,ns::Integer,sourceType::Symbol)
    t = sourceType == :Galvanic ? t0 + [0;cumsum(dt)] : t0 + cumsum(dt)
    #obsIdx = [searchsortedfirst(t,ti) for ti in ObsTimes]
    nt = length(ObsTimes)
    interpWeights = zeros(S,nt)
    mat = spzeros(S,nt*nr*ns,length(t)*nr*ns)
    for i in 1:nt
        obsIdx = searchsortedfirst(t,ObsTimes[i])
        if obsIdx > length(t) ; error("obsIdx > length(t)") ; end
        if obsIdx == 1
            interpWeights[i] = 1
            mat[1:nr*ns,1:nr*ns] = speye(nr*ns)
        else
            interpWeights[i] = (t[obsIdx]-ObsTimes[i])/(t[obsIdx]-t[obsIdx-1])
            w = interpWeights[i]
            mat[(i-1)*nr*ns+1:i*nr*ns,(obsIdx-2)*nr*ns+1:obsIdx*nr*ns] =
                spdiagm((w*ones(nr*ns),(1-w)*ones(nr*ns)),(0,nr*ns),nr*ns,2*nr*ns)
        end
    end
    return mat
end
