export MaxwellTimeParam, getMaxwellTimeParam

# MaxwellTime forward problem type. Symbols are used to select categorical options.
# Inner and outer constructors are required to allow instantiating the type with
# last three fields left uninitialized. Use getMaxwellTimeParam function
# for instantiation. It validates input and sets things to default values.
type MaxwellTimeParam{S<:Real,T<:AbstractSolver,U<:AbstractSolver} <: ForwardProbType
    Mesh::AbstractMesh
    Sources::AbstractArray{S}
    Obs::AbstractArray{S}
    ObsTimes::AbstractArray{S}
    dt::Vector{S}
    wave::Vector{S}
    sourceType::Symbol
    storageLevel::Symbol
    sensitivityMethod::Symbol
    timeIntegrationMethod::Symbol
    EMsolvers::Vector{T}
    DCsolver::U
    K::SparseMatrixCSC
    Matrices::Vector{SparseMatrixCSC}
    fields::Array{S}
    Sens::Array{S}
    AuxFields::Array{S} # Stores extra electric fields outside
                        # of regular time-stepping. Currently used
                        # to initialize BDF2 time stepping without additional BE
                        # factorization. Needs to be stored for use in sensitivity
                        # computations

    # Use custom constructors to allow creation of objects with
    # Matrices, fields, explicit sensitivity matrix, and possibly K
    # left uninitialized
    MaxwellTimeParam(Mesh::AbstractMesh,Sources::AbstractArray{S},Obs::AbstractArray{S},
                     ObsTimes::AbstractArray{S},dt::Vector{S},wave::Vector{S},sourceType::Symbol,
                     storageLevel::Symbol,sensitivityMethod::Symbol,
                     timeIntegrationMethod::Symbol,EMsolvers::Vector{T},
                     DCsolver::U) = new(Mesh,Sources,Obs,ObsTimes,dt,wave,sourceType,
                                     storageLevel,sensitivityMethod,
                                     timeIntegrationMethod,
                                     EMsolvers,DCsolver)

    MaxwellTimeParam(Mesh::AbstractMesh,Sources::AbstractArray{S},Obs::AbstractArray{S},
                     ObsTimes::AbstractArray{S},dt::Vector{S},wave::Vector{S},sourceType::Symbol,
                     storageLevel::Symbol,sensitivityMethod::Symbol,
                     timeIntegrationMethod::Symbol,EMsolvers::Vector{T},
                     DCsolver::U,K::SparseMatrixCSC) = new(
                                     Mesh,Sources,Obs,ObsTimes,dt,wave,
                                     sourceType,storageLevel,sensitivityMethod,
                                     timeIntegrationMethod,
                                     EMsolvers,DCsolver,K)
end

# Unfortunately parametric types need these matching outer and inner constructors
MaxwellTimeParam{S,T,U}(Mesh::AbstractMesh,Sources::AbstractArray{S},Obs::AbstractArray{S},
	                 ObsTimes::AbstractArray{S},dt::Vector{S},wave::Vector{S},sourceType::Symbol,
	                 storageLevel::Symbol,sensitivityMethod::Symbol,
	                 timeIntegrationMethod::Symbol,EMsolvers::Vector{T},
	                 DCsolver::U) = MaxwellTimeParam{S,T,U}(
	                                                 Mesh,Sources,Obs,ObsTimes,dt,wave,
	                                                 sourceType,storageLevel,
	                                                 sensitivityMethod,
                                                     timeIntegrationMethod,
                                                     EMsolvers,DCsolver)

MaxwellTimeParam{S,T,U}(Mesh::AbstractMesh,Sources::AbstractArray{S},Obs::AbstractArray{S},
	                 ObsTimes::AbstractArray{S},dt::Vector{S},wave::Vector{S},sourceType::Symbol,
	                 storageLevel::Symbol,sensitivityMethod::Symbol,
	                 timeIntegrationMethod::Symbol,EMsolvers::Vector{T},
	                 DCsolver::U,K::SparseMatrixCSC) = MaxwellTimeParam{S,T,U}(
	                                                 Mesh,Sources,Obs,ObsTimes,dt,wave,
	                                                 sourceType,storageLevel,
	                                                 sensitivityMethod,
                                                     timeIntegrationMethod,
                                                     EMsolvers,DCsolver,K)

# Supported options for categorical settings. See getMaxwellTimeParam
# docstring below for documentation.
supportedIntegrationMethods = [:BE; :BDF2; :BDF2Const; :TRBDF2]
supportedSourceTypes        = [:InductiveDiscreteWire;
                               :InductiveAnalyticLoop; :Galvanic]
supportedStorageLevels      = [:Factors; :Matrices; :None]
supportedSensitivityMethods = [:Implicit; :Explicit]

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
        dt  - vector holding size of each time-step. For direct solvers,
              the number of factorizations = length(unique(dt))
        wave - vector of length 1+length(dt) giving the magnitude of the
               current at each time in [0;cumsum(dt)]
        sourceType::Symbol - Options are :Galvanic, :InductiveDiscreteWire,
                             and :InductiveAnalyticLoop.
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

        EMSolverType::Symbol - Linear solver for EM time-stepping solves.
                               Currently supported options
                               are :MUMPS, :Pardiso and :JuliaSolver (to use
                               Julia's built in sparse Cholesky routine,
                               which I believe is cholmod, from SuiteSparse).
                               Default is :MUMPS

        DCSolverType::Symbol - Solver for DC problem for galvanic sources.
                               Same options as EMSolverType.



"""
function getMaxwellTimeParam{S<:Real}(Mesh::AbstractMesh,
                                                 Sources::AbstractArray{S},
                                                 Obs::AbstractArray{S},
                                                 ObsTimes::Vector{S},
                                                 dt::Vector{S},
                                                 wave::Vector{S},
                                                 sourceType::Symbol;
			                         storageLevel::Symbol=:Factors,
			                         sensitivityMethod::Symbol=:Implicit,
			                         timeIntegrationMethod::Symbol=:BE,
			                         EMsolverType::Symbol=:MUMPS,
			                         DCsolverType::Symbol=:MUMPS)

    # Check that user has chosen valid settings for categorical options
    in(timeIntegrationMethod,supportedIntegrationMethods) || error("Unsupported integration method")
    in(sourceType,supportedSourceTypes) || error("Unsupported source type")
    in(storageLevel,supportedStorageLevels) || error("Unknown storageLevel selection")
    in(sensitivityMethod,supportedSensitivityMethods) || error("Invalid sensitivity method")

    # Check consistency of dt and wave
    if length(wave) != (length(dt)+1)
        error("length(wave) != length(dt)+1")
    end

    # Setup solvers based on input options
    if DCsolverType == :MUMPS
        DCsolver = getMUMPSsolver([],1,0,1)
    elseif DCsolverType ==:Pardiso
        error("Not set up yet")
        #warn("DC problem can be unstable w Pardiso for large problems")
    elseif DCsolverType == :juliaSolver
        DCsolver = getJuliaSolver(sym=1)
    end

    if EMsolverType == :MUMPS
        baseEMSolver = getMUMPSsolver([],1,0,1)
    elseif EMsolverType ==:Pardiso
        error("Pardiso solver will be supported soon!")
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

    # ObsTimeMat interpolates observations from step times
    # to observation times
    ObsTimeMat = getObsTimeMatrix(ObsTimes,dt,size(Obs,2),size(s,2),sourceType)

    if sourceType == :InductiveAnalyticLoop
        K = getMaxwellCurlCurlMatrix(Mesh,fill(pi*4e-7,Mesh.nc))
        s = -K*s
        return MaxwellTimeParam(Mesh,s,Obs,ObsTimeMat,dt,wave,sourceType,storageLevel,
                                sensitivityMethod,timeIntegrationMethod,
                                EMsolvers,DCsolver,K)
    else
        K = spzeros(0,0)
        return MaxwellTimeParam(Mesh,s,Obs,ObsTimeMat,dt,wave,sourceType,storageLevel,
                                sensitivityMethod,timeIntegrationMethod,
                                EMsolvers,DCsolver,K)
    end
end

function getObsTimeMatrix{S<:Real}(ObsTimes::Vector{S},dt::Vector{S},
                                   nr::Integer,ns::Integer,sourceType::Symbol)
    t = sourceType == :Galvanic ? [0;cumsum(dt)] : cumsum(dt)
    #obsIdx = [searchsortedfirst(t,ti) for ti in ObsTimes]
    nt = length(ObsTimes)
    interpWeights = zeros(S,nt)
    mat = spzeros(S,nt*nr*ns,length(t)*nr*ns)
    for i in 1:nt
        obsIdx = searchsortedfirst(t,ObsTimes[i])
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
