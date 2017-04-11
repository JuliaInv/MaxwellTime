module MaxwellTime
	
using jInv.Mesh.AbstractMesh
using JOcTree
using jInv.Utils
using jInv.LinearSolvers
using KrylovMethods
		
export MaxwellTimeParam, getMaxwellTimeParam
import jInv.ForwardShare.ForwardProbType
export ForwardProbType

# Supported options for various settings
supportedIntegrationMethods = [:BE; :BDF2; :BDF2Const; :TRBDF2]
supportedSourceTypes        = [:Inductive; :Galvanic]
supportedStorageLevels      = [:Factors; :Matrices; :None]
supportedSensitivityMethods = [:Implicit; :Explicit]

# MaxwellTime forward problem type. Symbols are used to select categorical options.
# Inner and outer constructors are required to allow instantiating the type with
# last three fields left uninitialized. Use getMaxwellTimeParam function 
# for instantiation. It validates input and sets default values.
type MaxwellTimeParam{S<:Real,T<:AbstractSolver,U<:AbstractSolver} <: ForwardProbType
    M::AbstractMesh
    Sources::AbstractArray{S}
    Obs::AbstractArray{S}
    dt::Vector{S}
    wave::Vector{S}
    sourceType::Symbol
    storageLevel::Symbol
    sensitivityMethod::Symbol
    timeIntegrationMethod::Symbol
    EMsolvers::Vector{T}
    DCsolver::U
    Matrices::Vector{SparseMatrixCSC}
    fields::Array{S}
    Sens::Array{S}
    AuxFields::Array{S} # Stores extra electric fields outside
                        # of regular time-stepping. Currently used
                        # to initialize BDF2 time stepping without additional BE
                        # factorization. Needs to be stored for use in sensitivity
                        # computations
    
    # Use inner constructor to allow creation of objects with
    # Matrices, fields, and explicit sensitivity matrix left uninitialized
    MaxwellTimeParam(M::AbstractMesh,Sources::AbstractArray{S},Obs::AbstractArray{S},
                     dt::Vector{S},wave::Vector{S},sourceType::Symbol,
                     storageLevel::Symbol,sensitivityMethod::Symbol,
                     timeIntegrationMethod::Symbol,EMsolvers::Vector{T},
                     DCsolver::U) = new(M,Sources,Obs,dt,wave,sourceType,
                                     storageLevel,sensitivityMethod,
                                     timeIntegrationMethod,
                                     EMsolvers,DCsolver)
end
MaxwellTimeParam{S,T,U}(M::AbstractMesh,Sources::AbstractArray{S},Obs::AbstractArray{S},
	                 dt::Vector{S},wave::Vector{S},sourceType::Symbol,
	                 storageLevel::Symbol,sensitivityMethod::Symbol,
	                 timeIntegrationMethod::Symbol,EMsolvers::Vector{T},
	                 DCsolver::U) = MaxwellTimeParam{S,T,U}(
	                                                 M,Sources,Obs,dt,wave,
	                                                 sourceType,storageLevel,
	                                                 sensitivityMethod,
                                                         timeIntegrationMethod,
                                                         EMsolvers,DCsolver)
"""
function param = getMaxwellTimeParam(M,Sources,Obs,dt,wave,sourceType;kwargs)

Input:  Mandatory arguments:

        M::AbstractMesh
        Sources - Array of size ne X ns where ne is number of mesh edges 
                  and ns number of sources. Each column contains integral
                  of source wire path approximated onto the mesh. Normally
                  computed using getEdgeIntegralofPolygonalChain.
        Obs - Linear measurement matrix. Data at time tn = Obs'*efield(tn)
        dt  - vector holding size of each time-step. For direct solvers,
              the number of factorizations = length(unique(dt))
        wave - vector of length 1+length(dt) giving the magnitude of the
               current at each time in [0;cumsum(dt)]
        sourceType::Symbol - Options are :Inductive and :Galvanic.
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
function getMaxwellTimeParam{S<:Real,T<:Integer}(M::AbstractMesh,
                                                 Sources::AbstractArray{S},
                                                 Obs::SparseMatrixCSC{S,T},
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
    in(sourceType,supportedSourceTypes) || error("Source type must be :Inductive or :Galvanic")
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
    
    #Factor when there's a new time step size
    nFacs = length(unique(dt))
    if (storageLevel == :Factors) && (timeIntegrationMethod != :BDF2Const)
        EMsolvers = [copySolver(baseEMSolver) for i=1:nFacs]
    else
        EMsolvers = [baseEMSolver] 
    end
    
    # The following is a no-op for conformal meshes
    # For non-conformal meshes, restrict source to active edges.
    # This is done here to save having to do it every iteration
    # in an inversion.
    Ne,Qe, = getEdgeConstraints(M)
    s      = Ne'*Sources
    
    return MaxwellTimeParam(M,s,Obs,dt,wave,sourceType,storageLevel,
                            sensitivityMethod,timeIntegrationMethod,
                            EMsolvers,DCsolver)
end


# Define earth model type
immutable MaxwellTimeModel{S<:Real}
    sigma::Vector{S}
    mu::Vector{S}
    invertSigma::Bool
    invertmu::Bool  
    
    MaxwellTimeModel(sigma::Vector{S},mu::Vector{S}) = 
      new(sigma,mu,true,false)
end

MaxwellTimeModel{S}(sigma::Vector{S},mu::Vector{S}) = 
  MaxwellTimeModel{S}(sigma,mu)

# Define time-stepping functions and map integration method symbols to
# the appropriate functions defined in getFields.jl
include("getFields.jl")
integrationFunctions = Dict(zip(supportedIntegrationMethods,
                                [getFieldsBE;getFieldsBDF2;
                                 getFieldsBDF2ConstDT;
                                 getFieldsTRBDF2]))

# Add MaxwellTime specific methods to following three jInv generic functions
import jInv.ForwardShare.getData
import jInv.ForwardShare.getSensMatVec
import jInv.ForwardShare.getSensTMatVec
                                
include("getData.jl")
include("getSensMatVec.jl")
include("getSensTMatVec.jl")

# Interface MaxwellTime to jInv.LinearSolvers to solve linear systems
# of equations
include("solverFunctions.jl")

end #End module MaxwellTime
