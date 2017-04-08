export getData

# For backward compatibility and convenience if first argument to getData
# is a vector of real numbers and not an object of type MaxwellTimeModel
# then assume first argument is conductivity and set permeability mu=mu0
function getData{S<:Real}(sigma::Vector{S},param::MaxwellTimeParam)
    m = MaxwellTimeModel(sigma,fill(4*pi*1e-7,param.M.nc))
    return getData(m,param)
end


"""
function getData(model::MaxwellTimeModel,param::MaxwellTimeParam)

Given electrical conductivity and magnetic permeability models, compute the
electric field at requested time steps, store the fields in param, and
return data and param.

Input:

      model::MaxwellTimeModel--Composite type containing conductivity 
                               and permeability. Conductivity can be isotropic,
                               diagonally anisotropic or generally anisotropic.
                               Permeability may be isotropic or diagonally
                               anisotropic.
      param::MaxwellTimeParam--Contains structures and information needed to 
                               construct and solve Maxwell's equations such
                               as a mesh, a source, and time steps. For a full
                               description see the MaxwellTimeParam documentation.
                               
      
      Output: 
      
             D::Array--Data, computed by applying the observation matrix
                       to the electric fields at each time step. For more
                       information on the observation matrix see the
                       MaxwellTimeParam documentation.
                       
             param::MaxwellTimeParam--param is returned modified. It may store
                    matrix factorizations and the electric fields for further
                    use in sensitivity computations.
"""
function getData(model::MaxwellTimeModel,param::MaxwellTimeParam)
#function getData(sigma,param)   
        
    # Unpack model into conductivity and magnetic permeability
    sigma = model.sigma
    mu    = model.mu
    
    #Unpack param
    M             = param.M
    sourceType    = param.sourceType
    storageLevel  = param.storageLevel
    EMsolver      = param.EMsolvers   
    dt            = param.dt
    
    # Check model input 
    in(length(sigma),[M.nc; 3*M.nc; 6*M.nc]) || error("MaxwellTime.getData: Invalid length of sigma")
    if length(mu) == 6*M.nc
        error("MaxwellTime.getData: Generally anisotropic permeability not supported")
    elseif ~in(length(mu),[M.nc; 3*M.nc])
        error("MaxwellTime.getData: Invalid length of mu")
    end

    # If explicit sensitivities are being used, clear them
    if param.sensitivityMethod == :Explicit
        param.Sens = Array{eltype(sigma),2}() #Empty 2D array
    end
    
    # Form conductivity edge mass matrix
    Ne,Qe, = getEdgeConstraints(M)
    Msig   = getEdgeMassMatrix(M,sigma)
    Msig   = Ne'*Msig*Ne
    
    # Form K = Curl'*Mmu*Curl
    K = getMaxwellCurlCurlMatrix(M,mu)
    
    # Restrict source to active (not hanging) edges
    s = Ne'*param.Sources
    
    # Initialize electric field storage
    ne            = size(Ne,2) #Number of active edges
    ns            = size(s,2)
    nt            = length(dt)+1
    param.fields  = zeros(eltype(s),ne,ns,nt)
    e             = param.fields
    
    # Initialize matrix storage if needed. Note that if a direct solver is
    # being used and factorizatons are being stored, then matrices don't 
    # need to be stored
    if storageLevel == :Matrices
        param.Matrices = Vector{SparseMatrixCSC{eltype(G.nzval),eltype(G.colptr)}}()
        Matrices       = param.Matrices
    end
    
    # Get initial fields. They're zero for inductive sources
    if sourceType == :Galvanic
        param = getFieldsDC(Msig,s,param)
    end
        
    # Compute the electric fields at times [cumsum(dt)]
    # param.timeIntegrationMethod is a symbol. See supportedIntegrationMethods
    # in MaxwellTime.jl for valid options.
    # integrationFunctions is a dictionary that maps timeIntegrationMethod
    # symbols to their corresponding getTransientFields<IntegrationMethod> 
    # functions.
    getTransientFields  = integrationFunctions[param.timeIntegrationMethod]
    param               = getTransientFields(K,Msig,s,param)
    
    
    # Compute the data from the fields.
    # Note (for grounded sources) that DC data
    # Are computed as the integral of electric field over a receiver
    # dipole and not as a potential difference
    Pt = param.Obs'*Ne
    if sourceType == :Galvanic
        D = zeros(size(Pt,1),ns,length(dt)+1)
        D[:,:,1] = Pt*e[:,:,1]
        offset = 1
    else
        D = zeros(size(Pt,1),ns,length(dt))
        offset = 0
    end
    for i=1:length(dt)
        D[:,:,i+offset] = Pt*e[:,:,i+1]
    end
    
    return D, param
end
