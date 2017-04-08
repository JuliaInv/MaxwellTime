# Compute electric fields in time using various time-stepping strategies
# All methods assume steady state (including zero-field) initial conditions

# Available time stepping methods are
# 1) Backward Euler (timeIntegrationMethod=:BE)
# 2) BDF-2 (timeIntegrationMethod=:BDF2)
# 3) TRBDF-2 (timeIntegrationMethod=:TRBDF2)
#    TRBDF-2 seems to be unstable. Only for experimentation.

function getFieldsBE{T,N}(K::SparseMatrixCSC{T,N},
                          Msig::SparseMatrixCSC{T,N},
                          s::AbstractArray{T},
                          param::MaxwellTimeParam)
# Solve the Backward Euler system
# Curl*e + b_t = 0
# Curl'*(Mmuinv*b) - Msig*e = s(t)
#
# By eliminating b and solving for e
# b_t =  -Curl*e
# (Curl'*Mmuinv*Curl + 1/dt*Msig)*e_{n+1} = - 1/dt*(Msig*e_{n} + s_n - s_{n+1})

    # Unpack param
    storageLevel  = param.storageLevel
    EMsolvers     = param.EMsolvers   
    dt            = param.dt
    wave          = param.wave
    M             = param.M
    
    # Clear any factorizations still being stored
    for solver in EMsolvers
      solver.doClear = 1
    end
    
    # Do the time-stepping
    iStepSize = 0
    A         = []
    ew        = param.fields
    for i=1:length(dt)
        dtinv = 1.0/dt[i]
        rhs = dtinv*(Msig*ew[:,:,i]+(wave[i]-wave[i+1])*s)
        # Matrix only changes when the step-size changes
        if ( (i==1) || (dt[i] != dt[i-1]) )
            A = K + dtinv*Msig
            iStepSize += 1
            if storageLevel == :Matrices
                push!(Matrices,A)
            end
        end
        # Solve the e-field update system. Msig and M left as inputs
        # in case iterative solver is used in the future. Currently, this code
        # only supports direct solvers
        ew[:,:,i+1],EMsolvers = solveMaxTimeBE!(A,rhs,Msig,M,dt,i,iStepSize,
                                               storageLevel,EMsolvers)
    end
    
    return param
end

#------------------------------------------------------------------------

function getFieldsBDF2{T,N}(K::SparseMatrixCSC{T,N},
                            Msig::SparseMatrixCSC{T,N},
                            s::AbstractArray{T},
                            param::MaxwellTimeParam)

    error("BDF2 with variable step-size not implemented.\nUse timeIntegrationMethod :BDF2Const to use constant stepsize")
                            
    # Unpack param
    storageLevel = param.storageLevel
    EMsolver     = param.EMsolvers[1]
    dt           = param.dt[1]
    nt           = length(param.dt)
    wave         = param.wave
    M            = param.M
    ew           = param.fields
    
    # BDF2 with step size dt has same matrix as backward Euler with
    # stepsize 2dt/3. To save a factorization (we're mostly concerned
    # with direct solvers here) we initialize the time-stepping by taking
    # two BE steps to get to t=4dt/3 and then we use element-wise linear 
    # interpolation to compute electric field at t=dt
    EMsolver.doClear = 1
    A                = K + 3/(2*dt)*Msig
    rhs              = 3/(2*dt)*( Msig*ew[:,:,1] + (wave[1]-wave[2])*s )
    ehat,EMsolver    = solveMaxTime!(A,rhs,Msig,Msh,2/(3*dt),EMsolver)
    EMsolver.doClear = 0
    rhs              = 3/(2*dt)*( Msig*ehat )
    ehat2,EMsolver   = solveMaxTime!(A,rhs,Msig,Msh,2/(3*dt),EMsolver)
    ew[:,:,2]        = 0.5*(ehat+ehat2)
    param.AuxFields  = ehat
    
    #Continue time-stepping using BDF2 with constant step-size
    for i=2:nt
      rhs = -3/(2*dt)*( (wave[i+1]-(4/3)*wave[i]+wave[i-1]/3)*s + 
                        Msig*(-4/3*ew[:,:,i] + ew[:,:,i-1]/3) )
      ew[:,:,i+1],EMsolver = solveMaxTime!(A,rhs,Msig,M,dt,EMsolver)
    end
    if storageLevel ~= :Factors
      clear!(EMsolver)
      EMsolver.doClear = 1
    end   
    if storageLevel == :Matrices
        push!(param.Matrices,A)
    end
    return param
end


#------------------------------------------------------------------------

function getFieldsBDF2ConstDT{T,N}(K::SparseMatrixCSC{T,N},
                                   Msig::SparseMatrixCSC{T,N},
                                   s::AbstractArray{T},
                                   param::MaxwellTimeParam)

    # Unpack param
    storageLevel = param.storageLevel
    EMsolver     = param.EMsolvers[1]
    dt           = param.dt[1]
    nt           = length(param.dt)
    wave         = param.wave
    M            = param.M
    ew           = param.fields
    
    # BDF2 with step size dt has same matrix as backward Euler with
    # stepsize 2dt/3. To save a factorization (we're mostly concerned
    # with direct solvers here) we initialize the time-stepping by taking
    # two BE steps to get to t=4dt/3 and then we use element-wise linear 
    # interpolation to compute electric field at t=dt
    EMsolver.doClear = 1
    A                = K + 3/(2*dt)*Msig
    rhs              = 3/(2*dt)*( Msig*ew[:,:,1] + (wave[1]-wave[2])*s )
    ehat,EMsolver    = solveMaxTimeBDF2ConstDT!(A,rhs,Msig,M,2/(3*dt),EMsolver)
    EMsolver.doClear = 0
    rhs              = 3/(2*dt)*( Msig*ehat )
    ehat2,EMsolver   = solveMaxTimeBDF2ConstDT!(A,rhs,Msig,M,2/(3*dt),EMsolver)
    ew[:,:,2]        = 0.5*(ehat+ehat2)
    param.AuxFields  = ehat
    
    #Continue time-stepping using BDF2 with constant step-size
    for i=2:nt
      rhs = -3/(2*dt)*( (wave[i+1]-(4/3)*wave[i]+wave[i-1]/3)*s + 
                        Msig*(-4/3*ew[:,:,i] + ew[:,:,i-1]/3) )
      ew[:,:,i+1],EMsolver = solveMaxTimeBDF2ConstDT!(A,rhs,Msig,M,dt,EMsolver)
    end
    if ~(storageLevel == :Factors)
      clear!(EMsolver)
      EMsolver.doClear = 1
    end   
    if storageLevel == :Matrices
        push!(param.Matrices,A)
    end
    return param
end

#------------------------------------------------------------------------

# A curiosity that seems to be unstable. May tinker with it at some point.
function getFieldsTRBDF2{S<:Real}(sigma::Vector{S},mu::Vector{S},param::MaxwellTimeParam)
    
    # Unpack param initialize electric field, and setup source
    M             = param.M
    sourceType    = param.sourceType
    e             = param.fields
    EMsolver      = param.EMsolvers    
    Matrices      = param.Matrices
    
    # Form conductivity edge mass matrix and 
    Ne,Qe,        = getEdgeConstraints(M)
    Msig          = getEdgeMassMatrix(M,sigma)
    Msig          = Ne'*Msig*Ne
    
    # Restrict source to active (not hanging) edges
    s             = Ne'*param.Sources
    
    # Initialize electric field
    ne            = size(Ne,2) #Number of active edges
    ns            = size(s,2)
    nt            = length(dt)+1
    e             = zeros(eltype(s),ne,ns,nt)
    
    # Initialize matrix storage if needed. Note that matrices aren't reused and
    # thus don't need to be stored if factorizations are being stored.
    if storageLevel == :Matrices
        Matrices = Vector{SparseMatrixCSC{eltype(G.nzval),eltype(G.colptr)}}()
    end
    
    # Get initial fields. They're zero for inductive sources
    if sourceType == :Galvanic
        param = getFieldsDC(Msig,param)
    end
    
    # Form K = Curl'*Mmu*Curl
    K = getMaxwellCurlCurlMatrix(M,mu)
    
    # Clear any factorizations still being stored
    for solver in EMsolvers
      solver.doClear = 1
    end
    
    #Setup gamma related stuff
    gma     = 2-sqrt(2)
    gmaFctr = 2/(gma*dt)
    f1      = 1/(gma*(2-gma))
    f2      = ((1-gma)^2)/(gma*(2-gma))
    
    # Do the time-stepping
    for i=1:nt
        #Trapezoidal rule to get e_(n+gma)
        if i > 1
            rhs = -K*e[:,:,i] + gmaFctr*( Msig*e[:,:,i] + (wave[i]-wave[i+1])*s )
        else
            rhs = gmaFctr*( Msig*e[:,:,i] + (wave[i]-wave[i+1])*s )
        end
        eGma,mySolver = solveMaxTime(A,rhs,Msig,Msh,dt,mySolver)
        mySolver.doClear = 0
        
        #BDF-2 to get e_(n+1)
        wvGma = 0.0
        rhs = gmaFctr*( f1*Msig*eGma - f2*Msig*e[:,:,i] + (wvGma-wave[i]-wave[i+1])*s )
        e[:,:,i+1],mySolver = solveMaxTimeBDF2!(A,rhs,Msig,Msh,dt,mySolver)
    end
    mySolver.doClear = 1
    
    return param
end

#And now stuff that's called by the above getFields methods

#-----------------------------------------------------------------

function getFieldsDC{T,N}(Msig::SparseMatrixCSC{T,N},
                          s::AbstractArray{T},
                          param::MaxwellTimeParam)
                          
    M            = param.M
    solver       = param.DCsolver
    storageLevel = param.storageLevel
    ew           = param.fields
    wave         = param.wave
    
    solver.doClear = 1
    
    G      = getNodalGradientMatrix(M)
    Nn,Qn, = getNodalConstraints(M)
    Ne,Qe, = getEdgeConstraints(M)
    G      = Qe*G*Nn 
    Adc    = G'*Msig*G
    
    phi0,solver = solveDC!(Adc,wave[1]*G'*s,solver)
    ew[:,:,1]   = -G*phi0
    if storageLevel == :Factors
        solver.doClear = 0
    elseif storageLevel == :Matrices
        clear!(solver)
        push!(param.Matrices,Adc)
    else
        clear!(solver)
    end
    return param
end

function getMaxwellCurlCurlMatrix(M::AbstractMesh,mu)	
    Curl   = getCurlMatrix(M)
    Mmu    = getFaceMassMatrix(M,1./mu)
    Nf,Qf, = getFaceConstraints(M)
    Ne,    = getEdgeConstraints(M)
    Curl   = Qf*Curl*Ne
    Mmu    = Nf'*Mmu*Nf
    K      = Curl'*Mmu*Curl
    return K
end