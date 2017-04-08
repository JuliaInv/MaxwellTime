export solveMaxTimeBE!, solveMaxTimeBDF2ConstDT!, solveDC!

typealias directSolver Union{MUMPSsolver,jInvPardisoSolver,JuliaSolver}

#Handles direct solvers for BE
function solveMaxTimeBE!{T<:directSolver}(A,rhs,Msig,
                      M::AbstractMesh,dt::Vector{Float64},it::Int64,iFacIn::Int64,
                      storageLevel::Symbol,linSolParamList::Array{T},flag=0)
#= 
Solve the maxwell system using MUMPS or Pardiso

en, = solveMaxTimeBE!{T<:directSolver}(A,rhs,Msig,
                    M::AbstractMesh,dt::Vector{Float64},it::Int64,iFacIn::Int64,
                    linSolParamList::Array{T},flag=0)

=#
    iFac = iFacIn
    if storageLevel != :Factors
        iFac = 1
        if ( (it==1) || (dt[it] != dt[it-1]) )
            linSolParamList[iFac].doClear = 1
        end
    end
    X                             = zeros(size(rhs))
    X, linSolParamList[iFac]      = solveLinearSystem!(A,rhs,X,linSolParamList[iFac],
                                                       flag)
    linSolParamList[iFac].doClear = 0

return X, linSolParamList
end

#Handles direct solvers for BDF-2
function solveMaxTimeBDF2ConstDT!(A,rhs,Msig,
                           M::AbstractMesh,dt::Real,
                           linSolParam::directSolver,flag=0)
#= 
 	Solve the maxwell system using MUMPS or Pardiso

 	en, = solveMaxTimeBDF2ConstDT!(A,rhs,Msig,M::AbstractMesh,w::Real,linSolParam::Union{MUMPSsolver,jInvPardisoSolver},flag=0)

=#
  X              = zeros(size(rhs))
  X, linSolParam = solveLinearSystem!(A,rhs,X,linSolParam,flag)
  

return X, linSolParam
end

function solveDC!(A,rhs,linSolParam::directSolver,flag=0)
#= 
 	Solve the maxwell system using MUMPS

 	en, = solveDC!(A,rhs,linSolParam::directSolver,flag=0)

=#
  X              = zeros(size(rhs))
  X, linSolParam = solveLinearSystem!(A,rhs,X,linSolParam,flag)
  

  return X, linSolParam
end

# Commented out functions are old, probably broken code to solve e-field system
# using PCG or block-PCG preconditioned with SSOR decomposition of A-phi system
# Kept here so someone wanting to resurrect the functionality in the future doesn't
# have to start from scratch

# function solveMaxTime(A,rhs,Msig,M::AbstractMesh,dt::Real,linSolParam::CGSolver,flag=0)
#   # setup preconditioner using Aphi system
#   mu = 4*pi*1e-7
#   Grad = getNodalGradientMatrix(M)
#   Ace  = getEdgeAverageMatrix(M)
#   Ane  = abs(Grad)
#   V    = getVolume(M)
#   v     = diag(M.V)
#   Ne,Qe = getEdgeConstraints(M)
#   Nn,   = getNodalConstraints(M)
#   Grad  = Qe*Grad*Nn
#   
#   muInvCells = Ane'*(Ace'*v)
#   Mmuinvn    = Nn'*sdiag(muInvCells)*Nn
#   
#   STBa = Grad*Mmuinvn*Grad'
#   Aap = [A + STBa          1/dt*Msig*Grad; 
#   		1/dt*Grad'*Msig    1/dt*Grad'*Msig*Grad];
#   
#   Map(x) = ssor(Aap,x,out=-1)[1]
#   
#   P  = [speye(size(A,2)); Grad']
#   MM(x) = P'*Map(P*x)
#   en, = solveLinearSystem(A,rhs,linSolParam,flag)
#   return en, linSolParam
# end

# function solveMaxTime!(A,rhs,Msig,M::AbstractMesh,dt::Real,linSolParam::blockPCGSolver,flag=0)
#   # setup preconditioner using Aphi system
#   mu = 4*pi*1e-7
#   Grad = getNodalGradientMatrix(M)
#   Ace  = getEdgeAverageMatrix(M)
#   Ane  = abs(Grad)
#   V    = getVolume(M)
#   v     = diag(M.V)
#   Ne,Qe = getEdgeConstraints(M)
#   Nn,   = getNodalConstraints(M)
#   Grad  = Qe*Grad*Nn
#   
#   muInvCells = Ane'*(Ace'*v)
#   Mmuinvn    = Nn'*sdiag(muInvCells)*Nn
#   
#   STBa = Grad*Mmuinvn*Grad'
#   Aap = [A + STBa          1/dt*Msig*Grad; 
#   		1/dt*Grad'*Msig    1/dt*Grad'*Msig*Grad];
#   
#   Map(x) = ssor(Aap,x,out=-1)[1]
#   
#   P  = [speye(size(A,2)); Grad']
#   MM(x) = P'*Map(P*x)
#   en, = solveLinearSystem!(A,rhs,linSolParam,flag)
#   return en, linSolParam
# end


