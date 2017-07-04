export solveMaxTimeBE!, solveMaxTimeBDF2ConstDT!, solveDC!

const directSolver = Union{MUMPSsolver,jInvPardisoSolver,JuliaSolver}

"""
function solveMaxTimeBE!(A,rhs,Msig,M::AbstractMesh,dt::Vector{Real},it::Integer,
                      storageLevel::Symbol,linSolParam::directSolver,flag=0)

Solve the maxwell system using selected direct solver
"""
function solveMaxTimeBE!{T<:Real}(A,rhs,Msig,M::AbstractMesh,dt::Vector{T},it::Integer,
                      storageLevel::Symbol,linSolParam::directSolver,flag=0)

    X                   = zeros(eltype(rhs),size(rhs))
    X, linSolParam      = solveLinearSystem!(A,rhs,X,linSolParam,flag)
    linSolParam.doClear = 0

return X, linSolParam
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

"""
en, = solveDC!(A,rhs,linSolParam::directSolver,flag=0)

Solve the divSigmaGrad system using selected direct solver
"""
function solveDC!(A,rhs,linSolParam::directSolver,flag=0)

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
