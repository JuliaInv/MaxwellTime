export solveMaxTime!

# function solveMaxTime(A,rhs,Msig,M::AbstractMesh,dt::Real,linSolParam::iterativeSolver,flag=0)
# 	if linSolParam.doClear == 1
# 		iterTol = 1e-10
# 		mu = 4*pi*1e-7
# 		# setup preconditioner using Aphi system
# 
# 		Grad = getNodalGradientMatrix(M)
# 		Ace  = getEdgeAverageMatrix(M)
# 		Ane  = abs(Grad)
# 		V    = getVolume(M)
# 		v    = diag(M.V)
# 
# 		#println(size(Ane'),"   ",size(Ace'),"    ",size(v))
# 		muInvCells = Ane'*(Ace'*v)
# 		Mmuinvn    = sdiag(muInvCells)
# 
# 		STBa = Grad*Mmuinvn*Grad'
# 		Aap = [A + STBa          1/dt*Msig*Grad; 
# 				1/dt*Grad'*Msig    1/dt*Grad'*Msig*Grad];
#  
# 		Map(x) = sor(Aap,x,out=-1)[1]
# 
# 		P  = [speye(size(A,2)); Grad']
# 		MM(x,v) = P'*Map(P*x)
# 		linSolParam.Ainv = MM
# 	end
# 	en, = solveLinearSystem(A,rhs,linSolParam,flag)
# 
# 	return en, linSolParam
# end

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

function solveMaxTime!(A,rhs,Msig,M::AbstractMesh,dt::Real,linSolParam::Union{MUMPSsolver,jInvPardisoSolver},flag=0)
#= 
 	Solve the maxwell system using MUMPS

 	en, = solveMaxTime!(A,rhs,Msig,M::AbstractMesh,w::Real,linSolParam::Union{MUMPSsolver,jInvPardisoSolver},flag=0)

=#
  X              = zeros(size(rhs))
  X, linSolParam = solveLinearSystem!(A,rhs,X,linSolParam,flag)
  

return X, linSolParam
end