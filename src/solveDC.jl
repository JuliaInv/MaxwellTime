export solveDC!

function solveDC!(A,rhs,linSolParam::Union{MUMPSsolver,jInvPardisoSolver},flag=0)
#= 
 	Solve the maxwell system using MUMPS

 	en, = solveDC!(A,rhs,linSolParam::Union{MUMPSsolver,jInvPardisoSolver},flag=0)

=#
  X              = zeros(size(rhs))
  X, linSolParam = solveLinearSystem!(A,rhs,X,linSolParam,flag)
  

  return X, linSolParam
end
