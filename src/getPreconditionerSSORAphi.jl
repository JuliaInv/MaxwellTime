export SSORAphi,MaxwellTimePrecon

type MaxwellTimePreconParam <: preconParam
  sigma
  timeStepIdx::Integer
  pForParam::MaxwellTimeParam
end

function SSORAphi(A::SparseMatrixCSC,paramPrecon::MaxwellTimePreconParam)

  M     = paramPrecon.pForParam.M
  it    = paramPrecon.timeStepIdx
  sigma = paramPrecon.sigma
  dt    = paramPrecon.pForParam.dt[it]
  
  # setup preconditioner using Aphi system
  mu = 4*pi*1e-7
  Grad = getNodalGradientMatrix(M)
  Ace  = getEdgeAverageMatrix(M)
  Ane  = abs(Grad)
  V    = getVolume(M)
  v     = diag(M.V)
  Ne,Qe = getEdgeConstraints(M)
  Nn,   = getNodalConstraints(M)
  Grad  = Qe*Grad*Nn
  Msig  = Ne'*Msig*Ne
  
  muInvCells = Ane'*(Ace'*v)
  Mmuinvn    = Nn'*sdiag(muInvCells)*Nn
  
  STBa = Grad*Mmuinvn*Grad'
  Aap = [A + STBa          1/dt*Msig*Grad; 
  		1/dt*Grad'*Msig    1/dt*Grad'*Msig*Grad];
  
  Map(x) = ssor(Aap,x,out=-1)[1]
  
  P  = [speye(size(A,2)); Grad']
  MM(x) = P'*Map(P*x)
  return MM
  
end