function C(u,sigma::Array{Float64,1},param::MaxwellTimeParam)
	dt       = param.dt
	nt       = length(param.dt)
	wave     = param.wave
	Msh      = param.Mesh

	mu     = 4*pi*1e-7
	Curl   = getCurlMatrix(Msh)
	Msig   = getEdgeMassMatrix(Msh,vec(sigma))
	Mmu    = getFaceMassMatrix(Msh,1/mu*ones(size(sigma)))
        Ne,Qe, = getEdgeConstraints(Msh)
        G      = getNodalGradientMatrix(Msh)
        Nn,Qn, = getNodalConstraints(Msh)
        Nf,Qf, = getFaceConstraints(Msh)

        G     = Qe*G*Nn
        Curl  = Qf*Curl*Ne
        Mmu   = Nf'*Mmu*Nf
        Msig  = Ne'*Msig*Ne
        s     = param.Sources

#         K     = Curl'*Mmu*Curl
#         A     = K + 3/(2*dt)*Msig
#         #Ainv  = inv(A)
#         Ainv   = A\Msig
# #         println(size(Ainv))
# #         println(typeof(Ainv))
# #         println(size(dt))
# #         println(typeof(dt))
#
#         nen  = size(Ne,2)
#         nnn  = size(Nn,2)
#         blnkn = spzeros(nnn,4*nen)
#         blnkee = spzeros(nen,nen)
#         blnken = spzeros(nen,nnn)
#         dCdu  = [G'*Msig*G blnkn;
#                 (3/(4*dt)*Msig+9/(8*dt^2)*Msig*Ainv)*G A blnkee blnkee blnkee;
#                 -1/(2*dt)*Msig*G -2/dt*Msig A blnkee blnkee;
#                 blnken 1/(2*dt)*Msig -2/dt*Msig A blnkee;
#                 blnken blnkee 1/(2*dt)*Msig -2/dt*Msig A]
nen  = size(Ne,2)
nnn  = size(Nn,2)
blnkn = spzeros(nnn,5*nen)
blnkee = spzeros(nen,nen)
blnken = spzeros(nen,nnn)
K     = Curl'*Mmu*Curl
tau   = dt[2:end]./dt[1:end-1]
g12   = (1 + 2*tau[1])/(1+tau[1])
g13   = (1 + 2*tau[2])/(1+tau[2])
g14   = (1 + 2*tau[3])/(1+tau[3])
g22   = 1 + tau[1]
g23   = 1 + tau[2]
g24   = 1 + tau[3]
g32   = (tau[1]^2)/(1+tau[1])
g33   = (tau[2]^2)/(1+tau[2])
g34   = (tau[3]^2)/(1+tau[3])
A1    = K + 3/(2*dt[1])*Msig
A2    = K + 3/(2*dt[2])*Msig
A3    = K + 3/(2*dt[3])*Msig
A4    = K + 3/(2*dt[4])*Msig
dCdu  = [G'*Msig*G blnkn;
         3/(2*dt[1])*Msig*G A1 blnkee blnkee blnkee blnkee;
         3/(4*dt[1])*Msig*G -3/(4*dt[1])*Msig A1 blnkee blnkee blnkee;
         -g32/(dt[2])*Msig*G blnkee -g22/dt[2]*Msig A2 blnkee blnkee;
         blnken blnkee g33/(dt[3])*Msig -g23/dt[3]*Msig A3 blnkee;
         blnken blnkee blnkee g34/(dt[4])*Msig -g24/dt[4]*Msig A4]


        q = [G'*vec(s);zeros(5*size(Curl,2))]
        return dCdu*u-q
end
