using MaxwellTime

export checkDerivativeMax

import jInv.Utils.getRandomTestDirection
import jInv.ForwardShare: interpGlobalToLocal,interpLocalToGlobal

function checkDerivativeMax(f::Function,df::Function,x0;kwargs...)
	function testFun(x,v=[])
		if v != []
			return f(x), df(v,x)
		else
			return f(x)
		end
	end
	return checkDerivativeMax(testFun,x0;kwargs...)
end

function checkDerivativeMax(f::Function,df::Function,x0,xbg;kwargs...)
	function testFun(x,v=[])
		if v!=[]
			return f(x), df(v,x)
		else
			return f(x)
		end
	end
	return checkDerivativeMax(testFun,x0,xbg;kwargs...)
end



function checkDerivativeMax(f::Function,x0;out::Bool=true,tol::Float64=1.9,nSuccess::Int=3,v=getRandomTestDirection(x0))
# checkDerivative(f::Function,x0;out::Bool=true,tol::Float64=1.9,nSuccess::Int=3)
	if out
		println(@sprintf("%9s\t%9s\t%9s\t%9s\t%9s\t%5s","h","E0","E1","O1","O2","OK?"))
	end

	f0,dvf  = f(x0,v)
	f0      = vec(f0)
	nf0     = norm(f0)
	#dvf    = real(dvf)
	Error   = zeros(10,2)
	Order   = zeros(10,2)
	Success = zeros(10)
	for j=1:10
		ft = f(x0+2.0^(-j)*v)                # function value
		ft = vec(ft)
		Error[j,1] = norm(f0-ft)/nf0          # Error TaylorPoly 0
		Error[j,2] = norm(f0 .+2.0^(-j)*dvf .- ft)/nf0 # Error TaylorPoly 1
		if j>1
			Order[j,:] = log2.(Error[j-1,:]./Error[j,:]);
		end
		if (Order[j,2]>tol) || (Error[j,1]/Error[j,2] > 4); Success[j]=1; end
		if out
			println(@sprintf("%1.3e\t%1.3e\t%1.3e\t%1.3e\t%1.3e\t%5d",
							2.0^(-j), Error[j,1],Error[j,2], Order[j,1],Order[j,2],Success[j]))
		end
	end
	pass = sum(Success) > nSuccess
	return  pass,Error,Order
end

function checkDerivativeMaxMu(f::Function,x0::MaxwellTimeModel,pFor;out::Bool=true,tol::Float64=1.9,nSuccess::Int=3)
    chi  = 5*(10.^getRandomTestDirection(x0.mu))
    #chi  = 1.;
    v    = x0.mu .* chi
    vmod = MaxwellTimeModel(zeros(length(v)),v,false,true)
#     if any((x0.mu-0.1*v) .< pi*4e-7)
#         error("mu < mu0, aborting test")
#     end
    f0     = vec(f(x0.mu))
    dvf    = getSensMatVec(vmod,x0,pFor)
    nf0    = norm(f0)
    Error   = zeros(10,2)
    Order   = zeros(10,2)
    Success = zeros(10)
    for j=1:10
    	ft = f(x0.mu+10.0^(-j)*v)                # function value
    	ft = vec(ft)
    	Error[j,1] = norm(f0-ft)/nf0          # Error TaylorPoly 0
    	Error[j,2] = norm(f0 .+10.0^(-j)*dvf .- ft)/nf0 # Error TaylorPoly 1
    	if j>1
    		Order[j,:] = log10(Error[j-1,:]./Error[j,:]);
    	end
    	if (Order[j,2]>tol) || (Error[j,1]/Error[j,2] > 100); Success[j]=1; end
    	if out
    		println(@sprintf("%1.3e\t%1.3e\t%1.3e\t%1.3e\t%1.3e\t%5d",
    						10.0^(-j), Error[j,1],Error[j,2], Order[j,1],Order[j,2],Success[j]))
    	end
    end
    pass = sum(Success) > nSuccess
    return  pass,Error,Order
end

function checkDerivativeMax(f::Function,x0in::MaxwellTimeModel,x0bg;
	                        out::Bool=true,tol::Float64=1.9,
							nSuccess::Int=3,base=10.0)
# checkDerivative(f::Function,x0;out::Bool=true,tol::Float64=1.9,nSuccess::Int=3)
	mu0 = pi*4e-7
	invertSigma = in("sigmaCell",x0in.activeInversionProperties)
	invertMu    = in(   "muCell",x0in.activeInversionProperties)
	if invertSigma & ~invertMu
	    x0 = x0in.values["sigmaCell"]
	    v  = getRandomTestDirection(x0in.sigma)
	    modfun = x-> begin
                             m  = MaxwellTimeModel(Dict("sigmaCell"=>exp(x)),["sigmaCell"])
                             dm = MaxwellTimeModelDerivative(Dict("sigmaCell"=>spdiagm(exp.(x))),["sigmaCell"])
                             return m,dm
                         end
	elseif ~invertSigma & invertMu
	    x0 = x0in.values["muCell"]
	    v  = 100*rand(length(x0))
	    modfun = x-> begin
                             m  = MaxwellTimeModel(Dict("muCell"=>mu0*(1+x)),["muCell"])
                             dmudm = spdiagm(fill(pi*4e-7,length(x)))
							 dm = MaxwellTimeModelDerivative(Dict("muCell"=>dmudm),["muCell"])
                             return m,dm
                         end
	elseif invertSigma & invertMu
	    x0 = [x0in.values["sigmaCell"];x0in.values["muCell"]]
	    v1 = getRandomTestDirection(x0in.values["sigmaCell"])
	    v2 = getRandomTestDirection(x0in.values["muCell"])
	    v  = [v1;v2]
	    modfun = x-> begin
	                     n = length(x0in.values["sigmaCell"])
                             m = MaxwellTimeModel(Dict("sigmaCell"=>exp.(x[1:n]),
							       "muCell"=>mu0*(1+x[n+1:end])),["sigmaCell","muCell"])
                             dm = MaxwellTimeModelDerivative(Dict("sigmaCell"=>[spdiagm(exp.(x[1:n])) spzeros(n,n)],
							        "muCell"=>[spzeros(n,n) spdiagm(fill(pi*4e-7,n))]),["sigmaCell","muCell"])
                             return m,dm
                         end
	end
	if out
		println(@sprintf("%9s\t%9s\t%9s\t%9s\t%9s\t%5s","h","E0","E1","O1","O2","OK?"))
	end
	sig,dsig = modfun(x0)
	sigloc   = interpGlobalToLocal(sig,1.0,x0bg)
	vloc     = dsig*v
	f0,dvf  = f(sigloc,vloc)
	f0      = vec(f0)
	nf0     = norm(f0)
	#dvf    = real(dvf)
	Error   = zeros(10,2)
	Order   = zeros(10,2)
	Success = zeros(10)
	for j=1:10
	        sigd,   = modfun(x0+base^(-j)*v)
	        sigdloc = sigd+x0bg #interpGlobalToLocal(sigd,speye(length(sigd.sigma)),x0bg)
		ft = f(sigdloc)                # function value
		ft = vec(ft)
		Error[j,1] = norm(f0-ft)/nf0          # Error TaylorPoly 0
		Error[j,2] = norm(f0 .+base^(-j)*dvf .- ft)/nf0 # Error TaylorPoly 1
		if j>1
			Order[j,:] = log.(base,Error[j-1,:]./Error[j,:]);
		end
		if (Order[j,2]>tol) || (Error[j,1]/Error[j,2] > base^2); Success[j]=1; end
		if out
			println(@sprintf("%1.3e\t%1.3e\t%1.3e\t%1.3e\t%1.3e\t%5d",
							base^(-j), Error[j,1],Error[j,2], Order[j,1],Order[j,2],Success[j]))
		end
	end
	pass = sum(Success) > nSuccess
	return  pass,Error,Order
end
