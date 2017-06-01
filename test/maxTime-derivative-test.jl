using MaxwellTime

export checkDerivativeMax

import jInv.Utils.getRandomTestDirection

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
		ft = f(x0+10.0^(-j)*v)                # function value
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

function checkDerivativeMax(f::Function,x0in::MaxwellTimeModel,x0bg;out::Bool=true,tol::Float64=1.9,nSuccess::Int=3)
# checkDerivative(f::Function,x0;out::Bool=true,tol::Float64=1.9,nSuccess::Int=3)
	mu0 = pi*4e-7
	if x0in.invertSigma & ~x0in.invertMu
	    x0 = x0in.sigma
	    v  = getRandomTestDirection(x0in.sigma)
	    modfun = x-> begin
                             m = MaxwellTimeModel(exp(x),zeros(length(x)),true,false)
                             Dsig = Vector();
                             push!(Dsig,spdiagm(exp(x)))
                             push!(Dsig,UniformScaling(1.0))
                             return m,Dsig
                         end
	elseif ~x0in.invertSigma & x0in.invertMu
	    x0 = x0in.mu
	    v  = 500*getRandomTestDirection(x0in.mu)
	    modfun = x-> begin
                             m = MaxwellTimeModel(zeros(length(x)),mu0*(1+x),false,true)
                             Dsig = Vector();
                             push!(Dsig,UniformScaling(1.0))
                             dmudm = spdiagm(fill(pi*4e-7,length(x)))
                             push!(Dsig,dmudm)
                             return m,Dsig
                         end
	elseif x0in.invertSigma & x0in.invertMu
	    x0 = [x0in.sigma;x0in.mu]
	    v1 = getRandomTestDirection(x0in.sigma)
	    v2 = getRandomTestDirection(x0in.mu)
	    v  = [v1;v2]
	    modfun = x-> begin
	                     n = length(x0in.sigma)
                             m = MaxwellTimeModel(exp(x[1:n]),mu0*(1+x[n+1:end]),true,true)
                             Dsig = Vector();
                             push!(Dsig,[spdiagm(exp(x[1:n])) spzeros(n,n)])
                             push!(Dsig,[spzeros(n,n) spdiagm(fill(pi*4e-7,n))])
                             return m,Dsig
                         end
	end
	if out
		println(@sprintf("%9s\t%9s\t%9s\t%9s\t%9s\t%5s","h","E0","E1","O1","O2","OK?"))
	end
	sig,dsig = modfun(x0)
	sigloc   = interpGlobalToLocal(sig,speye(length(sig.sigma)),x0bg)
	vloc     = interpGlobalToLocal(dsig,sig,v,speye(length(sig.sigma)))
	f0,dvf  = f(sigloc,vloc)
	f0      = vec(f0)
	nf0     = norm(f0)
	#dvf    = real(dvf)
	Error   = zeros(10,2)
	Order   = zeros(10,2)
	Success = zeros(10)
	for j=1:10
	        sigd,   = modfun(x0+10.0^(-j)*v)
	        sigdloc = interpGlobalToLocal(sigd,speye(length(sigd.sigma)),x0bg)
		ft = f(sigdloc)                # function value
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
 
