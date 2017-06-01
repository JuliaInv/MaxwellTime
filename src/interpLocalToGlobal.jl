export interpLocalToGlobal, interpGlobalToLocal

# Below methods rely on MaxwellTimeModel * methods defined in MaxwellTimeModel.jl

function interpGlobalToLocal{T}(x::MaxwellTimeModel{T}, P, 
                                y0::MaxwellTimeModel{T})::MaxwellTimeModel{T}
    ym    = P'*x
    sigma = x.invertSigma ? y0.sigma + ym.sigma : y0.sigma
    mu    = x.invertMu    ? y0.mu    + ym.mu    : y0.mu
    return MaxwellTimeModel(sigma,mu,x.invertSigma,x.invertMu)
end

function interpGlobalToLocal{T}(Dsig,sigdum::MaxwellTimeModel{T},
                                x::Vector{T}, P)::MaxwellTimeModel{T}
    if size(Dsig)!=(2,)
        error("interpGlobalToLocal: Must input Dsig/dm and Dmu/dm in vector of length 2")
    end
    dsigdm = sigdum.invertSigma ? Dsig[1]*x : Vector{T}()
    dmudm  = sigdum.invertMu    ? Dsig[2]*x : Vector{T}()
    dm     = MaxwellTimeModel(dsigdm,dmudm,sigdum.invertSigma,sigdum.invertMu)
    return P'*dm
end

# function interpGlobalToLocal{T,N}(x::MaxwellTimeModel{T}, P)::MaxwellTimeModel{T}
#     return xloc = P'*x
# end




function interpLocalToGlobal{T}(Dsig,x::MaxwellTimeModel{T}, P)::Vector{T}
    if size(Dsig)!=(2,)
        error("interpGlobalToLocal: Must input Dsig/dm and Dmu/dm as Vector{MatrixOrScaling}")
    end
    px = P*x
    JTxSig = x.invertSigma ? Dsig[1]'*px.sigma : 0
    JTxMu  = x.invertMu    ? Dsig[2]'*px.mu    : 0
    return JTxSig + JTxMu
end