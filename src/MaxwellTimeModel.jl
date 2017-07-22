export MaxwellTimeModel

# Import stuff from base in order to define multiplication and addition
# for MaxwellTimeModels.
import Base.*
import Base.+

# Define earth model type
struct MaxwellTimeModel{S<:Real} <: AbstractModel
    sigma::Vector{S}
    mu::Vector{S}
    invertSigma::Bool
    invertMu::Bool

    MaxwellTimeModel{S}(sigma::Vector{S},mu::Vector{S}) where {S<:Real} =
      new(sigma,mu,true,false)

    MaxwellTimeModel{S}(sigma::Vector{S},mu::Vector{S},invertSigma::Bool,invertMu::Bool) where{S<:Real} = 
            new(sigma,mu,invertSigma,invertMu)
end

MaxwellTimeModel{S}(sigma::Vector{S},mu::Vector{S}) =
  MaxwellTimeModel{S}(sigma,mu)

MaxwellTimeModel{S}(sigma::Vector{S},mu::Vector{S},invertSigma::Bool,invertMu::Bool) =
  MaxwellTimeModel{S}(sigma,mu,invertSigma,invertMu)

# Define addition of two models
function +(m1::MaxwellTimeModel,m2::MaxwellTimeModel)
    invSig = m1.invertSigma || m2.invertSigma
    invMu  = m1.invertMu    || m2.invertMu
    return MaxwellTimeModel(m1.sigma+m2.sigma,m1.mu+m2.mu,invSig,invMu)
end

# Define multiplication with a MaxwellTimeModel
function *{T1<:Real,T2}(a::T1,m::MaxwellTimeModel{T2})
    aT2 = convert(T2,a)
    return MaxwellTimeModel(aT2*m.sigma,aT2*m.mu,m.invertSigma,m.invertMu)
end

function *{T,N}(A::SparseMatrixCSC{T,N},m::MaxwellTimeModel{T})
    sigma = m.invertSigma ? A*m.sigma : Vector{T}()
    mu    = m.invertMu    ? A*m.mu    : Vector{T}()
    return MaxwellTimeModel(sigma,mu,m.invertSigma,m.invertMu)
end


function *{T}(D::Vector,m::MaxwellTimeModel{T})::MaxwellTimeModel{T}
    #T = promote_type(T1,T2)
    if length(D) == 1
        return MaxwellTimeModel(D[1]*m.sigma,D[1]*m.mu,m.invertSigma,m.invertMu)
    elseif length(D) == 2
        return MaxwellTimeModel(D[1]*m.sigma,D[2]*m.mu,m.invertSigma,m.invertMu)
    else
        error("*: MaxwellTimeModel cannot be multiplied with vector of matrices with length > 2")
    end
end

# Below code is for experimenting with type stability

# function *(D::Vector,m::MaxwellTimeModel)
#     #T = promote_type(T1,T2)
#     if length(D) == 1
#         #return MaxwellTimeModel(D[1]*m.sigma,D[1]*m.mu,m.invertSigma,m.invertMu)
#         return D_mul_m(D[1],m)
#     elseif length(D) == 2
#         #return MaxwellTimeModel(D[1]*m.sigma,D[2]*m.mu,m.invertSigma,m.invertMu)
#         return D_mul_m(D[1],D[2],m)
#     else
#         error("*: MaxwellTimeModel cannot be multiplied with vector of matrices with length > 2")
#     end
# end
#
# function D_mul_m{T,N}(D1::SparseMatrixCSC{T,N},m::MaxwellTimeModel{T})
#     return MaxwellTimeModel(D1*m.sigma,D2*m.sigma,m.invertSigma,m.invertMu)
# end
#
# function D_mul_m{T}(D1::UniformScaling{T},m::MaxwellTimeModel{T})
#     return MaxwellTimeModel(D1*m.sigma,D1*m.sigma,m.invertSigma,m.invertMu)
# end
#
#
# function D_mul_m{T,N}(D1::SparseMatrixCSC{T,N},D2::SparseMatrixCSC{T,N},m::MaxwellTimeModel{T})
#     return MaxwellTimeModel(D1*m.sigma,D2*m.sigma,m.invertSigma,m.invertMu)
# end
#
# function D_mul_m{T,N}(D1::SparseMatrixCSC{T,N},D2::UniformScaling{T},m::MaxwellTimeModel{T})
#     return MaxwellTimeModel(D1*m.sigma,D2*m.sigma,m.invertSigma,m.invertMu)
# end
#
# function D_mul_m{T,N}(D1::UniformScaling{T},D2::SparseMatrixCSC{T,N},m::MaxwellTimeModel{T})
#     return MaxwellTimeModel(D1*m.sigma,D2*m.sigma,m.invertSigma,m.invertMu)
# end
