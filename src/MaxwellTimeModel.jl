export MaxwellTimeModel, MaxwellTimeModelDerivative, MaxwellTimeModelTransform

# Import stuff from base in order to define multiplication and addition
# for MaxwellTimeModels.
import Base.*
import Base.+
import Base.-
import Base.Ac_mul_B

import jInv.InverseSolve: AbstractModel,AbstractModelDerivative,
                          AbstractModelTransform

# Define earth model type
validProperties = ["sigmaCell","muCell"]
struct MaxwellTimeModel{S<:Real} <: AbstractModel
    values::Dict{String,Vector{S}}
    activeInversionProperties::Vector{String}

    MaxwellTimeModel{S}(sigma::Vector{S},mu::Vector{S}) where {S<:Real} =
      new(Dict{String,Vector{S}}("sigmaCell"=>sigma,"muCell"=>mu),["sigmaCell"])

    MaxwellTimeModel{S}(values::Dict{String,Vector{S}},
      activeProps::Vector{String}) where{S<:Real} =
            new(values,activeProps)
end

MaxwellTimeModel{S}(sigma::Vector{S},mu::Vector{S}) =
  MaxwellTimeModel{S}(sigma,mu)

MaxwellTimeModel{S}(values::Dict{String,Vector{S}},activeProps::Vector{String}) =
  MaxwellTimeModel{S}(values,activeProps)
function MaxwellTimeModel()
    return MaxwellTimeModel(Dict{String,Vector{Float64}}(),Vector{String}())
end

struct MaxwellTimeModelDerivative <: AbstractModelDerivative
    mats::Dict{String,AbstractArray}
    activeInversionProperties::Vector{String}
end

struct MaxwellTimeModelTransform{T<:Real,N<:Integer} <: AbstractModelTransform
    PCell::SparseMatrixCSC{T,N}
    PEdge::SparseMatrixCSC{T,N}
end

import Base.isempty
isempty(m::MaxwellTimeModel) = isempty(m.values)

# Define addition of two models
function +(m1::MaxwellTimeModel,m2::MaxwellTimeModel)
    values = merge(+,m1.values,m2.values)
    active = unique([m1.activeInversionProperties;m2.activeInversionProperties])
    return MaxwellTimeModel(values,active)
end

function -(m1::MaxwellTimeModel,m2::MaxwellTimeModel)
    values = merge(-,m1.values,m2.values)
    active = unique([m1.activeInversionProperties;m2.activeInversionProperties])
    return MaxwellTimeModel(values,active)
end

# Define multiplication with a MaxwellTimeModel
function *(a::Union{Real,UniformScaling}, m::MaxwellTimeModel)
    for key in keys(m.values)
        m.values[key] = a*m.values[key]
    end
  return m
end

function *{T<:Real}(D::MaxwellTimeModelDerivative, x::Vector{T})
    values = Dict{String,Vector{T}}()
    for key in D.activeInversionProperties
        values[key] = D.mats[key]*x
    end
    MaxwellTimeModel(values,D.activeInversionProperties)
end

function Ac_mul_B(D::MaxwellTimeModelDerivative, y::MaxwellTimeModel)
    key1 = first(keys(D.mats))
    T    = eltype(D.mats[key1])
    xOut = zeros(T,size(D.mats[key1],2))
    for key in D.activeInversionProperties
        xOut += D.mats[key]'*y.values[key]
    end
    return xOut
end

# Define multiplication by interpolation matrix
function *{T}(P::MaxwellTimeModelTransform, x::MaxwellTimeModel{T})
    values = Dict{String,Vector{T}}()
    for key in x.activeInversionProperties
        if in(key,["sigmaCell","muCell"])
            values[key] = P.PCell*x.values[key]
        elseif in(key,["sigmaEdge","muEdge"])
            values[key] = P.PEdge*x.values[key]
        end
    end
    return MaxwellTimeModel(values, x.activeInversionProperties)
end

function Ac_mul_B{T,N}(P::MaxwellTimeModelTransform{T,N}, y::MaxwellTimeModel{T})
    values = Dict{String,Vector{T}}()
    for key in y.activeInversionProperties
        if in(key,["sigmaCell","muCell"])
            values[key] = P.PCell'*y.values[key]
        elseif in(key,["sigmaEdge","muEdge"])
            values[key] = P.PEdge'*y.values[key]
        end
    end
    return MaxwellTimeModel(values, y.activeInversionProperties)
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
