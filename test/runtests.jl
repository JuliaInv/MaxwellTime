using MaxwellTime
using Base.Test

# write your own tests here
println(" ")
println("==========  Test BE sensitivities Ind Src ======================")
println(" ")
include("testMaxwellTimeIndSrc.jl")
println(" ")
println("==========  Test BDF2 sensitivities Ind Src ======================")
println(" ")
include("testMaxwellTimeIndSrc-BDF2.jl")
println(" ")
println("==========  Test BE sensitivities Grd Src ======================")
println(" ")
include("testMaxwellTimeGrdSrc.jl")
println(" ")
println("==========  Test BDF2 sensitivities Grd Src ======================")
println(" ")
include("testMaxwellTimeGrdSrc-BDF2.jl")
