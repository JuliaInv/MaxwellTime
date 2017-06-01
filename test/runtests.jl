using MaxwellTime
using Base.Test

include("maxTime-derivative-test.jl")

@testset "MaxwellTime tests" begin
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
end