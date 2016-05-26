module MaxwellTime
	
	using jInv.Mesh.AbstractMesh
	using JOcTree
	using jInv.Utils
	using jInv.LinearSolvers
	using KrylovMethods
			
	export MaxwellTimeParam, getMaxwellTimeParam,MaxwellTimeBDF2Param,
	       getMaxwellTimeBDF2Param,MaxwellTimeTRBDF2Param,getMaxwellTimeTRBDF2Param
	import jInv.ForwardShare.ForwardProbType
	export ForwardProbType

	#Maxwell time domain with backward Euler
	type MaxwellTimeParam <: ForwardProbType
		M::AbstractMesh
		Sources::AbstractArray{Float64}
		Obs::AbstractArray{Float64}
		dt::Vector{Float64}
		wave::Vector{Float64}
		EMsolver::Union{AbstractSolver,Array}
		DCsolver::AbstractSolver
		storeDCfactors::Bool
		fields::Array{Float64}
	end
	function getMaxwellTimeParam(M::AbstractMesh,
				  Sources,
				  Obs,
				  dt,
				  wave,
				  EMsolver::Union{AbstractSolver,Array}=getMUMPSsolver(),
				  DCsolver::AbstractSolver=getMUMPSsolver(),
				  storeDCfactors::Bool=true)

		return MaxwellTimeParam(M,Sources,Obs,dt,wave,EMsolver,
		                        DCsolver,storeDCfactors,Float64[])
	end

	#Maxwell Time domain with BDF-2 time stepping
	type MaxwellTimeBDF2Param <: ForwardProbType
		M::AbstractMesh
		Sources::AbstractArray{Float64}
		Obs::AbstractArray{Float64}
		dt::Float64
		nt::Int64
		wave::Vector{Float64}
		EMsolver::AbstractSolver
		storeEMfactors::Bool
		DCsolver::AbstractSolver
		storeDCfactors::Bool
		fields::Array{Float64}
		ehat::Array{Float64}
	end
	
	function getMaxwellTimeBDF2Param(M::AbstractMesh,
					Sources::AbstractArray{Float64},
					Obs::AbstractArray{Float64},
					dt::Float64,
					nt::Int64,
					wave::Vector{Float64},
					EMsolver::AbstractSolver=getMUMPSsolver(),
					storeEMfactors::Bool=false)

		return MaxwellTimeBDF2Param(M,Sources,Obs,dt,nt,wave,
		         EMsolver,storeEMfactors,getMUMPSsolver(),
		         false,Float64[],Float64[])
	end
	
	function getMaxwellTimeBDF2Param(M::AbstractMesh,
				Sources::AbstractArray{Float64},
				Obs::AbstractArray{Float64},
				dt::Float64,
				nt::Int64,
				wave::Vector{Float64},
				EMsolver::AbstractSolver=getMUMPSsolver(),
				storeEMfactors::Bool=true,
				DCsolver::AbstractSolver=getMUMPSsolver(),
				storeDCfactors::Bool=true)

	return MaxwellTimeBDF2Param(M,Sources,Obs,dt,nt,wave,
	         EMsolver,storeEMfactors,DCsolver,storeDCfactors,
	         Float64[],Float64[])
	end
	
	#Maxwell Time domain with TR-BDF2 time stepping
	type MaxwellTimeTRBDF2Param <: ForwardProbType
		M::AbstractMesh
		Sources::AbstractArray{Float64}
		Obs::AbstractArray{Float64}
		dt::Float64
		nt::Int64
		wave::Vector{Float64}
		fields::Array{Float64}
		fname::AbstractString
		solver::AbstractSolver
	end
	
	function getMaxwellTimeTRBDF2Param(M::AbstractMesh,
										  Sources,
										  Obs,
										  dt,
										  nt,
										  wave,
										  solver;
										  fields=Array(Float64,0,0),
							          	  fname="")

		return MaxwellTimeTRBDF2Param(M,Sources,Obs,dt,nt,wave,fields,fname,solver)
	end
	
	# Maxwell Time domain with explicit sensitivity storage
	export MaxwellTimeSEParam, getMaxwellTimeSEParam
	type MaxwellTimeSEParam <: ForwardProbType
		M::AbstractMesh
		Sources::AbstractArray{Float64}
		Obs::AbstractArray{Float64}
		dt::Vector{Float64}
		wave::Vector{Float64}
		Sens::Array{Float64}
		fields::Array{Float64}
		fname::AbstractString
		solver::AbstractSolver
	end # type MaxwellTimeSEParam
	
	function getMaxwellTimeSEParam(M::AbstractMesh,
										  Sources,
										  Obs,
										  dt,
										  wave,
										  solver;
										  Sens=Array(Float64,0,0),
										  fields=Array(Float64,0,0),
							          	  fname="")

		return MaxwellTimeSEParam(M,Sources,Obs,dt,wave,Sens,fields,fname,solver)
	end
	
	import jInv.ForwardShare.getData
	import jInv.ForwardShare.getSensMatVec
	import jInv.ForwardShare.getSensTMatVec
	
	include("getData.jl")
	include("getSensMatVec.jl")
	include("getSensTMatVec.jl")
	include("magnetostaticsCurrentPoly.jl")
        include("getDataAlt.jl")
        include("solveMaxTime.jl")
        include("solveDC.jl")
	
	
end
