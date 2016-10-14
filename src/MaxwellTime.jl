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
	type MaxwellTimeParam{S<:Real} <: ForwardProbType
		M::AbstractMesh
		Sources::AbstractArray{S}
		Obs::AbstractArray{S}
		dt::Vector{S}
		wave::Vector{S}
		EMsolver::Union{AbstractSolver,Array}
		DCsolver::AbstractSolver
		storeDCfactors::Bool
		fields::Array{S}
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
	type MaxwellTimeBDF2Param{S<:Real} <: ForwardProbType
		M::AbstractMesh
		Sources::AbstractArray{S}
		Obs::AbstractArray{S}
		dt::S
		nt::Int64
		wave::Vector{S}
		EMsolver::AbstractSolver
		storeEMfactors::Bool
		DCsolver::AbstractSolver
		storeDCfactors::Bool
		fields::Array{S}
		ehat::Array{S}
	end
	
	
	function getMaxwellTimeBDF2Param{S<:Real}(M::AbstractMesh,
				Sources::AbstractArray{S},
				Obs::AbstractArray{S},
				dt::S,
				nt::Int64,
				wave::Vector{S},
				EMsolver::AbstractSolver=getMUMPSsolver(),
				storeEMfactors::Bool=true,
				DCsolver::AbstractSolver=getMUMPSsolver(),
				storeDCfactors::Bool=true)

	return MaxwellTimeBDF2Param(M,Sources,Obs,dt,nt,wave,
	         EMsolver,storeEMfactors,DCsolver,storeDCfactors,
	         S[],S[])
	end
	
	#Maxwell Time domain with TR-BDF2 time stepping
	type MaxwellTimeTRBDF2Param{S<:Real} <: ForwardProbType
		M::AbstractMesh
		Sources::AbstractArray{S}
		Obs::AbstractArray{S}
		dt::S
		nt::Int64
		wave::Vector{S}
		fields::Array{S}
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
	type MaxwellTimeSEParam{S<:Real} <: ForwardProbType
		M::AbstractMesh
		Sources::AbstractArray{S}
		Obs::AbstractArray{S}
		dt::Vector{S}
		wave::Vector{S}
		Sens::Array{S}
		fields::Array{S}
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
        include("solveMaxTime.jl")
        include("solveDC.jl")
	
	
end
