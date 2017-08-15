module MaxwellTime

using jInv.Mesh.AbstractMesh
using JOcTree
using jInv.Utils
using jInv.LinearSolvers
using KrylovMethods

hasMUMPS = false
try
    using MUMPS
    hasMUMPS = true
catch
end

hasPardiso = false
try
    using Pardiso
    hasPardiso = true
catch
end

# We will add MaxwellTime specific subtypes/methods
# to the following jInv generics
import jInv.ForwardShare.ForwardProbType
import jInv.ForwardShare.getData
import jInv.ForwardShare.getSensMatVec
import jInv.ForwardShare.getSensTMatVec
import jInv.ForwardShare.interpGlobalToLocal
import jInv.ForwardShare.interpLocalToGlobal
import jInv.InverseSolve.AbstractModel
#export ForwardProbType

# Import MatrixOrScaling typealias
import jInv.ForwardShare.MatrixOrScaling

# Define MaxwellTime forward problem param type and associated setup functions
include("MaxwellTimeParam.jl")

# Define MaxwellTimeModel earth model type and associated methods
include("MaxwellTimeModel.jl")

# Define time-stepping (integration) functions and map
# integration method symbols to the appropriate
# functions defined in getFields.jl
# getFields.jl also contains matrix construction helper
# functions used in both forward modelling and sensitivity
# computations
include("getFields.jl")
integrationFunctions = Dict(zip(supportedIntegrationMethods,
                                [getFieldsBE;getFieldsBDF2;
                                 getFieldsBDF2ConstDT;
                                 getFieldsTRBDF2]))

# Forward solve, sensitivity, and mesh2mesh interpolation routines
include("getData.jl")
include("getSensMatVec.jl")
include("getSensTMatVec.jl")
include("interpLocalToGlobal.jl") # Specific methods for MaxwellTimeModel
                                  # mesh 2 mesh interpolation

# Interface MaxwellTime to jInv.LinearSolvers to solve linear systems
# of equations
include("solverFunctions.jl")

end #End module MaxwellTime
