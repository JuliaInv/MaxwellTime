export getNorandaModel
# generate a survey to invert and we can scale
#
using MAT
using Mesh

function getNorandaModel()

# load the noranda example
fid = matopen("noranda512X512X128i.mat");
mm  = read(fid, "model");
close(fid)

# prepare an OcTree that captures this model
mm = reshape(mm,512,512,128)

# find the surface
surface  = zeros(UInt8,512,512)

for i=1:512;
	for j=1:512;
		surface[i,j] = 1; 
		for k=1:128; 
			if mm[i,j,k] == 0
				 surface[i,j] = k; 
				 break; 
			 end; 
		 end
	 end
 end

S  = regularizeOcTree(createOcTreeFromImage(uint8(mm),10));
h  = [50  50  25.0]
n  = size(S)

# project the model onto the mesh and change to conductivity
ii    = find(S)
model = mm[ii]

isactive = (model .> 0) 
inactive = (model .== 0) 

sigma = 10.^(mm/255*(-0.6990+8)-8)
sigAir = 1e-8

sigmaBack = zeros(size(sigma))
sigmaBack[find(inactive.==true)] = sigAir

return sigma, sigmaBack, isactive, h, n, surface
end







 