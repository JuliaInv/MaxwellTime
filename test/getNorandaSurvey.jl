export getNorandaSurvey

function getNorandaSurvey(nsx,nsy,Msh,surface)

   n1 = 512; n2 = 512; n3 = 128
	h1 = 50.0;  h2 = 50.0;  h3 = 25.0
	L1 = n1*h1; L2 = n2*h2; L3 = n3*h3
	
	loopSize = 50
	
	xmin = 1000
	xmax = L1-1000
	ymin = 1000
	ymax = L2-1000
	msh2D    = getRegularMesh([xmin,xmax,ymin,ymax],[nsx,nsy])
	msh2Dall = getRegularMesh([0,L1,0,L2],[n1,n2])
	XY       = getNodalGrid(msh2D)
   XYall    = getNodalGrid(msh2Dall)

	# source locations (center of dipole)
   SrcLoc = zeros(Float64,size(XY,1),3)
	for i=1:size(SrcLoc,1)
		p = XY[i,1]; q = XY[i,2]
		tmp = (p-XYall[:,1]).^2 + (q-XYall[:,2]).^2
		isurf  = indmin(tmp)
		flightHeight = surface[isurf]
		SrcLoc[i,:] = [XY[i,1] XY[i,2] flighHeight]
	end

	# receiver locations (center of dipole)
	RecLoc = SrcLoc
	

	SrcRecMap = speye(size(SrcLoc,1))

	# setup source dipoles
	Srcs   = Array(Array{Float64},size(SrcLoc,1))
	for i=1:size(SrcLoc,1)
		Srcs[i]  = [
							SrcLoc[i,1]-loopSize/2 SrcLoc[i,2]-loopSize/2  SrcLoc[i,3];
							SrcLoc[i,1]+loopSize/2 SrcLoc[i,2]-loopSize/2  SrcLoc[i,3];
							SrcLoc[i,1]+loopSize/2 SrcLoc[i,2]+loopSize/2  SrcLoc[i,3];
							SrcLoc[i,1]-loopSize/2 SrcLoc[i,2]+loopSize/2  SrcLoc[i,3];
							SrcLoc[i,1]-loopSize/2 SrcLoc[i,2]-loopSize/2  SrcLoc[i,3];
						]
	end

	# setup dipole receivers
	Recs = Srcs

	return SrcLoc,Srcs,RecLoc,Recs,SrcRecMap,freq
	end

end