export EFDmagDipole,EDCDipole,ETDDipole

#Electric field due to unit a vertical magnetic dipole
#in a fullspace
function EFDmagDipole(Ex,Ey,Ez,x0::Vector{Float64},freq::Float64,
                      sigma::Float64)
  @assert freq > 0.0 "EFDmagDipole: Frequency must be > 0"
  w = 2*pi*freq; mu = pi*4e-7
  k = sqrt(-im*mu*sigma*w)
  x = Ex[:,1]-x0[1]; y = Ex[:,2]-x0[2]; z = Ex[:,3]-x0[3]
  r = sqrt(x.^2 + y.^2 + z.^2)
  @assert all(r .> 0.1) "EFDmagDipole: Dipole must be located at face centre, >0.1m from any edge"
  efx = im*w*mu./(4*pi*r.^2).*(im*k*r + 1).*exp(-im*k*r).*y./r

  x = Ey[:,1]-x0[1]; y = Ey[:,2]-x0[2]; z = Ey[:,3]-x0[3]
  r = sqrt(x.^2 + y.^2 + z.^2)
  @assert all(r .> 0.1) "EFDmagDipole: Dipole must be located at face centre, >0.1m from any edge"
  efy = -im*w*mu./(4*pi*r.^2).*(im*k*r + 1).*exp(-im*k*r).*x./r


  x = Ez[:,1]-x0[1]; y = Ez[:,2]-x0[2]; z = Ez[:,3]-x0[3]
  r = sqrt(x.^2 + y.^2 + z.^2)
  @assert all(r .> 0.1) "EFDmagDipole: Dipole must be located at face centre, >0.1m from any edge"
  efz = complex(zeros(length(r)))

  return [efx;efy;efz]
end

#Electric field due to a DC electric dipole in a fullspace
function EDCDipole(Ex,Ey,Ez,x0::Vector{Float64},sigma,
                   ds::Float64,I::Float64)
  mu = 4*pi*1e-7
  x = Ex[:,1]-x0[1]; y = Ex[:,2]-x0[2]; z = Ex[:,3]-x0[3]
  r = sqrt.(x.^2 + y.^2 + z.^2)
  I0 = find(r .< 0.1); r[I0] = 0.1
  dcx = (I*ds./(4*pi*sigma*r.^3)).*(3*(x.^2)./r.^2 -1.0)

  x = Ey[:,1]-x0[1]; y = Ey[:,2]-x0[2]; z = Ey[:,3]-x0[3]
  r = sqrt.(x.^2 + y.^2 + z.^2)
  I0 = find(r .< 0.1); r[I0] = 0.1
  dcy = (I*ds./(4*pi*sigma*r.^3)).*(3*x.*y./r.^2)

  x = Ez[:,1]-x0[1]; y = Ez[:,2]-x0[2]; z = Ez[:,3]-x0[3]
  r = sqrt.(x.^2 + y.^2 + z.^2)
  I0 = find(r .< 0.1); r[I0] = 0.1
  dcz = (I*ds./(4*pi*sigma*r.^3)).*(3*x.*z./r.^2)
  return [dcx;dcy;dcz]
end
#Electric field due to a transient electric dipole in a fullspace
function ETDDipole(Ex,Ey,Ez,x0::Vector{Float64},sigma,
                   ds::Float64,I::Float64,t::Vector{Float64})
  mu = 4*pi*1e-7
  if (t[1] == 0.0)
    nt = length(t)
  else
    nt = length(t) + 1
  end
  efx = zeros(length(Ex[:,1]),nt)
  efy = zeros(length(Ey[:,1]),nt)
  efz = zeros(length(Ez[:,1]),nt)
  for i=2:length(t)
    theta = sqrt.(mu*sigma./(4*t))
    x = Ex[:,1]-x0[1]; y = Ex[:,2]-x0[2]; z = Ex[:,3]-x0[3]
    r = sqrt.(x.^2 + y.^2 + z.^2)
    I0 = find(r .< 0.1); r[I0] = 0.1
    efx[:,i] = (I*ds./(4*pi*sigma*r.^3)).*( ((4/sqrt.(pi)*theta[i]^3*r.^3 +
               6/sqrt.(pi)*theta[i]*r).*exp.(-(theta[i]^2)*r.^2) +
               3*erfc.(theta[i]*r)).*((x.^2)./r.^2) -
               ((4/sqrt.(pi)*theta[i]^3*r.^3 +
               2/sqrt.(pi)*theta[i]*r).*exp.(-(theta[i]^2)*r.^2) +
               erfc.(theta[i]*r)) )

    x = Ey[:,1]-x0[1]; y = Ey[:,2]-x0[2]; z = Ey[:,3]-x0[3]
    r = sqrt.(x.^2 + y.^2 + z.^2)
    I0 = find(r .< 0.1); r[I0] = 0.1
    efy[:,i] = (I*ds./(4*pi*sigma*r.^3)).*( ((4/sqrt.(pi)*theta[i]^3*r.^3 +
               6/sqrt.(pi)*theta[i]*r).*exp.(-(theta[i]^2)*r.^2) +
               3*erfc.(theta[i]*r)).*(x.*y./r.^2) )

    x = Ez[:,1]-x0[1]; y = Ez[:,2]-x0[2]; z = Ez[:,3]-x0[3]
    r = sqrt.(x.^2 + y.^2 + z.^2)
    I0 = find(r .< 0.1); r[I0] = 0.1
    efz[:,i] = (I*ds./(4*pi*sigma*r.^3)).*( ((4/sqrt.(pi)*theta[i]^3*r.^3 +
               6/sqrt.(pi)*theta[i]*r).*exp.(-(theta[i]^2)*r.^2) +
               3*erfc.(theta[i]*r)).*(x.*z./r.^2) )
  end
  eTrans = [efx;efy;efz]
  edc = EDCDipole(Ex,Ey,Ez,x0,sigma,ds,I)
  return edc.-eTrans
end

# Vertical component of dhdt due to a horizontal loop of radius rLoop
# at the surface of a halfspace of conductivity sigma
function VMDdhdtz(X,x0,rLoop,sigma,t)
  mu    = pi*4e-7
  theta = sqrt.(mu*sigma./(4*t))
  rho   = sqrt.((X[:,1]-x0[1]).^2 + (X[:,2] - x0[2]).^2)
  inSide = find(rho .<= rLoop)
  rho[inSide] = rLoop + 1.0
  dhdz  = zeros(length(rho),length(t))
  for i=1:length(t)
    dhdz[:,i] = -1./(2*pi*mu*sigma*rho.^5).*( 9*erf.(theta[i]*rho) -
                     2*theta[i]*rho/sqrt(pi).*(9+6*theta[i]^2*rho.^2 +
                     4*theta[i]^4*rho.^4).*exp.(-theta[i]^2.*rho.^2) )
  end
  return dhdz
end

function VMDdhdtzCentLoop(radius, t, sigma)
    mu = pi*4e-7
    ta = radius*sqrt.((sigma*mu)./(4*t))
    return (3*erf.(ta) - 2./sqrt(pi).*ta.*(3 + 2*ta.^2).*exp.(-ta.^2))./(mu*sigma*radius^3)
end

function hzAnalyticCentLoop(radius, t, sigma)
    mu    = pi*4e-7
    theta = sqrt.((sigma*mu)./(4*t))
    ta = radius*theta
    eta = erf.(ta)
    t1 = (3./(sqrt(pi)*ta)).*exp.(-ta.^2)
    t2 = (1 - (3./(2*ta.^2))).*eta
    hz = (t1 + t2)/(2*radius)
    return hz
end
