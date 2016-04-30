export magnetostaticsCurrentPoly

function magnetostaticsCurrentPoly(x, p)
# A = magnetostaticsCurrentPoly(x, p)
# Compute magnetic vector potential of a unit line current flowing along
# straight segments given in p.

m = size(x, 1)
n = size(p, 1)

A = zeros(m, 3)
# Sum over all line segments
for i = 1:n-1
  a   = p[i,   :]
  b   = p[i+1, :]
  xa  = x - repmat(a, m, 1)
  nxa = vec(pythag(xa))

  xb  = x - repmat(b, m, 1)
  nxb = vec(pythag(xb))

  t   = b - a
  nt  = pythag(t); nt = nt[1]
  if nt == 0
    warning("Polygon segment of zero length detected."); #ok<WNTAG>
  end

  t   = t / nt
  txa = vec(xa * t')
  txb = vec(xb * t')

  jn = (nxa .> 0) & (nxb .> 0)
  jb = txb .> 0

  j  = vec(jn &  jb & (nxa .!= -txa))
  A[j, :] = A[j, :] + log((nxa[j] + txa[j]) ./ (nxb[j] + txb[j])) * t
  j  = vec(jn & !jb & (nxa .!=  txa))
  A[j, :] = A[j, :] - log((nxa[j] - txa[j]) ./ (nxb[j] - txb[j])) * t

end

A = A * 1e-7; # factor mu / 4 / pi = 4e-7 * pi / 4 / pi

return A

end

function pythag(x)
n = size(x, 2)
z = maximum(abs(x), 2)
y = z .* sqrt(sum((x./repmat(z, 1, n)).^2, 2))
y[z .== 0] = 0

return y

end
