function ReducedMeshVectorPotential{T<:Real}(mesh::AbstractMesh,
                                             a::T,x0::Vector{T})
    Ex,Ey,Ez = getEdgeGrids(mesh)
    Ne,Qe,Ce,p,Qx,Qy,Qz = getEdgeConstraints(eltype(mesh.h),mesh.S)
    Ex = Qx*Ex
    Ey = Qy*Ey
    Ez = Qz*Ez
    Ax     = zeros(size(Ex,1))
    Ay     = zeros(size(Ey,1))
    mu0    = pi*4e-7

    for (Ej,j) in zip((Ex,Ey),1:2)

        x = Ej[:,1]-x0[1]; y = Ej[:,2]-x0[2]; z = Ej[:,3]-x0[3]

        rho = sqrt.(x.^2 + y.^2)
        m   = 4*a*rho./((a+rho).^2 + z.^2)

        m[m .> 1.0] = 1.0 # m might be slightly larger than 1 due to rounding errors
                      # but ellipke requires 0 <= m <= 1

        # Compute the elliptic integrals
        K,E = ellipke(m)

        # 1/r singular at r = 0 and K(m) singular at m = 1
        i = find((rho .> 0.0) .& (m .< 1.0))

        # Compute Aphi (cylindrical coordinates)
        Aphi = zeros(mesh.ne[j])
        # Aphi[i] = (a*mu0/pi)./sqrt.((a + rho).^2 + z.^2).*
        #            (((2-m).*K[i] - 2*E[i])./m)
        Aphi[i] = (4e-7)./sqrt.(m[i]).*sqrt.(a./rho[i]).*
                  ((1-m[i]/2).*K[i] - E[i])

        # Convert to cartesian coordinates
        if j == 1
            Ax[i]   = Aphi[i] .* (-y[i] ./ rho[i])
        else
            Ay[i]   = Aphi[i] .* ( x[i] ./ rho[i])
        end
    end
    return [Ax;Ay;zeros(size(Ez,1))]
end
