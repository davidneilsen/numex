
function init_data!(fields)
    nx = fields.grid.nx
    ny = fields.grid.ny
    nz = fields.grid.nz

    Ex = fields.u[1]
    Ey = fields.u[2]
    Bz = fields.u[3]
    x = fields.grid.x
    y = fields.grid.y
    z = fields.grid.z

    amp1 = 1.0
    lambda1 = 1.0

    for k=1:nz, j=1:ny, i=1:nx
        r = sqrt(x[i]*x[i] + y[j]*y[j] + z[k]*z[k])
        Bz[i,j,k] = - 8.0*amp1*lambda1*lambda1*exp(-lambda1*r*r)
    end

end

function Xinit_data!(fields)
    # Ex and Ey are set at t = 0.0
    # Bz is set at  t = -dt/2

    md = 2.0
    nd = 2.0
    c = 1.0

    dimx = fields.grid.xmax - fields.grid.xmin
    dimy = fields.grid.ymax - fields.grid.ymin
    nx = fields.grid.nx
    ny = fields.grid.ny
    dx = fields.grid.dx
    dy = fields.grid.dy
    dt = fields.grid.dt
    omega = c * sqrt((md*pi/dimx)^2+(nd*pi/dimy)^2)

    Ex = fields.u[1]
    Ey = fields.u[2]
    Bz = fields.u[3]

    println("init_data:  nx = ",nx," ny = ",ny)
    for j=1:ny, i=1:nx
        Bz[i,j] = (- cos(md*pi*((i-0.5)*dx/dimx)) 
                        * cos(nd*pi*((j-0.5)*dy/dimy))
                        * cos(omega*(-0.5*dt)) )
#        @printf("i=%d, j=%d, Ex=%g, Ey=%g, Bz = %g",i,j,Ex[i,j],Ey[i,j],Bz[i,j])
    end  
end

function maxwellEqs!(dtu, u, dxu, dyu, x, y, dx, dy, time)
    dtEx = dtu[1]
    dtEy = dtu[2]
    dtEz = dtu[3]
    dtHx = dtu[4]
    dtHy = dtu[5]
    dtHz = dtu[6]

    dxEx = dxu[1]
    dxEy = dxu[2]
    dxEz = dxu[3]
    dxBx = dxu[4]
    dxBy = dxu[5]
    dxBz = dxu[6]

    dyEx = dyu[1]
    dyEy = dyu[2]
    dyEz = dyu[3]
    dyBx = dyu[4]
    dyBy = dyu[5]
    dyBz = dyu[6]

    @. diff22_x!(dxu, u, dx)
    @. diff22_y!(dyu, u, dy)
    @. diff22_z!(dzu, u, dz)

    @. dtEx = dyHz - dzHy
    @. dtEy = dzHx - dxHz
    @. dtEy = dxHy - dyHx

    @. dtHx = dzEy - dyEz
    @. dtHy = dxEz - dzEx
    @. dtHz = dyEx - dxEy

    @. sommerfeld_bcs(dtu, u, dxu, dyu, dzu, x, y, z)

end

function sommerfeld_bcs(dtu, u, dxu, dyu, dzu, x, y, z)

    nx, ny, nz = size(u)

    u0::Float64 = 0.0
    u_falloff::Float64 = 2.0

    ############  j = 1
    j = 1
    for k = 1:ny
        for i = 1:nx
            dtu[i,j,k] = - (x[i]*dxu[i,j,k] + y[j]*dyu[i,j,k] z[k]*dzu[i,j,k] + u_falloff*(u[i,j,k] - u0))/sqrt(x[i]*x[i] + y[j]*y[j] + z[k]*z[k])
        end
    end

    ############  j = ny
    j = ny
    for k = 1:nz
        for i = 1:nx
            dtu[i,j,k] = - (x[i]*dxu[i,j,k] + y[j]*dyu[i,j,k] z[k]*dzu[i,j,k] + u_falloff*(u[i,j,k] - u0))/sqrt(x[i]*x[i] + y[j]*y[j] + z[k]*z[k])
        end
    end

    ############  i = 1
    i = 1
    for k = 1:nz
        for j = 1:ny
            dtu[i,j,k] = - (x[i]*dxu[i,j,k] + y[j]*dyu[i,j,k] z[k]*dzu[i,j,k] + u_falloff*(u[i,j,k] - u0))/sqrt(x[i]*x[i] + y[j]*y[j] + z[k]*z[k])
        end
    end

    ############  i = nx
    i = nx
    for k = 1:nz
        for j = 1:ny
            dtu[i,j,k] = - (x[i]*dxu[i,j,k] + y[j]*dyu[i,j,k] z[k]*dzu[i,j,k] + u_falloff*(u[i,j,k] - u0))/sqrt(x[i]*x[i] + y[j]*y[j] + z[k]*z[k])
        end
    end

    ############  k = 1
    k = 1
    for j = 1:ny
        for i = 1:nx
            dtu[i,j,k] = - (x[i]*dxu[i,j,k] + y[j]*dyu[i,j,k] z[k]*dzu[i,j,k] + u_falloff*(u[i,j,k] - u0))/sqrt(x[i]*x[i] + y[j]*y[j] + z[k]*z[k])
        end
    end

    ############  k = nz
    k = nz
    for j = 1:ny
        for i = 1:nx
            dtu[i,j,k] = - (x[i]*dxu[i,j,k] + y[j]*dyu[i,j,k] z[k]*dzu[i,j,k] + u_falloff*(u[i,j,k] - u0))/sqrt(x[i]*x[i] + y[j]*y[j] + z[k]*z[k])
        end
    end
    
end
