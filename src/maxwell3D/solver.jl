using Printf
using WriteVTK

include("maxwell.jl")
import .Maxwell

function evolve!(fields, nt, VTKOutFreq)
    time = [0.0]
    vtkFileCount::Int64 = 0
    screenOutFreq = VTKOutFreq
    filename = @sprintf("maxwell_%05d",vtkFileCount)
    vtk_grid(filename, fields.grid.x, fields.grid.y) do vtk
        vtk["Ex",VTKPointData()] = fields.u[1]
        vtk["Ey",VTKPointData()] = fields.u[2]
        vtk["Hz",VTKPointData()] = fields.u[3]
        vtk["time",VTKFieldData()] = time[1]
    end
    for i = 1:nt
        Maxwell.rk2_step!(Maxwell.maxwell_TE!, fields, time)
        if (mod(i,screenOutFreq)==0)
            @printf("Step=%d, time=%g, |Bz|=%g\n",i,time[1],Maxwell.l2norm(fields.u[3]))
        end
        if (mod(i,vtkOutFreq)==0)
            vtkFileCount += 1
            filename = @sprintf("maxwell_%05d",vtkFileCount)
            vtk_grid(filename, fields.grid.x, fields.grid.y) do vtk
                vtk["Ex",VTKPointData()] = fields.u[1]
                vtk["Ey",VTKPointData()] = fields.u[2]
                vtk["Hz",VTKPointData()] = fields.u[3]
                vtk["time",VTKFieldData()] = time[1]
            end
        end
    end
end

function main(nt, cfl, n, VTKOutFreq)
    nx = ny = n
    bbox = [-10.0, 10.0, -10.0, 10.0]
    println("main:  nx = ",nx, " ny = ",ny)
    grid = Maxwell.Grid(nx, ny, bbox, cfl)
    println("main grid:  nx = ",grid.nx, " ny = ",grid.ny)

    fields = Maxwell.GridFields(3, grid)
    Maxwell.init_data!(fields)
    
    evolve!(fields, nt, VTKOutFreq)
end

n = 401
nt = 5000
cfl = 0.10
vtkOutFreq = 10

main(nt,cfl,n,vtkOutFreq)


