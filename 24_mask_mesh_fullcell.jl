using Oceananigans
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid, GridFittedBottom, PartialCellBottom
using Printf
using Plots


function show_mask(grid)

    #print("grid = ", grid, "\n")
    c = CenterField(CPU(), grid)
    c .= 1

    mask_immersed_field!(c)

    x, y, z = nodes(c)

    return x, y, z, c
end

### Define grid with a seamount
h0, L = 0.5, 0.25
seamount(x, y, z) = z < - 1 + h0*exp(-y^2/L^2)
grid = RectilinearGrid(size=(16, 8), y=(-1, 1), z=(-1, 0),
                       topology=(Flat, Periodic, Bounded), halo=(3,3)
                       )
grid_with_seamount = ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(seamount_field.data))
#grid_with_seamount = ImmersedBoundaryGrid(grid, GridFittedBoundary(seamount))

### Plot masked region: mask.png
x, y, z, cc = show_mask(grid_with_seamount)
plt = heatmap(y, z, interior(cc)[1,:,:]', xlabel = "y", ylabel = "z", title = "Masked Region_Full_cell")
savefig(plt,"mask_fullcell.png")
