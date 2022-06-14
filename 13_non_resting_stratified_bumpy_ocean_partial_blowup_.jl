using Printf
using Oceananigans
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid, GridFittedBottom, PartialCellBottom
#using Oceananigans.ImmersedBoundaries: solid_node
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Plots
using Oceananigans.ImmersedBoundaries: mask_immersed_field!

arch = CPU()

function show_mask(grid)

    print("grid = ", grid, "\n")
    c = CenterField(CPU(), grid)
    c .= 1

    mask_immersed_field!(c)

    x, y, z = nodes(c)

    return x, y, z, c
end

grid =  RectilinearGrid(arch,
                                  size=(160, 20), halo=(3, 3), 
                                  y = (-4, 4),
                                  z = (-1, 0),
                                  topology=(Flat, Periodic, Bounded))

# Gaussian seamount
h₀ = 0.1 # bump height
L = 1 # bump width
@inline h(y) = h₀ * exp(- y^2 / L^2)
@inline seamount(x, y) = - 1 + h(y)

seamount_field = Field{Center, Center, Nothing}(grid)
#grid_with_seamount = ImmersedBoundaryGrid(grid, GridFittedBoundary(seamount))
grid_with_seamount = ImmersedBoundaryGrid(grid, GridFittedBottom(seamount_field.data))
#grid = ImmersedBoundaryGrid(underlying_grid, PartialCellBottom(seamount_field.data,minimum_fractional_Δz=0.2))


###

#Ny = grid_with_seamount.data.Ny
Ny = grid_with_seamount.Ny

topography_index = zeros(Int64, Ny)
for (iyC, yC) in enumerate(grid_with_seamount.yC[1:Ny])
    print("y = ", yC, " z = ", grid_with_seamount.zC[1], " ")

    #izC = 1
    #test = true
    #while test && izC <= grid_with_seamount.Nz
    #	  print(izC, " ")
	 # test =  immersed_cell(1, iyC, izC, grid_with_seamount)
	  #izC += 1
    #end
    #print(" ")
    
    #print(immersed_cell(1, iyC, 1, grid_with_seamount)," ")
    #print(immersed_cell(1, iyC, 2, grid_with_seamount),"\n")
    #print(solid_node(Center(), Center(), Center(), 1,iyC,1, grid_with_seamount)," ")
    #print(solid_node(Center(), Center(), Center(), 1,iyC,2, grid_with_seamount),"\n")
    #topography_index[iyC] = izC - 1 
end

print(topography_index, "\n")

# suppose we have a tracer and we are calling it Theta

Theta_topography_slice = zeros(Ny)
for (iyC, yC) in enumerate(grid_with_seamount.grid.yC[1:Ny])
    Theta_topography_slice[iyC] = Theta[1, iyC, topography_index[iyC]]
end

plt = plot(yC[1:Ny], Theta_topography_slice)
display(plt)
