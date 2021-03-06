using Oceananigans
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid, GridFittedBottom, PartialCellBottom
using Printf
using GLMakie
using Oceananigans.ImmersedBoundaries: mask_immersed_field!


arch = CPU()
#tracer_advection = CenteredSecondOrder()
#momentum_advection = CenteredSecondOrder()
momentum_advection = WENO5()
tracer_advection = WENO5()
#momentum_advection = CenteredFourthOrder()
#tracer_advection = CenteredFourthOrder()
#momentum_advection = UpwindBiasedFirstOrder()
#tracer_advection = UpwindBiasedFirstOrder()
#momentum_advection = UpwindBiasedThirdOrder()
#tracer_advection = UpwindBiasedThirdOrder()
#momentum_advection = UpwindBiasedFifthOrder()
#tracer_advection = UpwindBiasedFifthOrder()


underlying_grid = RectilinearGrid(arch,
                                  size=(128, 64), halo=(3, 3), 
                                  y = (-1, 1),
                                  z = (-1, 0),
                                  topology=(Flat, Periodic, Bounded))

# A bump
h₀ = 0.5 # bump height
L = 0.25 # bump width
@inline h(y) = h₀ * exp(- y^2 / L^2)
@inline seamount(x, y) = - 1 + h(y)

seamount_field = Field{Center, Center, Nothing}(underlying_grid)
set!(seamount_field, seamount)
fill_halo_regions!(seamount_field)

minimum_fractional_Δz = 1

immersed_boundaries = [
                       PartialCellBottom(seamount_field.data;
                                         minimum_fractional_Δz),
                       GridFittedBottom(seamount_field.data)
                      ]

b = []
v = []

function progress(sim)
    vmax = maximum(abs, sim.model.velocities.v)
    @info @sprintf("Iter: %d, time: %.2e, max|v|: %.2e",
                   iteration(sim), time(sim), vmax)

    return nothing
end

for ib in immersed_boundaries
    grid = ImmersedBoundaryGrid(underlying_grid, ib)
    #grid = ImmersedBoundaryGrid(underlying_grid, PartialCellBottom(seamount_field.data))


    @show grid

    





### In this example of the paper the velocities were prescribe
### which means that they won't evolve during the simulation
### and using the function below we prescribe the velocities
#velocities = PrescribedVelocityFields( v=V, w=W)
  
  model = HydrostaticFreeSurfaceModel(; grid,
                                        tracer_advection,
                                        momentum_advection,
                                        #coriolis = FPlane(f=0.1),
                                        tracers = :b,
                                        #velocities = velocities,
                                        buoyancy = BuoyancyTracer())

B= Field{Center, Center, Center}(grid)   
N² = 1
#bᵢ(x, y, z) = N² * z
#B(x, y, z) = N² + z                                   
set!(B, (x, y, z) -> (N² + z) ) 
mask_immersed_field!(B)

U, V, W = model.velocities

Ψ = Field{Center, Face, Face}(grid)

h(y)    = h₀*exp(-y^2/L^2)
ζ(y, z) = z/(h(y) - 1)
set!(Ψ, (x, y, z) -> (1 - ζ(y, z))^2)
fill_halo_regions!(Ψ, arch)


### We mask psi
mask_immersed_field!(Ψ)

### V is (Center, Face, Center)
V = YFaceField(grid)
### W is (Center, Centere, Face)
W = ZFaceField(grid)

### V and W are deravatives of psi 
### which in this code we use exact expression 
### which we computed in maple
### we are calling this method 'analytical'
V.=  ∂z(Ψ)
W.= -∂y(Ψ)

### We mask V and W        
mask_immersed_field!(V)
mask_immersed_field!(W)

### We fill the halo regions of V and W

fill_halo_regions!(V, arch)
fill_halo_regions!(W, arch)


  
    set!(model, b = B,  w=W, v=V)

    simulation = Simulation(model; Δt=1e-4, stop_time=0.2)
    simulation.callbacks[:p] = Callback(progress, IterationInterval(10))


    run!(simulation)

    push!(b, Array(interior(model.tracers.b, 1, :, :)))
    push!(v, Array(interior(model.velocities.v, 1, :, :)))
end

b_partial = b[1]
b_full    = b[2]
Δb = b_full .- b_partial

v_partial = v[1]
v_full    = v[2]
Δv = v_full .- v_partial

fig = Figure(resolution=(1200, 1800))

partial_cell_title = @sprintf("PartialCellBottom with ϵ = %.1f", minimum_fractional_Δz)
ax_bp = Axis(fig[1, 2], title=partial_cell_title)
ax_bf = Axis(fig[2, 2], title="GridFittedBottom")
ax_bd = Axis(fig[3, 2], title="Difference (GridFitted - PartialCell)")

# ax_vp = Axis(fig[1, 3])
# ax_vf = Axis(fig[2, 3])
# ax_vd = Axis(fig[3, 3])

color = (:black, 0.5)
linewidth = 3
levels = 15

hmbp = heatmap!(ax_bp, b_partial)
contour!(ax_bp, b_partial; levels, color, linewidth)
Colorbar(fig[1, 1], hmbp, label="Buoyancy", flipaxis=false)

hmbf = heatmap!(ax_bf, b_full)
contour!(ax_bf, b_full; levels, color, linewidth)
Colorbar(fig[2, 1], hmbf, label="Buoyancy", flipaxis=false)

hmbd = heatmap!(ax_bd, Δb)
Colorbar(fig[3, 1], hmbd, label="Buoyancy", flipaxis=false)

#=
hmvp = heatmap!(ax_vp, v_partial)
#contour!(ax_vp, v_partial, levels=15)
#Colorbar(fig[1, 1], hmvp)
hmvf = heatmap!(ax_vf, v_full)
#contour!(ax_vf, v_full, levels=15)
#Colorbar(fig[2, 1], hmvf)
hmvd = heatmap!(ax_vd, Δv)
#Colorbar(fig[3, 1], hmvd)
=#

display(fig)
