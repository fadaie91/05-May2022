ENV["GKSwstype"] = "nul"
using Oceananigans
using Oceananigans.BoundaryConditions: fill_halo_regions!
using Oceananigans.ImmersedBoundaries: ImmersedBoundaryGrid, GridFittedBottom, PartialCellBottom
using Oceananigans.ImmersedBoundaries: mask_immersed_field!
using Printf
using Plots


arch = CPU()

###Testing different advection schemes

tracer_advection = CenteredSecondOrder()
#tracer_advection = UpwindBiasedFirstOrder()
#tracer_advection = UpwindBiasedThirdOrder()
#tracer_advection = CenteredFourthOrder()
#tracer_advection = WENO5()
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

grid = ImmersedBoundaryGrid(underlying_grid, PartialCellBottom(seamount_field.data,minimum_fractional_Δz=0.1))
#grid = ImmersedBoundaryGrid(underlying_grid, GridFittedBottom(seamount_field.data))

# Terrain following coordinate
ζ(y, z) = z / (h(y) - 1)

# Calculate streamfunction
Ψᵢ(x, y, z) = (1 - ζ(y, z))^2
Ψ = Field{Center, Face, Face}(grid)
set!(Ψ, Ψᵢ)
fill_halo_regions!(Ψ, arch)
mask_immersed_field!(Ψ)

# Set velocity field from streamfunction
v = YFaceField(grid)
w = ZFaceField(grid)
v .= + ∂z(Ψ)
w .= - ∂y(Ψ)

fill_halo_regions!(v, arch)
fill_halo_regions!(w, arch)
mask_immersed_field!(v)
mask_immersed_field!(w)

D = compute!(Field(∂y(v) + ∂z(w)))
@info @sprintf("Maximum divergence is %.2e.", maximum(D))

## Set up Model
velocities = PrescribedVelocityFields(; v, w)
    
model = HydrostaticFreeSurfaceModel(; grid,
                                    tracer_advection,
                                    velocities,
                                    tracers = :b,
                                    buoyancy = BuoyancyTracer(),
                                    #coriolis = FPlane(f=1e-4)
                                          coriolis = nothing)

N² = 1
bᵢ(x, y, z) = N² * z
#bᵢ(x, y, z) = 1 + z
set!(model, b = bᵢ)

# Simulation                             
stop_time = 1
Δy = grid.Δyᵃᶜᵃ
@show Δt = 1e-2 * Δy
simulation = Simulation(model; Δt, stop_time)
#simulation = Simulation(model; Δt, stop_iteration=10)

run!(simulation)

b = model.tracers.b
#u, v, w = model.velocities

xb, yb, zb = nodes((Center, Center, Center), grid)
bplot = contourf(yb, zb, interior(b)[1, :, :]', title="Buoyancy", xlabel="y", ylabel="z")
savefig(bplot, "tracer.png")

#xv, yv, zv = nodes((Face, Center, Center), grid)
#vplot = contourf(yv, zv, interior(v)[1, :, :]', title="velocity_V", xlabel="y", ylabel="z")
#savefig(vplot, "Velocity_V.png")
