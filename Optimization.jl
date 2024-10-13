include("SimCore.jl")
include("SimObjects.jl")

G = 6.67430 * 10^(-11)
M_earth = 5.972e24
mu = G * M_earth
R_earth = 6371e3
laser_power = 0.0 #21000.0

int_method = 2 #1 is mdpt, 2 is yoshida, 3 is rk4

#setting sim initial_conditions
time = 10*60*60
dt = 1.0/56.0
t_span = (0.0,time)
#place holder until scaling happens
n_lasers = 10
n_debris = 10
n_obj = n_lasers + n_debris

#setting laser spec
laser_range = 900000.0
cm = 99.0

laser_list = []
debris_list = []


using Evolutionary
using Random




#baseline(debris_list)
debris_altitudes = []

Random.seed!(42) # Set seed for reproducibility

min_alt = 900000 #meters
max_alt = 110000 #meters
min_inc = 0 #degrees
max_inc = 90 #degrees
min_tracks = 1
max_tracks = 5
min_along_track = 1
max_along_track = 10

for i in 1:n_lasers
    altitude = rand() * (1100000.0 - 900000.0) + 900000.0 + R_earth
    println(altitude)
    inclination = rand() * 90.0
    TA = rand() * 360.0
    e = 0.0
    raan_deg = rand() * 360.0
    omega_deg = rand() * 360.0
    power = laser_power
    laz = LADR(altitude, e, inclination, omega_deg, raan_deg, TA, power)
    push!(laser_list, laz)
end
# Generate the debris objects
for i in 1:n_debris
    altitude = rand() * (1100000.0 - 900000.0) + 900000.0 + R_earth
    inclination = rand() * 90.0
    TA = rand() * 360.0
    e = 0.0
    raan_deg = rand() * 360.0
    omega_deg = rand() * 360.0
    mass = 1600.0 # kg
    deb = Debris(altitude, e, inclination, omega_deg, raan_deg, TA, mass)
    push!(debris_list,deb)
end

#baselined based off initial conditions altitude
for i in 1:length(debris_list)
    push!(debris_altitudes,norm(debris_list[i].r))
end
# Function to create a 3D plot of orbit matrix and Earth
using Plots  # Only use Plots.jl, no need for PlotlyJS directly
plotlyjs()   # Set the backend to PlotlyJS for interactive 3D plots

# Function to create a 3D plot of orbit matrix and Earth
function plot_orbit_and_earth(orbit_matrix::Matrix{Float64})
    # Extract dimensions of the orbit matrix
    n, m = size(orbit_matrix)
    
    # Extract coordinates for 3D plotting (using orbit matrix values as z-coordinates)
    x = repeat(collect(1:n), outer = m)  # X-axis: Satellite position or time
    y = repeat(reshape(collect(1:m), m, 1), n)  # Y-axis: Reshape to match dimensions, no transpose needed
    z = vec(orbit_matrix)                 # Z-axis: The matrix values represent altitude or position

    # Plot the orbit data as a 3D scatter plot
    p = Plots.scatter3d(x, y, z, label = "Orbit Points", marker = :circle, markersize = 5, color=:blue)

    # Plot the Earth as a sphere
    # Parameters for the Earth (let's assume radius = 1 for simplicity)
    sphere_radius = 1.0
    θ = LinRange(0, 2π, 50)  # Theta (longitude)
    φ = LinRange(0, π, 50)   # Phi (latitude)

    # Create a 3D grid for the sphere (Earth)
    X = [sphere_radius * sin(φ[i]) * cos(θ[j]) for i in 1:length(φ), j in 1:length(θ)]
    Y = [sphere_radius * sin(φ[i]) * sin(θ[j]) for i in 1:length(φ), j in 1:length(θ)]
    Z = [sphere_radius * cos(φ[i]) for i in 1:length(φ), j in 1:length(θ)]

    # Plot the Earth sphere
    plot3d!(X, Y, Z, st = :surface, color = :green, opacity = 0.5, label = "Earth")
    
    # Customize plot appearance
    xlabel!("X-axis")
    ylabel!("Y-axis")
    zlabel!("Z-axis")
    title!("3D Orbit and Earth Plot")

    display(p)
end

@time begin

laser_pos = sim_wrapper(laser_list,debris_list,5,0)
sim_test_kepler = sim_wrapper(laser_list,debris_list,4,0)
end

println("plotting started")

plot_kepler_variation(sim_test_kepler,0)
plot_kepler_variation(laser_pos,1)
# Create a 3D plot
plt = plot()

sim_alt_hybrid = []
for i in 1:n_debris
    final_pos = sim_test_kepler[i,size(sim_test_kepler,2)]
    push!(sim_alt_hybrid,norm(final_pos))
end

time_event_error = sim_alt_hybrid[1] - debris_altitudes[1]
perc_error_time = time_event_error/debris_altitudes[1]
println("Time event error over $time seconds is $time_event_error and $perc_error_time %")




# #RK4
# sim_alt_rk4 = sim_wrapper(laser_list,debris_list,1,0)
# #mdpt
# sim_alt_mdpt = sim_wrapper(laser_list,debris_list,2,0)
# #yoshida
# sim_alt_yoshi = sim_wrapper(laser_list,debris_list,3,0)


# rk4_error = sim_alt_rk4 - debris_altitudes
# mdpt_error = sim_alt_mdpt - debris_altitudes
# yoshi_error = sim_alt_yoshi - debris_altitudes 

# perc_error_rk4 = rk4_error/sim_alt_rk4 * 100
# perc_error_mdpt = mdpt_error/sim_alt_mdpt * 100
# perc_error_yoshi = yoshi_error/sim_alt_yoshi * 100

# println("RK4 Absolute error over $time seconds is $rk4_error and $perc_error_rk4 %")
# println("Midpoint Absolute error over $time seconds is $mdpt_error and $perc_error_mdpt %")
# println("Yoshida Absolute error over $time seconds is $yoshi_error and $perc_error_yoshi %")

function simulator(x)
    alt_prime = x[1]
    inc_prime = x[2]
    n_tracks = x[3]
    n_along_tracks = x[4]


    for i in 1:n_tracks
        for j in 1:n_along_tracks
            altitude = alt_prime
            inclination = inc_prime + 180/n_tracks
            TA = 360 / n_along_tracks
            e = 0.0
            raan_deg = rand() * 360.0
            omega_deg = rand() * 360.0
            power = laser_power
            laz = LADR(altitude, e, inclination, omega_deg, raan_deg, TA, power)
            push!(laser_list, laz)
        end
    end


    int_method = 2

    sim_outcome = sim_wrapper(laser_list,debris_list,int_method)
    return sim_outcome
end

function baseline(debris_list)
    n_lasers = 1
    for i in 1:n_lasers
        altitude = rand() * (1100000.0 - 900000.0) + 900000.0 + R_earth
        inclination = rand() * 90.0
        TA = rand() * 360.0
        e = 0.0
        raan_deg = rand() * 360.0
        omega_deg = rand() * 360.0
        power = laser_power
        laz = LADR(altitude, e, inclination, omega_deg, raan_deg, TA, power)
        push!(laser_list, laz)
    end

    baseline_results = sim_wrapper(laser_list,debris_list,int_method)
    las_pos,debris_pos = extract_positions(baseline_results)
    final_altitudes = final_positions_magnitude(debris_pos)

    return final_altitudes

end
