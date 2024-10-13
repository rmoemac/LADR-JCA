include("SimCore.jl")
include("SimObjects.jl")

G = 6.67430 * 10^(-11)
M_earth = 5.972e24
mu = G * M_earth
R_earth = 6371e3
laser_power = 21000.0

int_method = 2 #1 is mdpt, 2 is yoshida, 3 is rk4

#setting sim initial_conditions
time = 10*60*60
dt = 1.0/56.0
t_span = (0.0,time)
#place holder until scaling happens
n_lasers = 10
n_debris = 30
n_obj = n_lasers + n_debris

#setting laser spec
laser_range = 900000.0
cm = 99.0

laser_list = []
debris_list = []


using Evolutionary
using Random

# Generate the debris objects
for i in 1:n_debris
    altitude = rand() * (1100000.0 - 900000.0) + 900000.0 + R_earth
    inclination = rand() * 90.0
    TA = rand() * 360.0
    e = rand() * 0.2
    raan_deg = rand() * 360.0
    omega_deg = rand() * 360.0
    mass = 1600.0 # kg
    deb = Debris(altitude, e, inclination, omega_deg, raan_deg, TA, mass)
    push!(debris_list,deb)
end

baseline(debris_list)

Random.seed!(42) # Set seed for reproducibility

min_alt = 900000 #meters
max_alt = 110000 #meters
min_inc = 0 #degrees
max_inc = 90 #degrees
min_tracks = 1
max_tracks = 5
min_along_track = 1
max_along_track = 10

baseline_alts = baseline(debris_list)



result = Evolutionary.optimize(simulator, [min_alt,max_alt], [min_inc, max_inc], [min_tracks, max_tracks], [min_along_track,max_along_track], # Boundary for the parameter space
                               Evolutionary.SGA(; selection=Evolutionary.RouletteWheelSelection(),
                                                crossover=Evolutionary.BLXAlphaCrossover(0.5),
                                                mutation=Evolutionary.GaussianMutation(0.1),
                                                μ=50, λ=100, generations=200))

function cost_function()

end
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
