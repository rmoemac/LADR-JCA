include("SimObjects.jl")
include("SimFunctions.jl")
include("Rk4-Integrator.jl")
include("YoshIntegrator.jl")

using Base.Threads


G = 6.67430 * 10^(-11)
M_earth = 5.972e24
mu = G * M_earth
R_earth = 6371e3
laser_power = 21000.0

#setting sim initial_conditions
time = 10*60*60
t_span = (0.0,time)
dt = 1.0/56.0
#place holder until scaling happens
n_lasers = 10
n_debris = 10
n_obj = n_lasers + n_debris

#setting laser spec
laser_range = 900000.0
cm = 99.0

laser_list = []
debris_list = []

#generate laser platforms
# Generate laser platforms
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


# Function to simulate orbital encounter
# t_span, LADR1, debris1, laser_output_power
function simulate_orbital_encounter(t, y)
    # Unpack the state vector
    x_l, y_l, z_l, vx_l, vy_l, vz_l, x_d, y_d, z_d, vx_d, vy_d, vz_d = y

    # Calculate magnitude of the debris position vector
    r_mag = norm([x_d, y_d, z_d])

    # Unpack the state vector
    ladroit_position, ladroit_velocity = y[1:3], y[4:6]
    debris_position, debris_velocity = y[7:9], y[10:12]
    
    # Calculate gravitational acceleration for L'ADROIT and debris
    ax_l, ay_l, az_l = calculate_gravitational_acceleration(ladroit_position[1],ladroit_position[2],ladroit_position[3])
    ax_d, ay_d, az_d = calculate_gravitational_acceleration(debris_position[1],debris_position[2],debris_position[3])
    energy = calculate_orbital_energy(ladroit_velocity, ladroit_position, G * M_earth)
  
    # Calculate the relative distance between L'ADROIT and the debris
    relative_distance = norm([x_l - x_d, y_l - y_d, z_l - z_d])
    
    #println(relative_distance)
  
    if relative_distance < 2150000 && t < 30  # Assuming 250 km as effective laser range
          # Use small debris pulse energy
        delta_v = calculate_momentum_change(99.0, debris1.mass, LADR1.power)
        delta_a = calculate_accel_change(99.0, debris1.mass, LADR1.power)

        # Apply delta_v in the opposite direction of the debris's velocity vector
        velocity_vector = [vx_d, vy_d, vz_d]
        velocity_magnitude = norm(velocity_vector)
  
        if velocity_magnitude > 0  # To avoid division by zero
            # Normalize velocity vector and apply delta_v in the opposite direction
            normalized_velocity_vector = velocity_vector / velocity_magnitude
            velocity_change = -normalized_velocity_vector * delta_v * 0  # Multiplied by 0 as in original code
            vx_d += velocity_change[1]
            vy_d += velocity_change[2]
            vz_d += velocity_change[3]

            acceleration_vector = [ax_d, ay_d, az_d]
            acceleration_magnitude = norm(acceleration_vector)
            normalized_acceleration_vector = acceleration_vector / acceleration_magnitude
            acceleration_change = -normalized_velocity_vector * delta_a
            ax_d += acceleration_change[1]
            ay_d += acceleration_change[2]
            az_d += acceleration_change[3]
        end
    end

    return [vx_l, vy_l, vz_l, ax_l, ay_l, az_l, vx_d, vy_d, vz_d, ax_d, ay_d, az_d]
end

function yoshi_orbital_encounter(t,y, v)

    #parse inputs to set as laser and debris IC for the return
    for i in 0:(n_obj-1)
        if i < n_lasers
            # Lasers
            #println("Processing laser index: $i")
            laser_list[i+1].r[1] = y[i*3 + 1]
            laser_list[i+1].r[2] = y[i*3 + 2]
            laser_list[i+1].r[3] = y[i*3 + 3]
        else
            # Debris
            debris_index = i - n_lasers
            #println("Processing debris index: $debris_index")
            debris_list[debris_index+1].r[1] = y[i*3 + 1]
            debris_list[debris_index+1].r[2] = y[i*3 + 2]
            debris_list[debris_index+1].r[3] = y[i*3 + 3]
        end
    end

    #Calculate gravitational acceleration for Laser platforms and debris
    for i in 1:n_lasers
        x = laser_list[i].r[1]
        y = laser_list[i].r[2]
        z = laser_list[i].r[3]
        ax,ay,az = calculate_gravitational_acceleration(x,y,z)
        laser_list[i].accel = [ax,ay,az]
    
    end
    for i in 1:n_debris
        x = debris_list[i].r[1]
        y = debris_list[i].r[2]
        z = debris_list[i].r[3]
        ax,ay,az = calculate_gravitational_acceleration(x,y,z)
        debris_list[i].accel = [ax,ay,az]
    end



    # for i in 1:n_lasers
    #     #sorting debris based on distance
    #     sorted_debris= sort_debris_by_distance(laser_list[i],debris_list)
    #     # Calculate the relative distance between L'ADROIT and the debris
    #     for j in 1:n_debris
    #         relative_distance = norm(laser_list[i].r - debris_list[j].r)

    #         if relative_distance < laser_range  # Assuming 250 km as effective laser range
    #             # Use small debris pulse energy
    #             delta_a = calculate_accel_change(cm, debris_list[j].mass, laser_list[i].power)

    #             # Apply delta_v in the opposite direction of the debris's velocity vector
    #             velocity_vector = debris_list[j].v
    #             velocity_magnitude = norm(velocity_vector)
        
    #             if velocity_magnitude > 0   # To avoid division by zero
    #                 # Normalize velocity vector and apply delta_v in the opposite direction
    #                 normalized_velocity_vector = velocity_vector / velocity_magnitude

    #                 acceleration_change = -normalized_velocity_vector * delta_a * 0
    #                 debris_list[j].r[1] += acceleration_change[1]
    #                 debris_list[j].r[2] += acceleration_change[2]
    #                 debris_list[j].r[3] += acceleration_change[3]
    #             end
    #         else
    #             break
    #         end
    #     end
    # end

    ret = zeros(n_obj * 3)  # This will store the accelerations

    for i in 0:(n_obj-1)
        if i < n_lasers
            # Lasers
            #println("Processing laser index: $i")
            ret[i*3 + 1] = laser_list[i+1].accel[1]
            ret[i*3 + 2] = laser_list[i+1].accel[2]
            ret[i*3 + 3] = laser_list[i+1].accel[3]
        else
            # Debris
            debris_index = i - n_lasers + 1
            #println("Processing debris index: $debris_index")
            ret[i*3 + 1] = debris_list[debris_index].accel[1]
            ret[i*3 + 2] = debris_list[debris_index].accel[2]
            ret[i*3 + 3] = debris_list[debris_index].accel[3]
        end
    end
    
    #println("Acceleration vector (ret): ", ret)

    #print(ret)
    
    return ret
end

#Declare all objects
#debris
# sma_deb = 600e3
# deb_alt = sma_deb + R_earth
# debris1 = Debris(deb_alt,0.0,98.7,0.0,0.0,170.0,1600.0)
# r_deb = debris1.r
# v_deb = debris1.v

# #laser
# hp_las = 560e3
# ha_las = 960e3
# sma_laser = (hp_las + ha_las)/2
# laser_alt = sma_laser + R_earth
# LADR1 = LADR(laser_alt,0.028,90.1,-180.1,0.0,0.0,21000.0)
# r_laser = LADR1.r
# v_laser = LADR1.v'

initial_conditions = []

for i in 1:n_lasers
    pos = laser_list[i].r
    vel = laser_list[i].v
    push!(initial_conditions,pos[1])
    push!(initial_conditions,pos[2])
    push!(initial_conditions,pos[3])
    push!(initial_conditions,vel[1])
    push!(initial_conditions,vel[2])
    push!(initial_conditions,vel[3])
end
for i in 1:n_debris
    pos = debris_list[i].r
    vel = debris_list[i].v
    push!(initial_conditions,pos[1])
    push!(initial_conditions,pos[2])
    push!(initial_conditions,pos[3])
    push!(initial_conditions,vel[1])
    push!(initial_conditions,vel[2])
    push!(initial_conditions,vel[3])
end




#put into IC vector
#initial_conditions = [r_laser;v_laser;r_deb;v_deb]

#result = solve_ivp(simulate_orbital_encounter,t_span,initial_conditions,method="RK45",atol = 0.1,max_step = 1/56)
#results = rk4(t_span,initial_conditions,dt)

#Initialize accelerations
for i in 1:n_lasers
    x = laser_list[i].r[1]
    y = laser_list[i].r[2]
    z = laser_list[i].r[3]
    ax,ay,az = calculate_gravitational_acceleration(x,y,z)
    laser_list[i].accel = [ax,ay,az]

end
for i in 1:n_debris
    x = debris_list[i].r[1]
    y = debris_list[i].r[2]
    z = debris_list[i].r[3]
    ax,ay,az = calculate_gravitational_acceleration(x,y,z)
    debris_list[i].accel = [ax,ay,az]
end

@time begin
    results_yoshi = yoshida(t_span, initial_conditions, dt)
end
println("integration complete")


# Extract the positions over time from the integration results
function extract_positions(results_yoshi, n_lasers, n_debris)
    n_steps = length(results_yoshi)
    laser_positions = [zeros(3, n_steps) for _ in 1:n_lasers]
    debris_positions = [zeros(3, n_steps) for _ in 1:n_debris]

    for (step_idx, (_, state)) in enumerate(results_yoshi)
        for i in 1:n_lasers
            laser_positions[i][:, step_idx] = state[(i-1)*6+1:(i-1)*6+3]
        end
        for i in 1:n_debris
            debris_positions[i][:, step_idx] = state[(n_lasers+i-1)*6+1:(n_lasers+i-1)*6+3]
        end
    end
    return laser_positions, debris_positions
end

# Extract positions from the results
laser_positions, debris_positions = extract_positions(results_yoshi, n_lasers, n_debris)

# Generate Earth sphere for visualization
n = 100  # Resolution of the sphere
x_sphere, y_sphere, z_sphere = generate_sphere([0, 0, 0], R_earth, n)

#Plotting the trajectories with Plotly
plotlyjs()
p = plot3d(
    laser_positions[1][1, :], laser_positions[1][2, :], laser_positions[1][3, :],
    seriestype = :scatter,
    markersize = 1,
    label = "Laser 1",
    xlabel = "X (m)",
    ylabel = "Y (m)",
    zlabel = "Z (m)",
    title = "Orbital Trajectories"
)

for i in 2:n_lasers
    plot3d!(laser_positions[i][1, :], laser_positions[i][2, :], laser_positions[i][3, :], label = "Laser $i")
end

for i in 1:n_debris
    plot3d!(debris_positions[i][1, :], debris_positions[i][2, :], debris_positions[i][3, :], label = "Debris $i")
end

# Add the Earth sphere to the plot
surface!(p, x_sphere, y_sphere, z_sphere, opacity = 0.3, color = :blue, label = "Earth", colorbar = false)

# Add starting positions as markers (convert to arrays)
for i in 1:n_lasers
    scatter3d!([laser_positions[i][1, 1]], [laser_positions[i][2, 1]], [laser_positions[i][3, 1]], markercolor=:red, markershape=:cross, markersize=3, label="Laser Start $i")
end

for i in 1:n_debris
    scatter3d!([debris_positions[i][1, 1]], [debris_positions[i][2, 1]], [debris_positions[i][3, 1]], markercolor=:green, markershape=:cross, markersize=3, label="Debris Start $i")
end

# Display the plot
display(p)