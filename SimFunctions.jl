using LinearAlgebra

# Define constants
#const G = 6.67430e-11  # Gravitational constant in m^3 kg^-1 s^-2
#const M_earth = 5.972e24  # Mass of the Earth in kg

# Convert degrees to radians
function deg2rad(deg)
    return deg * pi / 180.0
end
function generate_sphere(center, radius, resolution)
    theta = LinRange(0, 2 * pi, resolution)
    phi = LinRange(0, pi, resolution)
    x = [center[1] + radius * sin(p) * cos(t) for p in phi, t in theta]
    y = [center[2] + radius * sin(p) * sin(t) for p in phi, t in theta]
    z = [center[3] + radius * cos(p) for p in phi, t in theta]
    return x, y, z
end
function plot_kepler_variation(sim_output,laser_yn) #1 designates plotting lasers, 0 designates plotting debris
    # Create a 3D plot
    plt = plot()

    # Iterate over each satellite
    num_sats = size(sim_output,1)
    num_timesteps = size(sim_output,2)
    for sat in 1:num_sats
        # Extract x, y, z coordinates for this satellite over time
        x_coords = [sim_output[sat, t][1] for t in 1:num_timesteps]
        y_coords = [sim_output[sat, t][2] for t in 1:num_timesteps]
        z_coords = [sim_output[sat, t][3] for t in 1:num_timesteps]
        
        # Plot the trajectory of each satellite
        if laser_yn == 1
            plot!(x_coords, y_coords, z_coords, label="Satellite $sat")
        else
            plot!(x_coords, y_coords, z_coords, label="debris $sat")
        end
    end

    # Show the plot
    display(plt)
    println("plotting complete")
end


# Function for converting orbital elements to ECI coordinates
function orbital_elements_to_eci(a, e, i_deg, omega_deg, raan_deg, nu_deg)
    i = deg2rad(i_deg)
    omega = deg2rad(omega_deg)
    raan = deg2rad(raan_deg)
    nu = deg2rad(nu_deg)
    
    # Standard gravitational parameter for Earth
    mu = G * M_earth

    # Calculate position and velocity in the perifocal frame
    p = a * (1 - e^2)  # Semi-latus rectum
    r_perifocal = (p / (1 + e * cos(nu))) * [cos(nu), sin(nu), 0]
    v_perifocal = sqrt(mu / p) * [-sin(nu), e + cos(nu), 0]

    # Rotation matrices
    R_raan = [cos(raan) -sin(raan) 0;
              sin(raan) cos(raan) 0;
              0 0 1]
    
    R_i = [1 0 0;
           0 cos(i) -sin(i);
           0 sin(i) cos(i)]
    
    R_omega = [cos(omega) -sin(omega) 0;
               sin(omega) cos(omega) 0;
               0 0 1]
    
    # Transform to ECI frame
    r_eci = R_raan * R_i * R_omega * r_perifocal
    v_eci = R_raan * R_i * R_omega * v_perifocal

    return r_eci, v_eci
end
function generate_orbit_matrix(n::Int, m::Int)
    # Initialize an n x m matrix with zeros
    orbit_matrix = zeros(Float64, n, m)
    
    # You can fill the matrix with actual data 
    # For example, let's assume the value represents a simple linear propagation
    for i in 1:n
        for j in 1:m
            orbit_matrix[i, j] = i * j  # replace with actual propagation logic
        end
    end
    
    return orbit_matrix
end

function remove_every_other_3(arr)
    n = length(arr)
    result = []
    
    # Iterate over the array in steps of 6, keeping the first 3 each time
    for i in 1:6:n
        # Determine the end of the block of three
        end_index = min(i+2, n)
        # Append the first block of three values to the result
        append!(result, arr[i:end_index])
    end
    
    return result
end


# Function for calculating gravitational acceleration
function calculate_gravitational_acceleration(x, y, z)
    r = sqrt(x^2 + y^2 + z^2)
    ax = -G * M_earth * x / r^3
    ay = -G * M_earth * y / r^3
    az = -G * M_earth * z / r^3
    return ax, ay, az
end

function print_progress(current, total)
    percent = (current / total) * 100
    print("\rProgress: ", round(percent, digits=2), "%")
    flush(stdout)  # Ensure the output is flushed
end

# Function for calculating orbital energy
function calculate_orbital_energy(v, r, mu)
    kinetic_energy = 0.5 * norm(v)^2
    potential_energy = r == 0 ? 0 : -mu / norm(r)
    return kinetic_energy + potential_energy
end
function extract_positions(results_yoshi, n_lasers, n_debris)
    n_steps = length(results_yoshi)
    laser_positions = [zeros(3, n_steps) for _ in 1:n_lasers]
    debris_positions = [zeros(3, n_steps) for _ in 1:n_debris]

    for (step_idx, (_, state)) in enumerate(results_yoshi)
        for i in 1:n_lasers
            laser_positions[i][:, step_idx] = state[(i-1)*3+1:(i-1)*3+3]
        end
        for i in 1:n_debris
            debris_positions[i][:, step_idx] = state[(n_lasers+i-1)*3+1:(n_lasers+i-1)*3+3]
        end
    end
    return laser_positions, debris_positions
end
function final_positions_magnitudes(debris_positions)
    n_debris = length(debris_positions)
    final_magnitudes = Float64[]
    
    for i in 1:n_debris
        final_position = debris_positions[i][:, end]
        magnitude = norm(final_position)
        println(magnitude)
        push!(final_magnitudes, magnitude)
    end
    
    return final_magnitudes
end
function extract_final_state(results)
    final_index = length(results)
    x,y,z = results[final_index][2][1],results[final_index][2][2],results[final_index][2][3]
    vx,vy,vz = results[final_index][2][4],results[final_index][2][5],results[final_index][2][6]
    final_pos = [x,y,z]
    final_vel = [vx,vy,vz]
    return final_pos, final_vel

end
function extract_position_single(results_yoshi,n_debris)
    n_steps = length(results_yoshi)
    #laser_positions = [zeros(3, n_steps) for _ in 1:n_lasers]
    debris_positions = [zeros(3, n_steps) for _ in 1:n_debris]

    for (step_idx, (_, state)) in enumerate(results_yoshi)
        # for i in 1:n_lasers
        #     laser_positions[i][:, step_idx] = state[(i-1)*3+1:(i-1)*3+3]
        # end
        for i in 1:n_debris
            debris_positions[i][:, step_idx] = state[(n_debris+i-1)*3+1:(n_debris+i-1)*3+3]
        end
    end
    return debris_positions
end

function extract_velocity_single(results_yoshi,n_debris)
    n_steps = length(results_yoshi)
    #laser_positions = [zeros(3, n_steps) for _ in 1:n_lasers]
    debris_velocities = [zeros(3, n_steps) for _ in 1:n_debris]

    for (step_idx, (_, state)) in enumerate(results_yoshi)
        # for i in 1:n_lasers
        #     laser_positions[i][:, step_idx] = state[(i-1)*3+1:(i-1)*3+3]
        # end
        for i in 1:n_debris
            debris_velocities[i][:, step_idx] = state[(n_debris+i-1)*3+4:(n_debris+i-1)*3+6]
        end
    end
    return debris_velocities
end
function final_position(debris_positions)
    n_debris = length(debris_positions)
    final_positions = Float64[]
    
    for i in 1:n_debris
        final_position = debris_positions[i][:, end]
        push!(final_positions, final_position)
    end
    
    return final_positions
end
function final_velocity(debris_velocities)
    n_debris = length(debris_positions)
    final_velocities = Float64[]
    
    for i in 1:n_debris
        final_velocity = debris_velocities[i][:, end]
        push!(final_velocities, final_velocity)
    end
    
    return final_velocities
end


# Function for calculate_momentum_change and delta_accel (placeholders, need proper implementation)
function calculate_momentum_change(Cm, target_mass, laser_output_power)
    laser_output_power = laser_output_power * 1e-6
    newtons = Cm * laser_output_power
    dt = 1e-6
    delta_momentum = newtons * dt
    
    delta_velocity = delta_momentum / target_mass
    return delta_velocity
end

function calculate_accel_change(Cm, target_mass, laser_output_power)
    newtons = Cm * (laser_output_power * 1e-6)
    thrust = newtons
    delta_accel = thrust / target_mass
    return delta_accel
end

function calculate_orbital_elements(r, v, mu)
    r_mag = norm(r)
    v_mag = norm(v)
    
    # Specific angular momentum
    h = cross(r, v)
    
    # Eccentricity vector
    e_vec = ((v_mag^2 - mu / r_mag) * r - dot(r, v) * v) / mu
    e = norm(e_vec)
    
    # Specific orbital energy
    epsilon = (v_mag^2 / 2) - mu / r_mag
    
    # Semi-major axis
    a = -mu / (2 * epsilon)
    
    return a, e, e_vec, h
end
function calculate_distance(pos1::Vector{Float64}, pos2::Vector{Float64})
    return norm(pos1 - pos2)
end
function sort_debris_by_distance(laser, debris_list)
    laser_position = laser.r

    # Compute distances and create a tuple of (debris, distance)
    distances = [(debris, calculate_distance(laser_position, debris.r)) for debris in debris_list]

    # Sort the tuple array based on distance (second element of each tuple)
    sorted_debris = sort(distances, by=x -> x[2])

    # Return only the sorted debris objects
    return [x[1] for x in sorted_debris]
end


