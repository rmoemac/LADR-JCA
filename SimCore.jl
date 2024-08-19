mu = G * M_earth
R_earth = 6371e3
laser_power = 21000

#setting sim initial_conditions
time = 1.5*60*60
t_span = (0,time)

# Function to simulate orbital encounter
function simulate_orbital_encounter(t, y, t_span, LADR1, debris1, laser_output_power)
    # Unpack the state vector
    x_l, y_l, z_l, vx_l, vy_l, vz_l, x_d, y_d, z_d, vx_d, vy_d, vz_d = y

    # Calculate magnitude of the debris position vector
    r_mag = norm([x_d, y_d, z_d])
    
    final_index = length(t_span) - 1
    percent_completion = (t / t_span[final_index]) * 100
    drag_force = 0

    # Unpack the state vector
    ladroit_position, ladroit_velocity = y[1:3], y[4:6]
    debris_position, debris_velocity = y[7:9], y[10:12]
    
    # Calculate gravitational acceleration for L'ADROIT and debris
    ax_l, ay_l, az_l = calculate_gravitational_acceleration(ladroit_position...)
    ax_d, ay_d, az_d = calculate_gravitational_acceleration(debris_position...)
    energy = calculate_orbital_energy(ladroit_velocity, ladroit_position, G * M_earth)
  
    # Calculate the relative distance between L'ADROIT and the debris
    relative_distance = norm([x_l - x_d, y_l - y_d, z_l - z_d])
    
    #println(relative_distance)
  
    if relative_distance < 2150000 && t < 30  # Assuming 250 km as effective laser range
        pulse_energy = LADR1.pulse_energy  # Use small debris pulse energy
        delta_v = calculate_momentum_change(LADR1.pulse_energy, drag_force, LADR1.Cm, debris1.mass, laser_output_power, 1e-10)
        delta_a = calculate_accel_change(LADR1.pulse_energy, LADR1.Cm, debris1.mass, laser_output_power, 1e-10)

        # Apply delta_v in the opposite direction of the debris's velocity vector
        velocity_vector = [vx_d, vy_d, vz_d]
        velocity_magnitude = norm(velocity_vector)
  
        if velocity_magnitude > 0  # To avoid division by zero
            # Normalize velocity vector and apply delta_v in the opposite direction
            normalized_velocity_vector = velocity_vector / velocity_magnitude
            velocity_change = -normalized_velocity_vector * delta_v * 0  # Multiplied by 0 as in original code
            println(delta_v)
            println(velocity_change)
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

#Declare all objects
#debris
sma_deb = 600e3
deb_alt = sma_deb + R_earth
debris1 = Debris(deb_alt,0,98.7,0.0,0.0,170.0)
r_deb = debris1.r
v_deb = debris1.v

#laser
hp_las = 560e3
ha_las = 960e3
sma_laser = (hp_las + ha_las)/2
laser_alt = sma_laser + R_earth
laser1 = LADR(laser_alt,0.028,90.1,-180.1,0,0)
r_laser = laser1.r
v_laser = laser1.v

#put into IC vector
initial_conditions = [r_laser,v_laser,r_deb,v_deb]