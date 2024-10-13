include("SimObjects.jl")
include("SimFunctions.jl")
include("Rk4-Integrator.jl")
include("YoshIntegrator.jl")
include("MidpointIntegrator.jl")
include("KeplersAlgorithm.jl")

using Base.Threads

using Plots
#using PlotlyJS

function single_propogator(laser_list,object_conditions,int_method,laser_power)


end

function sim_wrapper(laser_list,debris_list,int_method,laser_power)
    n_lasers = length(laser_list)
    n_debris = length(debris_list)
    n_obj = n_debris + n_lasers

    initial_conditions = zeros(n_obj*6)

    for i in 0:n_lasers-1
        pos = laser_list[i+1].r
        vel = laser_list[i+1].v
        initial_conditions[i*6 + 1] = pos[1]
        initial_conditions[i*6 + 2] = pos[2]
        initial_conditions[i*6 + 3] = pos[3]
        initial_conditions[i*6 + 4] = vel[1]
        initial_conditions[i*6 + 5] = vel[2]
        initial_conditions[i*6 + 6] = vel[3]
        # push!(initial_conditions,pos[1])
        # push!(initial_conditions,pos[2])
        # push!(initial_conditions,pos[3])
        # push!(initial_conditions,vel[1])
        # push!(initial_conditions,vel[2])
        # push!(initial_conditions,vel[3])
    end
    for i in n_lasers:n_lasers+n_debris-1
        pos = debris_list[i+1-n_lasers].r
        vel = debris_list[i+1-n_lasers].v
        initial_conditions[i*6 + 1] = pos[1]
        initial_conditions[i*6 + 2] = pos[2]
        initial_conditions[i*6 + 3] = pos[3]
        initial_conditions[i*6 + 4] = vel[1]
        initial_conditions[i*6 + 5] = vel[2]
        initial_conditions[i*6 + 6] = vel[3]
    end


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

    #results_yoshi = yoshida(t_span,initial_conditions,dt)
    if int_method == 1
        @time begin
            results_rk4 = rk4(t_span,initial_conditions,dt)
        end
        #println("plotting started")
        laser_pos, debris_pos = extract_positions(results_rk4,n_lasers,n_debris)
        #println(laser_pos)
        final_altitudes = final_positions_magnitudes(laser_pos)
        return final_altitudes
    end
    if int_method == 2
        @time begin
            results_mdpt = midpoint(t_span,initial_conditions,dt)
        end
        laser_pos, debris_pos = extract_positions(results_mdpt,n_lasers,n_debris)
        #println(results_mdpt)
        final_altitudes = final_positions_magnitudes(laser_pos)
        return final_altitudes
        # println("plotting started")
        # plot_mdpt()

    end
    if int_method == 3
        @time begin
            results_yoshi = yoshida(t_span,initial_conditions,dt,laser_power)
            #plot_yoshida()
        end
        laser_pos, debris_pos = extract_positions(results_yoshi,n_lasers,n_debris)
        final_altitudes = final_positions_magnitudes(laser_pos)
        #plot_yoshida()
        return final_altitudes
        # println("plotting started")
        # plot_yoshida()
    end
    if int_method == 4
        kepler_counter = 0.0
        int_counter = 0.0
        n_debris = length(debris_list)
        @time begin
            laser_positions = Array{Vector{Float64},2}(undef,n_lasers,Int(t_span[2]))
            #println(laser_positions)
            for i in 1:n_lasers
                kepler_method!(t_span[2],laser_list[i],laser_positions,i)           
            end
            debris_results = zeros(n_debris)
            debris_positions = Array{Vector{Float64},2}(undef,n_debris,Int(t_span[2] / 10.0))
            for i in 1:n_debris
                debris_positions[i,1] = [debris_list[i].r[1],debris_list[i].r[2],debris_list[i].r[3]]
            end
            @threads for i in 1:n_debris
                laser_counter = 0
                t_prop = 10
                for t in 2:Int(t_span[2]/10)
                    #println(t)
                    laser_counter = 0
                    if t_span[2] - t < t_prop
                        debris_positions[i,t] = debris_positions[i,t-1]
                        break
                    end

                    for j in 1:n_lasers
                        rel_range = norm(debris_positions[i,t-1] - laser_positions[j,t])
                        if rel_range < laser_range
                            laser_counter +=1
                        end
                            
                    end
                    if rand() < 1
                        laser_counter = 1
                    else
                        laser_counter = 0
                    end

                    initial_conditions = zeros(6)
                    #t_prop = 10.0 #time to propogate before checking again
                    if laser_counter > 0

                        pos = debris_list[i].r
                        vel = debris_list[i].v

                        initial_conditions[1] = pos[1]
                        initial_conditions[2] = pos[2]
                        initial_conditions[3] = pos[3]
                        initial_conditions[4] = vel[1]
                        initial_conditions[5] = vel[2]
                        initial_conditions[6] = vel[3]

                        t_interval = (0.0,Float64(t_prop))
                        dt = 1.0/56.0
                        temp_deb_results = midpoint(t_interval,initial_conditions,dt,laser_power)
                        final_position,final_velocity = extract_final_state(temp_deb_results)
                        # println(length(temp_deb_results))
                        # println(temp_deb_results[length(temp_deb_results)][2])
                        # pos_log = extract_position_single(temp_deb_results,n_debris)
                        # vel_log = extract_velocity_single(temp_deb_results,n_debris)

                        # final_position = final_position(pos_log)
                        # final_velocity = final_velocity(vel_log)
                        
                        #redefining object list
                        debris_list[i].r = final_position
                        debris_list[i].v = final_velocity
                        
                        debris_positions[i,t] = final_position
                        int_counter += 1
                        #println(debris_positions)

                        
                    else
                        original_position = debris_positions[i,t]
                        debris_positions[i,t] = kepler_method_event(t,t_prop,debris_list[i],debris_positions,i)
                        debris_list[i].r = debris_positions[i,t]
                        delta_pos = debris_list[i].r - original_position
                        #println(debris_positions[i,t])
                        # vel = sqrt((G*M_earth)/norm(debris_list[i].r))
                        # unit_vel = debris_list[i].v/norm(debris_list[i].v) #asumes prior velocity vector is roughly same as current, HORRIBLE ASSUMPTION NEED TO FIX
                        #repeat update for velocity in debris list somehow, to ensure numerical integration is accurate
                        debris_list[i].v = delta_pos/t_prop

                        kepler_counter += 1
                        
                        
                    end
                end
            end
        end
        perc_kepler = kepler_counter/(kepler_counter+int_counter) * 100
        println("Percent of modeling using kepler")
        println(perc_kepler)
        return debris_positions

    end

    if int_method == 5 #Keplers method test platform
        @time begin
        laser_positions = Array{Vector{Float64},2}(undef,n_lasers,Int(t_span[2]))
            for i in 1:n_lasers
                laser_interval = LinRange(0,time,time*56)
                kepler_method!(t_span[2],laser_list[i],laser_positions,i)
            end
        end
        return laser_positions
    end



end


# Function to simulate orbital encounter
# t_span, LADR1, debris1, laser_output_power
# Extract the positions over time from the integration results



function plot_yoshida(results_yoshi,n_lasers,n_debris)
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
end

function plot_mdpt()

    # Extract positions from the results
    laser_positions, debris_positions = extract_positions(results_yoshi, n_lasers, n_debris)

    # Generate Earth sphere for visualization
    n = 100  # Resolution of the sphere
    x_sphere, y_sphere, z_sphere = generate_sphere([0, 0, 0], R_earth, n)

    # Plotting the trajectories with Plotly
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

    for i in 1:n_lasers
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
end


function midpt_orbital_encounter(t,y,laser_power)

    #parse inputs to set as laser and debris IC for the return
    n_obj = 1
    for i in 0:(n_obj-1)
        if i < n_lasers
            # Lasers
            #println("Processing laser index: $i")
            laser_list[i+1].r[1] = y[i*6 + 1]
            laser_list[i+1].r[2] = y[i*6 + 2]
            laser_list[i+1].r[3] = y[i*6 + 3]
            laser_list[i+1].v[1] = y[i*6 + 4]
            laser_list[i+1].v[2] = y[i*6 + 5]
            laser_list[i+1].v[3] = y[i*6 + 6]
        else
            # Debris
            debris_index = i - n_lasers
            #println("Processing debris index: $debris_index")
            debris_list[debris_index+1].r[1] = y[i*6 + 1]
            debris_list[debris_index+1].r[2] = y[i*6 + 2]
            debris_list[debris_index+1].r[3] = y[i*6 + 3]
            debris_list[debris_index+1].v[1] = y[i*6 + 4]
            debris_list[debris_index+1].v[2] = y[i*6 + 5]
            debris_list[debris_index+1].v[3] = y[i*6 + 6]
        end
    end

    #Calculate gravitational acceleration for Laser platforms and debris
    for i in 1:n_lasers
        x = laser_list[i].r[1]
        y = laser_list[i].r[2]
        z = laser_list[i].r[3]
        #println(sqrt(x^2 + y^2 + z^2))
        ax,ay,az = calculate_gravitational_acceleration(x,y,z)
        #println(ax,ay,az)
        laser_list[i].accel = [ax,ay,az]
    
    end
    for i in 1:n_debris
        x = debris_list[i].r[1]
        y = debris_list[i].r[2]
        z = debris_list[i].r[3]
        #println(sqrt(x^2 + y^2 + z^2))
        ax,ay,az = calculate_gravitational_acceleration(x,y,z)
        debris_list[i].accel = [ax,ay,az]
    end



    for i in 1:n_lasers
        #sorting debris based on distance
        #sorted_debris= sort_debris_by_distance(laser_list[i],debris_list)
        # Calculate the relative distance between L'ADROIT and the debris
        for j in 1:n_debris
            relative_distance = norm(laser_list[i].r - debris_list[j].r)

            if relative_distance < laser_range  # Assuming 250 km as effective laser range
                # Use small debris pulse energy
                delta_a = calculate_accel_change(cm, debris_list[j].mass, laser_power)

                # Apply delta_v in the opposite direction of the debris's velocity vector
                velocity_vector = debris_list[j].v
                velocity_magnitude = norm(velocity_vector)
        
                if velocity_magnitude > 0   # To avoid division by zero
                    # Normalize velocity vector and apply delta_v in the opposite direction
                    normalized_velocity_vector = velocity_vector / velocity_magnitude

                    acceleration_change = -normalized_velocity_vector * delta_a * laser_power
                    debris_list[j].accel[1] += acceleration_change[1]
                    debris_list[j].accel[2] += acceleration_change[2]
                    debris_list[j].accel[3] += acceleration_change[3]
                end
            else
                break
            end
        end
    end

    ret = zeros(n_obj * 6)  # This will store the derivative state

    for i in 0:(n_obj-1)
        if i < n_lasers
            # Lasers
            #println("Processing laser index: $i")
            ret[i*6 + 1] = laser_list[i+1].v[1]
            ret[i*6 + 2] = laser_list[i+1].v[2]
            ret[i*6 + 3] = laser_list[i+1].v[3]
            ret[i*6 + 4] = laser_list[i+1].accel[1]
            ret[i*6 + 5] = laser_list[i+1].accel[2]
            ret[i*6 + 6] = laser_list[i+1].accel[3]
        else
            # Debris
            debris_index = i - n_lasers + 1
            #println("Processing debris index: $debris_index")
            ret[i*6 + 1] = debris_list[debris_index].v[1]
            ret[i*6 + 2] = debris_list[debris_index].v[2]
            ret[i*6 + 3] = debris_list[debris_index].v[3]
            ret[i*6 + 4] = debris_list[debris_index].accel[1]
            ret[i*6 + 5] = debris_list[debris_index].accel[2]
            ret[i*6 + 6] = debris_list[debris_index].accel[3]
        end
    end
    
    #println("Acceleration vector (ret): ", ret)

    #println(ret)
    return ret
end

function yoshi_orbital_encounter(t,y, v,laser_power)
    n_lasers = 0

    #parse inputs to set as laser and debris IC for the return
    for i in 0:(n_obj-1)
        # if i < n_lasers
        #     # Lasers
        #     #println("Processing laser index: $i")
        #     laser_list[i+1].r[1] = y[i*3 + 1]
        #     laser_list[i+1].r[2] = y[i*3 + 2]
        #     laser_list[i+1].r[3] = y[i*3 + 3]
        # else
            # Debris
            debris_index = i - n_lasers
            #println("Processing debris index: $debris_index")
            debris_list[debris_index+1].r[1] = y[i*3 + 1]
            debris_list[debris_index+1].r[2] = y[i*3 + 2]
            debris_list[debris_index+1].r[3] = y[i*3 + 3]
        #end
    end

    #Calculate gravitational acceleration for Laser platforms and debris
    # for i in 1:n_lasers
    #     x = laser_list[i].r[1]
    #     y = laser_list[i].r[2]
    #     z = laser_list[i].r[3]
    #     ax,ay,az = calculate_gravitational_acceleration(x,y,z)
    #     laser_list[i].accel = [ax,ay,az]
    
    # end
    for i in 1:n_debris
        x = debris_list[i].r[1]
        y = debris_list[i].r[2]
        z = debris_list[i].r[3]
        ax,ay,az = calculate_gravitational_acceleration(x,y,z)
        debris_list[i].accel = [ax,ay,az]
    end



    for i in 1:n_lasers
        #sorting debris based on distance
        sorted_debris= sort_debris_by_distance(laser_list[i],debris_list)
        # Calculate the relative distance between L'ADROIT and the debris
        for j in 1:n_debris
            relative_distance = norm(laser_list[i].r - debris_list[j].r)

            if relative_distance < laser_range  # Assuming 250 km as effective laser range
                # Use small debris pulse energy
                delta_a = calculate_accel_change(cm, debris_list[j].mass, laser_power)

                # Apply delta_v in the opposite direction of the debris's velocity vector
                velocity_vector = debris_list[j].v
                velocity_magnitude = norm(velocity_vector)
        
                if velocity_magnitude > 0   # To avoid division by zero
                    # Normalize velocity vector and apply delta_v in the opposite direction
                    normalized_velocity_vector = velocity_vector / velocity_magnitude

                    acceleration_change = -normalized_velocity_vector * delta_a * 0
                    debris_list[j].r[1] += acceleration_change[1]
                    debris_list[j].r[2] += acceleration_change[2]
                    debris_list[j].r[3] += acceleration_change[3]
                end
            else
                break
            end
        end
    end

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

