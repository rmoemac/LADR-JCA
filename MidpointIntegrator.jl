function print_progress(current, total)
    percent = (current / total) * 100
    print("\rProgress: ", round(percent, digits=2), "%")
    flush(stdout)  # Ensure the output is flushed
end

function midpoint(tspan, y0, h,laser_power)
    #println("Midpoint Integrator Started")
    t0, tf = tspan
    t = t0
    y = y0
    y_filtered = remove_every_other_3(y)
    results = [(t, y_filtered)]
    n_objects = length(y) / 6
    y_next = y
    #println(y)


    #println(y)
    
    while t < tf
        if t + h > tf
            h = tf - t
        end

        #print_progress(t,tf)
        
        k1 = midpt_orbital_encounter(t, y,laser_power)
        y_mid = y .+ (h/2) .* k1
        k2 = midpt_orbital_encounter(t + h/2, y_mid,laser_power)
        
        y = y + (h.*k2)
        #println(y)
        #println(h*k2)
        #println(y_next - y)
        
        t += h

        #y_filtered = remove_every_other_3(y)
        #println(y_filtered)
        
        push!(results, (t, y))
    end
    
    return results
end