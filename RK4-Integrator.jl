function print_progress(current, total)
    percent = (current / total) * 100
    print("\rProgress: ", round(percent, digits=2), "%")
    flush(stdout)  # Ensure the output is flushed
end

function rk4(tspan, y0, h)
    print("RK4 Started")
    #println(y0)
    t0, tf = tspan
    t = t0
    y = y0
    results = [(t, y)]
    
    while t < tf
        if t + h > tf
            h = tf - t
        end
        print_progress(t,tf)
        
        k1 = midpt_orbital_encounter(t, y)
        k2 = midpt_orbital_encounter(t + h/2, y .+ (h/2) .* k1)
        k3 = midpt_orbital_encounter(t + h/2, y .+ (h/2) .* k2)
        k4 = midpt_orbital_encounter(t + h, y .+ h .* k3)
        
        y += (h/6) .* (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4)
        #println((sqrt(y[1]^2 + y[2]^2 + y[3]^2) - R_earth)/1000)
        t += h
        
        push!(results, (t, y))

    end
    
    return results
end