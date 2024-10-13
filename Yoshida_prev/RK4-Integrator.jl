

function rk4(tspan, y0, h)
    print("RK4 Started")
    t0, tf = tspan
    t = t0
    y = y0
    results = [(t, y)]
    
    while t < tf
        if t + h > tf
            h = tf - t
        end
        
        k1 = simulate_orbital_encounter(t, y)
        k2 = simulate_orbital_encounter(t + h/2, y .+ (h/2) .* k1)
        k3 = simulate_orbital_encounter(t + h/2, y .+ (h/2) .* k2)
        k4 = simulate_orbital_encounter(t + h, y .+ h .* k3)
        
        y += (h/6) .* (k1 .+ 2 .* k2 .+ 2 .* k3 .+ k4)
        t += h
        
        push!(results, (t, y))

    end
    
    return results
end