# Define a function to solve Kepler's equation using the Newton-Raphson method
using LinearAlgebra

function kepler_E(M, e, tol=1e-8, max_iter=40)
    E = M
    for _ in 1:max_iter
        E_new = E - (E - e * sin(E) - M) / (1 - e * cos(E))
        if abs(E_new - E) < tol
            return E_new
        end
        E = E_new
    end
    # RuntimeError(f"Kepler's equation did not converge after {max_iter} iterations")
end
using LinearAlgebra

function true_anomaly(E, e)
    return 2 * atan(sqrt(1 + e) * sin(E / 2), sqrt(1 - e) * cos(E / 2))
end
function orbital_position(a, e, E)
    r = a * (1 - e * cos(E))
    theta = true_anomaly(E, e)
    return r, theta
end
function kepler_method_event(t,t_prop,object,output_positions,output_index)
    a = object.a
    e = object.e
    i = deg2rad(object.i_deg)
    raan = deg2rad(object.raan_deg)
    argp = deg2rad(object.omega_deg)
    temp_positions = []

    T = 2 * pi * sqrt(a^3 / mu)
    # for t in 1:t_prop
    #     M = sqrt(mu / a^3) * t
    #     M_t = M + 2 * pi * t / T
    #     E = kepler_E(M_t, e)
    #     r, theta = orbital_position(a, e, E)
    #     x = r * (cos(raan) * cos(argp + theta) - sin(raan) * sin(argp + theta) * cos(i))
    #     y = r * (sin(raan) * cos(argp + theta) + cos(raan) * sin(argp + theta) * cos(i))
    #     z = r * sin(i) * sin(argp + theta)
        
    #     temp_positions = [x,y,z]
    # end
    t = t_prop
    M = sqrt(mu / a^3) * t
    M_t = M + 2 * pi * t / T
    E = kepler_E(M_t, e)
    r, theta = orbital_position(a, e, E)
    x = r * (cos(raan) * cos(argp + theta) - sin(raan) * sin(argp + theta) * cos(i))
    y = r * (sin(raan) * cos(argp + theta) + cos(raan) * sin(argp + theta) * cos(i))
    z = r * sin(i) * sin(argp + theta)
        
    temp_positions = [x,y,z]
    #output_positions[output_index,(t+t_prop)] = temp_positions
    #println(temp_positions)
    return temp_positions

end

function kepler_method!(t_length,object,output_positions,output_index)
    a = object.a
    e = object.e
    i = deg2rad(object.i_deg)
    raan = deg2rad(object.raan_deg)
    argp = deg2rad(object.omega_deg)

    T = 2 * pi * sqrt(a^3 / mu)
    for t in 1:t_length
        M = sqrt(mu / a^3) * t
        M_t = M + 2 * pi * t / T
        E = kepler_E(M_t, e)
        r, theta = orbital_position(a, e, E)
        x = r * (cos(raan) * cos(argp + theta) - sin(raan) * sin(argp + theta) * cos(i))
        y = r * (sin(raan) * cos(argp + theta) + cos(raan) * sin(argp + theta) * cos(i))
        z = r * sin(i) * sin(argp + theta)
        
        output_positions[output_index,t] = [x,y,z]
    end


end

function solve_kepler(M, e; tol=1e-8, max_iter=1000)
    # Initial guess for E
    E = M
    
    for i in 1:max_iter
        f = E - e*sin(E) - M
        fp = 1 - e*cos(E)
        dE = -f / fp
        E += dE
        
        # Check for convergence
        if abs(dE) < tol
            return E
        end
    end
    
    error("Kepler's equation did not converge")
end

# Function to convert eccentric anomaly E to true anomaly ν
function eccentric_to_true_anomaly(E, e)
    ν = 2 * atan(sqrt((1 + e) / (1 - e)) * tan(E / 2))
    return ν
end