#YOSHIDA 4TH ORDER INTEGRATOR ATTEMPT 

# c1 = 0.6756
# c2 = -0.1756
# c3 = -0.1756
# c4 = 0.6756

w0 = -(2.0^(1/3))/(2.0-(2.0^(1/3)))
w1 = 1.0/(2.0-(2.0^(1/3)))
c1 = w1/2.0
c4 = c1
c2 = (w0+w1)/2.0
c3 = c2
d1 = w1
d3 = d1
d2 = w0

function yoshida(tspan,y0,h) #h is dt
    print("yoshida started")
    t0, tf = tspan
    t = t0
    y = copy(y0)
    results = [(t, y)]

    n_objects = length(y0)/6
    n_objects = Int64(n_objects)

    xic = Float64[]
    vic = Float64[]

    output = zeros(Float64,6*n_objects)

    #Loading initial conditions into xi and vi for integration
    for i in 0:n_objects - 1
        x = y0[i*6 + 1]
        y = y0[i*6 + 2]
        z = y0[i*6 + 3]
        vx = y0[i*6 + 4]
        vy = y0[i*6 + 5]
        vz = y0[i*6 + 6]

        push!(xic,x)
        push!(xic,y)
        push!(xic,z)
        push!(vic,vx)
        push!(vic,vy)
        push!(vic,vz)
    end
    xi = xic
    vi = vic

    
    print("integration started")

    while t < tf

        if t + h > tf
            h = tf - t
        end

        xi1 = xi .+ c1 .* vi .* h
        axi1 = yoshi_orbital_encounter(t, xi1, vi)
        vi1 = vi .+ (d1 .* axi1 .* h)

        xi2 = xi1 .+ c2 .* vi1 .* h
        axi2 = yoshi_orbital_encounter(t, xi2, vi1)
        vi2 = vi1 .+ (d2 .* axi2 .* h)

        xi3 = xi2 .+ c3 .* vi2 .* h
        axi3 = yoshi_orbital_encounter(t, xi3, vi2)
        vi3 = vi2 .+ (d3 .* axi3 .* h)

        # computing state at next timestep
        xi = xi3 .+ (c4 .* vi3 .* h)
        vi = vi3

        #repacking into acceptable output form
        output = zeros(Float64, 6 * n_objects)
        for i in 0:n_objects -1
            output[i*6 + 1] = xi[i*3 + 1]
            output[i*6 + 2] = xi[i*3 + 2]
            output[i*6 + 3] = xi[i*3 + 3]
            output[i*6 + 4] = vi[i*3 + 1]
            output[i*6 + 5] = vi[i*3 + 2]
            output[i*6 + 6] = vi[i*3 + 3]
        end

        t += h
        
        push!(results, (t, output))
    end
    
    return results
end