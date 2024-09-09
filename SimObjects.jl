mutable struct LADR
    a::Float64
    e::Float64
    i_deg::Float64
    omega_deg::Float64
    raan_deg::Float64
    nu_deg::Float64
    power::Float64
    r::Vector{Float64}
    v::Vector{Float64}
    accel::Vector{Float64}

    function LADR(a::Float64,e::Float64,i_deg::Float64,omega_deg::Float64,raan_deg::Float64,nu_deg::Float64,power::Float64)
        r, v = orbital_elements_to_eci(a, e, i_deg, omega_deg, raan_deg, nu_deg)
        new(a,e,i_deg,omega_deg,raan_deg,nu_deg,power,r,v)
    end
end

mutable struct Debris
    a::Float64
    e::Float64
    i_deg::Float64
    omega_deg::Float64
    raan_deg::Float64
    nu_deg::Float64
    mass::Float64
    r::Vector{Float64}
    v::Vector{Float64}
    accel::Vector{Float64}

    function Debris(a::Float64,e::Float64,i_deg::Float64,omega_deg::Float64,raan_deg::Float64,nu_deg::Float64,mass::Float64)
        r, v = orbital_elements_to_eci(a, e, i_deg, omega_deg, raan_deg, nu_deg)
        new(a,e,i_deg,omega_deg,raan_deg,nu_deg,mass,r,v)
    end
end
