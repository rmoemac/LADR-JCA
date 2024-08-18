using LinearAlgebra

# Define constants
const G = 6.67430e-11  # Gravitational constant in m^3 kg^-1 s^-2
const M_earth = 5.972e24  # Mass of the Earth in kg

# Convert degrees to radians
function deg2rad(deg)
    return deg * pi / 180.0
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

# Function for calculating gravitational acceleration
function calculate_gravitational_acceleration(x, y, z)
    r = sqrt(x^2 + y^2 + z^2)
    ax = -G * M_earth * x / r^3
    ay = -G * M_earth * y / r^3
    az = -G * M_earth * z / r^3
    return ax, ay, az
end

# Function for calculating orbital energy
function calculate_orbital_energy(v, r, mu)
    kinetic_energy = 0.5 * norm(v)^2
    potential_energy = r == 0 ? 0 : -mu / norm(r)
    return kinetic_energy + potential_energy
end

# Function for calculate_momentum_change and delta_accel (placeholders, need proper implementation)
function calculate_momentum_change(pulse_energy, drag_force, Cm, target_mass, laser_output_power, pulse_width)
    laser_output_power = laser_output_power * 1e-6
    newtons = Cm * laser_output_power
    dt = 1e-6
    delta_momentum = newtons * dt
    
    delta_velocity = delta_momentum / target_mass
    return delta_velocity
end

function calculate_accel_change(pulse_energy, Cm, target_mass, laser_output_power, pulse_width)
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