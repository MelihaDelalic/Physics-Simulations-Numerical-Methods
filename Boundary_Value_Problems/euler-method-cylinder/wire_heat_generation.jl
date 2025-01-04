########### Solving second order differential equation using Euler's Method #################
using Plots 

## Properties of the resistance wire are:
r_0 = 0.006              # m 
z_0 = 0.0                # First derivative of temperature function at the center is 0 
e_gen = 16.4 * 10^6      # W/(m^3)  # Heat generation in the rod 
k = 15.2                 # Thermal Conductivity in W/(m*K)
h = 3200.0               # Convective heat transfer coefficient 
Δr = 0.00005             # Step size 
T_inf = 373 # K          # Temperature of the water given initially 100°C

function TemperatureDistribution(T_0:: Int64, r_0::Float64, z_o::Float64, e_gen::Float64, k::Float64, h::Float64, Δr::Float64,  T_inf::Int64, n::Int64)
    
    # Initialize arrays for storing values 
    T = Float64[]      # Temperature values 
    Z = Float64[]      # First derivative of the temperature function
    R = Float64[]      # Radius of the rod
    
    # Push the initial conditions in the arrays, starting from radius r = 0 outwards
    push!(T, T_0)  # For the initial temperature we start by guessing the value
    push!(Z, z_0)  # Here we apply our first boundary condition - thermal symmetry at the center
    push!(R, 0.0)  # The calculations start from r = 0
    # Create the local variables for the arrays 
    t = T[end]
    z = Z[end]
    r = R[end]
    # Apply Euler's method to obtain new values of each variable 
    for i in 1:n
        
        t_new = t + Δr*z 
        
        if r > 0
            z_new = z + Δr * (-e_gen/k - (z/r))
        else
            z_new = z + Δr * (-e_gen/k)  # Handle r = 0 case
        end
        
        r_new = r + Δr
        println("Step: $i, t: $t, z: $z, r: $r")  # Debug output
        # Now push the new values into their array 
        push!(T, t_new)
        push!(Z, z_new)
        push!(R, r_new)
        # Update the local variables for next iteration:
        t = t_new
        z = z_new
        r = r_new
    end
    # Now we apply the conduction-convection boundary condition 
    t_surface = T[end]
    z_surface = Z[end]
    r_surface = r_0
    # The first derivative of the temperature at the surface of the rod
    z_r0 = (-h/k)*(t_surface-T_inf)
    t_surface_new = t_surface - Δr*z_r0  # We account for the heat loss by convection at the surface
    push!(T, t_surface_new)
    push!(R, r_surface)
    return T, R
end

temp, radius = TemperatureDistribution(398, r_0, z_0, e_gen, k, h, Δr, T_inf, 120)

# Extract the surface temperature
surface_temperature = temp[end]

# Visualize the results 
plot(radius, temp, xlabel = "Radius (m)", ylabel = "Temperature (K)", label = "Temperature Distribution", lw=2, color=:red, linestyle=:solid)
annotate!(r_0, temp[end], text("Surface Temp", 10, :blue))
temp_Celsius = temp.- 273
plot(radius, temp_Celsius, xlabel = "Radius (m)", ylabel = "Temperature (°C)", label = "Temperature Distribution", lw=2, color=:red, linestyle=:solid)
savefig("Temperature Distribution in a Cylindrical Rod - Euler Method.png")