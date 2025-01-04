using Plots

function HalleysCometTaylor(x0, y0, vx0, vy0, t0, h, n, m)
    # Constants (normalized for this case)
    k = 39.478428  # Gravitational constant (scaled) # it accounts for both 
    
    # Initialize the arrays for X, Y, velocities, etc.
    X = zeros(Float64, n)
    Y = zeros(Float64, n)
    Vx = zeros(Float64, n)
    Vy = zeros(Float64, n)
    T = zeros(Float64, n)
    R = zeros(Float64, n)
    Gx = zeros(Float64, n)
    Gy = zeros(Float64, n)
    
    # Initial conditions
    X[1] = x0
    Y[1] = y0
    Vx[1] = vx0
    Vy[1] = vy0
    R[1] = sqrt(X[1]^2 + Y[1]^2)
    Gx[1]   = -k * X[1] / (R[1]^3)
    Gy[1] = -k * Y[1] / (R[1]^3)
    T[1] = t0
    
    h2 = h^2 / 2  # Pre-compute h^2 / 2 for efficiency
    
    # Time-stepping loop using 2nd-order Taylor expansion
    for i in 1:(n-1)

        # Update time
        T[i + 1] = T[i] + h
        
        # Update position
        X[i + 1] = X[i] + Vx[i] * h + Gx[i] * h2
        Y[i + 1] = Y[i] + Vy[i] * h + Gy[i] * h2
        
        # Compute new distance and accelerations
        R2 = X[i + 1]^2 + Y[i + 1]^2 # square of the position vector
        R[i + 1] = sqrt(R2)   # the actual position vector
        R3 = R2 * R[i + 1]    # the cube of the position vector
        
        Gx[i + 1] = -k * X[i + 1] / R3
        Gy[i + 1] = -k * Y[i + 1] / R3
        
        # Update velocity with 2nd-order correction
        Vx[i + 1] = Vx[i] + h * Gx[i] + h2 * (3 * X[i] * (X[i] * Vx[i] + Y[i] * Vy[i]) / (R[i]^2) - Vx[i]) / (R[i]^3)
        Vy[i + 1] = Vy[i] + h * Gy[i] + h2 * (3 * Y[i] * (Y[i] * Vy[i] + X[i] * Vx[i]) / (R[i]^2) - Vy[i]) / (R[i]^3)
    end
    
    # Plot the trajectory of the comet (Y position vs X position)
    plot(X, Y, xlabel = "Position in X (m)", ylabel = "Position in Y (m)", 
         label = "Trajectory of Halley's Comet", linestyle = :solid, lw = 2, color = :blue)
    savefig("comet_position_taylor.png")
    
    # Plot velocity in X direction over time
    plot(T[1:m:end], Vx[1:m:end], xlabel = "Time (s)", ylabel = "Velocity in X (m/s)", 
         label = "X Velocity over time", color = :green, lw = 2)
    savefig("comet_velocity_x_taylor.png")
    
    # Plot velocity in Y direction over time
    plot(T[1:m:end], Vy[1:m:end], xlabel = "Time (s)", ylabel = "Velocity in Y (m/s)", 
         label = "Y Velocity over time", color = :red, lw = 2)
    savefig("comet_velocity_y_taylor.png")
    
    return X, Y, Vx, Vy, T
end

# Call the function to simulate Halley's comet with parameters
X, Y, Vx, Vy, T = HalleysCometTaylor(1.966843, 0.0, 0.0, 0.815795, 0.0, 2.0/199999, 200000, 2000)
