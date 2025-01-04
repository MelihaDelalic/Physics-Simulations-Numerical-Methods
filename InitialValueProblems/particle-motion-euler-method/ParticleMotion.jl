using Plots

function ParticleMotion(v0, t0, x0, m, n, tau)
X = Float64[]
V = Float64[]
T = Float64[]
push!(X, x0)
push!(V, v0)
push!(T, t0)
    for i in 1:n
    x = X[end]
    v = V[end]
    t = T[end]
    force = f(x)

    # Update the new values for position, velocity and time
    x_new = x + tau*v
    v_new = v + (tau/m)*force 
    t = t + tau

    # Add the new values to arrays
    push!(X, x_new)
    push!(V, v_new)
    push!(T, t)
    end
    return X, V, T
end

# Force function
f(x) = -1.0 *x  # spring force with k = 1
tau = (2*pi)/100000


X, V, T = ParticleMotion(1.0, 0.0, 0.0, 1.0, 100000.0, tau)

# Plot velocity vs time
plot(T, V, label = "Velocity of the particle", xlabel = "Time", ylabel = "Velocity", 
     linestyle = :dash, color = :red, lw = 3)

# Plot position vs time on the same graph
plot!(T, X, label = "Position of the particle", xlabel = "Time", ylabel = "Position", 
      lw = 3, linestyle = :solid, color = :blue)

# Save the combined plot
savefig("particle_motion_euler.png")

