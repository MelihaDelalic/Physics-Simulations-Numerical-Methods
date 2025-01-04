using Plots

# Euler Method implementation
function EulerMethod(f, x0, y0, h, n) # where x0 and y0 represent a point that is our initial condition
    X = Float64[] # Array to store x values
    Y = Float64[] # Array to store y values
    push!(X, x0)
    push!(Y, y0)
    for i in 1:n
        x = X[end] + h
        push!(X, x)
        y = Y[end] + h * f(X[end], Y[end])
        push!(Y, y)
    end
    return X, Y
end

# Define the function f(y, t) = y - t^2 + 1
f(x, y) = y - x^2 + 1

# Call EulerMethod
X, Y = EulerMethod(f, 0, 0.5, 0.05, 50)

# Exact solution function
exact(t) = (t + 1)^2 - 0.5 * exp(t)

# Calculate exact solution for the same X values
exa = exact.(X)

# Plot Euler's method and exact solution
plot(X, Y, label = "Euler's Method", xlabel = "x", ylabel = "y", title = "Numerical Solution", lw = 2)
plot!(X, exa, label = "Exact Solution", lw = 2, linestyle = :dash, color = :red)

# Calculate absolute error
absoluteError = abs.(Y .- exa)
println("Absolute error: ", absoluteError)

# Plot absolute error
plot(X, absoluteError, label = "Absolute Error", xlabel = "x values", ylabel = "Absolute error difference", lw = 2)
