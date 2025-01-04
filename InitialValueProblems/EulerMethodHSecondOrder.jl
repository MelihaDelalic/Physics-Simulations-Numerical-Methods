using Plots
# Define the function for the system of first-order ODEs
function f(y, z)
    return z, -y  # Returns the derivatives: y' = z and z' = -y
end
function EulerMethodSecondOrder(f::Function, y0::Float64, z0::Float64, h::Float64, n::Int64)
    
    # Initialize the arrays for storing values of y and z 
    Y = Float64[]   # This Y array represents our solution to the second order differential equation
    Z = Float64[]   # The Z array represents the first derivative of this equation
    X = Float64[]
    # Push the existing initial conditions in the arrays 
    push!(Y, y0)
    push!(Z, z0)
    push!(X, 0.0)
    
        for i in 1:n
         # Also initialize the existing values locally 
        y = Y[end]
        z = Z[end]
        x = X[end]
            # Calculate the next values of these arrays 
             y_new = y + h*z 
             z_new = z + h*(-y)
             x_new = x + h 
             # Push new values into the arrays
            push!(Y, y_new)
            push!(Z, z_new)
            push!(X, x_new)
        end
        return Y, Z, X
end

Y, Z, X = EulerMethodSecondOrder(f, 1.0, 0.0, 0.1, 16)
plot(X, Y, xlabel = "X values", ylabel = "Y values", label = "Solution to y''(x)=-y(x)", linestyle=:solid, lw=2, color=:blue)
plot!(X, Z, ylabel="Y values", label="First Derivative (Z)", linestyle=:dash, color=:red, lw=2)
savefig("Second Order ODE Solution.png")