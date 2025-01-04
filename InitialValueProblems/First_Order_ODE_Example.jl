using Plots 

function EulerMethod(f, x0, y0, h, n)
X = Float64[]
Y = Float64[]
push!(X, x0)
push!(Y, y0)
for i in 1:n
    x = X[end]+h
    push!(X, x)
    y = Y[end]+h*f(X[end])
    push!(Y, y) 
end
return X,Y
end
f(x) = 2x # this is the derivative or tangent of the function
X, Y = EulerMethod(f, 1, 1, 0.1, 5)

plot(X, Y, label = "Solution for y=x^2", xlabel = "x values", ylabel = "y values")

# find exact values by plugging x in real function 

function ExactSolution(X, n)
    exactY = Float64[]
    for i in 1:n
        yValue = (X[i])^2
        push!(exactY, yValue)
    end
    return exactY
end

exactY = ExactSolution(X, length(X))

# Plot numerical solution and exact solution
plot(X, Y, label = "Solution for y = x^2", xlabel = "x values", ylabel = "y values")
plot!(X, exactY, label = "Exact solution", lw = 1, linestyle = :solid, color = :black)

# Save the combined plot
savefig("euler_method_solution.png")

