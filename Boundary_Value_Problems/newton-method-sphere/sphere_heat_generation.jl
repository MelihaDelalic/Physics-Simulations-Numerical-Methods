using LinearAlgebra, Printf, Plots


e_gen = 5 * 10^7                        # W/m^3
k = 15                                  # W/m°C
M = e_gen/k                             # Calculating constant beforehand
Ts = 110                                # Surface Temperature (°C)
r1 = 1e-4                               # Initial radius (m)
R = 0.04                                # Sphere outer radius (m)
n = 1000                                # Number of intervals 
h = R/(n+1)                             # Dividing radius into steps
r = collect(r1:h:R)                     # Creating vector representing the radius
T = collect(LinRange(1050, Ts, n-1))    # Creating temperature vector of size n-1, with initial temperature guess being 1050 and final temperature Ts


# Define the system of equations
function Equations(r1, r, h, T, M, n, Ts)
 
    f = zeros(n-1)

     f[1] = r1 * (T[2]-T[1]) + h *(T[2]-T[1]) + M*r1*(h^2)

     for i in 2:n-2 
       f[i] = r[i] * (T[i-1] - 2 * T[i] + T[i+1]) + h * (T[i+1] - T[i-1]) + M * r[i] * (h^2)
     end

     f[n-1] = r[n-1] * (T[n-2] - 2 * T[n-1] + Ts) + h * (Ts - T[n-2]) + M * r[n-1] * (h^2)

     return f

end

# Create the Jacobian Function
function Jacobian(r, h, n)

    # Create a 2-dimensional array of n-1 rows and n-1 columns
    jacobian = zeros(n-1, n-1) 

    # Populate the first row of the jacobian
    jacobian[1,1] = - r[1] - h
    jacobian[1,2] = r[1] + h 

    # Populate the interior of the matrix 
    for i in 2:n-2
        jacobian[i, i-1] = r[i] - h 
        jacobian[i, i] = -2 * r[i] 
        jacobian[i, i+1] = r[i] + h
    end

    # Finally, populate the last row 
    jacobian[n-1, n-2] = r[n-1] - h
    jacobian[n-1, n-1] = -2 * r[n-1]

    return jacobian

end

# Apply the Newton method
function Newton_Raphson_Method(T::Vector, r::Vector, h, M, n, max_iter = 100, tol = 1e-6) 

    for i in 1:max_iter

        # Evaluate the function values at current guess 
        F = Equations(r1, r, h, T, M, n, Ts)

        # Compute the jacobian for the curret values 
        J = Jacobian(r, h, n)

        # Compute the temperature difference at current step
        ΔT = -J \ F
        
        # Update temperature 
        T = T + ΔT 
       
        # Evaluate the the system of equations
        Fval = Equations(r1, r, h, T, M, n, Ts)

        # Check for convergence (using the norm of F)
        if norm(Fval) < tol
            return T 
        end

    end

end

# Testing 
max_iter = 1000
tol = 1e-6

solution = Newton_Raphson_Method(T, r, h, M, n, max_iter, tol)
println(solution)
# length(r)

T_center = solution[1]
@printf("The center temperature is: %.2f degrees Celsius.\n", T_center)

# Plot the results 

# Plots are done with T = 1050 °C for initial guess 
plot(r, solution, xlabel = "Radius (m)", ylabel = "Temperature (°C)", title = "Temperature distribution across the radius.", legend=true, label="Temperature function\nn=1000")

# Create arrays for the center and end points
center_r = [r[1]]
center_T = [solution[1]]
end_r = [r[end]]
end_T = [solution[end]]

# Add red dots at the center (T(0)) and at the end (T(R)) with labels for the legend
scatter!(center_r, center_T, color=:red, label="Center (T(0))", marker=:circle, markersize=4)
scatter!(end_r, end_T, color=:red, label="End (T(R))", marker=:circle, markersize=4)
savefig("sphere_n=1000.png")
