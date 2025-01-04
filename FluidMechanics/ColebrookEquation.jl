function ColebrookEquation(eps, D, Re, L, g, V)
    ro=1000
    f = (-1.8*log((6.9/Re) + ((eps/D)/3.7)^1.11))^-2
    hf = f*(L/D)*((V^2)/(2*g))
    p=ro*g*hf
    return f, hf, p
end

# Pipe 1
f, hf, p = ColebrookEquation(1.5*10^-6, 0.032, 19311.8, 0.97, 9.81, 0.6047)

# Pipe 2
f, hf, p = ColebrookEquation(1.5*10^-6, 0.025, 24720.3, 0.99, 9.81, 0.99079)

# Pipe 3
f, hf, p = ColebrookEquation(6*10^-5, 0.025, 24720.3, 1.07, 9.81, 0.99079 )