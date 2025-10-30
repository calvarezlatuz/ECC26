# Kinetics
function μ(s, i) 
    i = trunc(Int, i)
    return (a[i]*s)/(b[i]+s)
end
# Perturbation function 
function  M(θ, ε, x, s, i)
    i = trunc(Int, i)
    return ε*sum(matrix_Z[i,:].*x)+(θ/sum(x))*sum(matrix_H[i,:].*x)*x[i]
end 
# Dynamics of s
function dS(x, s, u)
    mu = (a*s)./(b.+s)
    return (u*(s_in-s))-sum(mu.*(x./Y))
end
# Simpson index S=1-D
function D(x)
    N = sum(x)^2
    return sum((x.^2))/N
end
# Production 
function P(u, x)
    return u*sum(x)
end;
# Dynamics 
function f(y, p, t)
    x, s = y[1:n], y[n+1]
    u, ε, θ = p
    dy = zeros(n+1)
    dy[1:n] = [((μ(s,i) - u) * x[i]) + M(θ,ε,x,s,i) for i in 1:n]
    dy[n+1] = dS(x, s, u)
    return dy
end
# Simulate trajectories
function solvePCS(f, y0, tf, u, ε, θ)
    p = [u, ε, θ]
    tspan = (0.0, tf)
    prob = ODEProblem(f, y0, tspan, p)
    return solve(prob)
end;