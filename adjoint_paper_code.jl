begin
    
    using Plots
    using Printf
    
    Δt = .001
    equilibrium = [0, 0]

    f(x, y) = [-y * (1 - x - y);
	       x * (1 - x - y)]

    g(x, y) = [x * (1 - y);
	       y * (-1 + x)]
    
    Q_func(x,y) = x^2 + y^2
    H_func(x,y) = x - log(x) + y - log(y)
    
    function f_initials(N, pert)
	x = zeros(N, 2)
	Q = zeros(N)
	x[1,1] = pert
	x[1,2] = pert
	Q[1] = Q_func(pert, pert)
	return x, Q
    end

    function g_initials(N, pert)
	x = zeros(N, 2)
	Q = zeros(N)
	x[1,1] = 1.0 + pert
	x[1,2] = 1.0 + pert
	Q[1] = H_func(x[1,1], x[1,2])
	return x, Q
    end
    
    function forward_euler(Δt, T, dx, pert, initials, conserve)
	t = 0:Δt:T
	N = length(t)
	x, Q = initials(N, pert)
	for tstep in 2:N
	    x[tstep, :] .= x[tstep-1, :] .+ Δt * dx(x[tstep-1, 1], x[tstep-1, 2])
	    Q[tstep] = conserve(x[tstep, 1], x[tstep, 2])
	end
	return x, t, Q
    end

    function semi_implicit_euler(Δt, T, dx, pert, initials, conserve)
	t = 0:Δt:T
	N = length(t)
	x, Q = initials(N, pert)
	# loop over timesteps
	for tstep in 2:N
	    
	    x[tstep, 1] = x[tstep-1, 1] + Δt * dx(x[tstep-1, 1], x[tstep-1, 2])[1]
	    x[tstep, 2] = x[tstep-1, 2] + Δt * dx(x[tstep, 1], x[tstep-1, 2])[2]
	    Q[tstep] = conserve(x[tstep, 1], x[tstep, 2])
	end
	return x, t, Q
    end

    function linear_forward_euler(Δt, T, x, dx, initials, conserve)
	t = 0:Δt:T
	N = length(t)
	x = zeros(N, 2)
	x[1, :] .= initials		 
	for tstep in 2:N
	    J = dx(x[tstep, 1], x[tstep, 2])
	    x[tstep, :] = x[tstep - 1, :] + Δt * (J * x[tstep - 1, :])
	end
	return x, t
    end

    function linear_semi_implicit_euler(Δt, T, x, dx, initials)
	t = 0:Δt:T
	N = length(t)
	δ = zeros(N, 2)
	δ[1, :] .= initials			
	for tstep in 2:N
	    J = dx(x[tstep, 1], x[tstep, 2])
	    δ[tstep, 1] = δ[tstep - 1, 1] + Δt * (J[1,:]' * δ[tstep - 1, :])
	    δ[tstep, 2] = δ[tstep - 1, 2] + Δt * (J[2, :]' * δ[tstep, :])
	end
	return δ, t
    end	

    function linear_reverse_semi_implicit_euler(Δt, T, x, dx, initials)
	t = 0:Δt:T
	N = length(t)
	λ = zeros(N, 2)
	λ[N, :] .= initials
        
	for tstep in N-1:-1:1
	    J = dx(x[tstep, 1], x[tstep, 2])'
	    λ[tstep, 1] = λ[tstep + 1, 1] + Δt * (J[1,:]' * λ[tstep + 1, :])
	    λ[tstep, 2] = λ[tstep + 1, 2] + Δt * (J[2, :]' * λ[tstep, :])
	end
	return λ, t
    end	

    function linear_reverse_forward_euler(Δt, T, x, dx, initials)
	t = 0:Δt:T
	N = length(t)
	λ = zeros(N, 2)
	λ[N, :] .= initials
	for tstep in N-1:-1:1
	    J = dx(x[tstep, 1], x[tstep, 2])'
	    λ[tstep, :] = λ[tstep + 1, :] + Δt * (J * λ[tstep + 1, :])
	end
	return λ, t
    end	
    
    function plot_solution(dx, integrator, initials, conserve, pert, T)
	x, t, Q = integrator(Δt, T, dx, pert, initials, conserve)
	plt1 = plot(t, x[:, 1], ylabel="x(t)", xlabel="t", legend=false)
	plt2 = plot(t, x[:, 2], ylabel="y(t)", xlabel="t", legend=false)
	plt3 = plot(t, Q, ylabel="Q(x,y)", xlabel="t")
	plot(plt1, plt2, plt3, layout=(3,1))
    end
    
    function plot_contours(dx, integrator, initials, conserve, n_sol, pert_step, T)
	steps = collect(pert_step:pert_step:n_sol*pert_step)
	plt = plot()
	for (sn, step) in enumerate(steps)
	    x, t, Q = integrator(Δt, T, dx, step, initials, conserve)
	    plt = plot!(x[:,1], x[:,2], legend=false, xlabel = "x(t)", ylabel="y(t)")
	end
	return plt
    end
    
    function plot_surf_solution(dx, x_low, x_high, y_low, y_high, T,
				integrator, initials, conserve, 
		                n_sol, sol_step, dcolor)
	surf_x = x_low: .01 : x_high
	surf_y = y_low: .01 : y_high
	plt = plot(surf_x, surf_y, conserve, st=:surface, camera=(30,30), fillalpha=0.5)
	
	for step in sol_step:sol_step:sol_step*n_sol					
	    x, t, Q = integrator(Δt, T, dx, step, initials, conserve)
	    plt = scatter!(x[:,1], x[:,2], Q_func.(x[:,1], x[:2]),
                           color=dcolor, markersize=.1, legend=false,
                           markerstrokewidth=.01)
	end
	return plt
    end

    
    
    plot_solution(f, semi_implicit_euler, f_initials, Q_func, .01, 5000)
    plot_solution(f, forward_euler, f_initials, Q_func, .01, 5000)
    plot_contours(f, forward_euler, f_initials, Q_func, 10, .01, 100)
    plot_contours(f, semi_implicit_euler, f_initials, Q_func, 10, .01, 100)
    plot_surf_solution(f, -.2, .2, -.2, .2, 100,
		       semi_implicit_euler, f_initials, Q_func, 
		       10, .01, :green)
    plot_surf_solution(f, -.2, .2, -.2, .2, 100,
		       forward_euler, f_initials, Q_func, 
		       10, .01, :red)


    
    Jg(x,y) = [1-y (-x) ;
	       (y) x-1]
    Gequal = [1, 1]	

    w = [.01, 0.0]
    η = [0.01, 0.0]
    # linear forward equations
    x, t, Q = forward_euler(Δt, 2, g, .1, g_initials, H_func)
    # sensitivity
    δ, t = linear_semi_implicit_euler(Δt, 2, x, Jg, η)
    # symplactic adjoint 
    λ₁, t = linear_reverse_semi_implicit_euler(Δt, 2, x, Jg, w)
    # non-symplactic adjoint
    λ₂, t = linear_reverse_forward_euler(Δt, 2, x, Jg, w)
    
    @show w' * δ[end, :]
    @show λ₁[1,:]' * η
    @show λ₂[1,:]' * η
    @show abs(λ₁[1,:]' * η - w' * δ[end, :])/(w' * δ[end, :])
    @show abs(λ₂[1,:]' * η - w' * δ[end, :])/(w' * δ[end, :])
    
    
    plot(x[:, 1], x[:, 2], label = "solution")
    plot!(x[:, 1] .+ δ[:, 1], x[:, 2] .+ δ[:, 2], label= "x perturbed solution via sensativity")
    plot!(x[:, 1] .+ λ₁[:, 1], x[:, 2] .+ λ₁[:,2], label= "symplactic-adjoint")
    plot!(x[:, 1] .+ λ₂[:, 1], x[:, 2] .+ λ₂[:,2], label= "non-symplactic-adjoint")
        
    
end
