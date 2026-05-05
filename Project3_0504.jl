### A Pluto.jl notebook ###
# v0.20.24

using Markdown
using InteractiveUtils

# ╔═╡ d0c26027-876d-4426-8250-e3def024046d
begin
	import Pkg   # ← YOU ARE MISSING THIS LINE

	Pkg.activate(temp=true)

	Pkg.add("DifferentialEquations")
	Pkg.add("Plots")

	using DifferentialEquations
	using Plots
end

# ╔═╡ 0d0636e0-3225-4e3d-ba54-180e3f726f14
begin
	# Parameters
	m1 = 0.1      # kg
	m2 = 0.3      # kg
	k  = 10.0     # N/m
	g  = 9.81     # m/s^2
	L  = 0.4      # m

	J2 = (1/12)*m2*L^2  # slender bar inertia about its center
end

# ╔═╡ 984ec51b-75e6-4c0e-b8f5-78610983da1b
begin
	# Initial conditions
	x0 = 0.05              # block displacement, m
	xdot0 = 0.0            # block velocity, m/s
	θ0 = deg2rad(60.0)     # initial bar angle
	θdot0 = 0.0            # angular velocity, rad/s

	u0 = [x0, xdot0, θ0, θdot0]
	tspan = (0.0, 10.0)
end

# ╔═╡ c6ec4d9c-cb63-4c25-a5ac-5ae3af3ceb88
function eom!(du, u, p, t)
	m1, m2, k, g, L, J2 = p

	x = u[1]
	xdot = u[2]
	θ = u[3]
	θdot = u[4]

	A = [
		m1 + m2                 -(m2*L/2)*sin(θ)
		-(m2*L/2)*sin(θ)        J2 + m2*L^2/4
	]

	b = [
		(m2*L/2)*cos(θ)*θdot^2 - k*x
		-(m2*g*L/2)*cos(θ)
	]

	accel = A \ b

	xddot = accel[1]
	θddot = accel[2]

	du[1] = xdot
	du[2] = xddot
	du[3] = θdot
	du[4] = θddot
end

# ╔═╡ d25ee215-43c7-4b16-95eb-25703ae2aa59
begin
	p = (m1, m2, k, g, L, J2)

	prob = ODEProblem(eom!, u0, tspan, p)
	sol = solve(prob, Tsit5(), reltol=1e-8, abstol=1e-8)
end

# ╔═╡ 66f7cbe7-e9bd-42c6-9682-c8830595df43
plot(
	sol.t,
	[u[1] for u in sol.u],
	xlabel = "Time [s]",
	ylabel = "Block displacement x₁ [m]",
	label = "x₁(t)",
	linewidth = 2,
	grid = true
)

# ╔═╡ 81ef5fcb-d14f-4071-956b-c421783140c2
plot(
	sol.t,
	rad2deg.([u[3] for u in sol.u]),
	xlabel = "Time [s]",
	ylabel = "Bar angle θ₂ [deg]",
	label = "θ₂(t)",
	linewidth = 2,
	grid = true
)

# ╔═╡ 273d3f59-ed58-4bbb-82a4-4a716ec15a15
begin
	x1 = [u[1] for u in sol.u]
	θ2 = [u[3] for u in sol.u]

	x2 = x1 .+ (L/2).*cos.(θ2)
	y2 = (L/2).*sin.(θ2)

	plot(
		x2, y2,
		xlabel = "x₂ [m]",
		ylabel = "y₂ [m]",
		label = "bar center of mass path",
		linewidth = 2,
		aspect_ratio = :equal,
		grid = true
	)
end

# ╔═╡ 464909cb-d8cc-4c11-bd50-36dcd536ea17
begin
	anim = @animate for i in 1:5:length(sol.t)
		x = sol.u[i][1]
		θ = sol.u[i][3]

		pin_x = x
		pin_y = 0.0

		tip_x = pin_x + L*cos(θ)
		tip_y = pin_y + L*sin(θ)

		plot(
			[-0.5, 0.7], [0, 0],
			label = "",
			linewidth = 2,
			xlims = (-0.5, 0.7),
			ylims = (-0.25, 0.55),
			aspect_ratio = :equal,
			grid = true,
			title = "t = $(round(sol.t[i], digits=2)) s"
		)

		plot!([-0.4, x - 0.05], [0, 0], linewidth = 3, label = "spring")
		plot!([pin_x, tip_x], [pin_y, tip_y], linewidth = 5, label = "bar")
		scatter!([pin_x], [pin_y], markersize = 6, label = "pin")

		# block outline
		plot!(
			[x-0.05, x+0.05, x+0.05, x-0.05, x-0.05],
			[-0.05, -0.05, 0.05, 0.05, -0.05],
			linewidth = 2,
			label = "block"
		)
	end

	gif(anim, "spring_bar_motion.gif", fps = 20)
end

# ╔═╡ Cell order:
# ╠═0d0636e0-3225-4e3d-ba54-180e3f726f14
# ╠═984ec51b-75e6-4c0e-b8f5-78610983da1b
# ╠═c6ec4d9c-cb63-4c25-a5ac-5ae3af3ceb88
# ╠═d0c26027-876d-4426-8250-e3def024046d
# ╠═d25ee215-43c7-4b16-95eb-25703ae2aa59
# ╠═66f7cbe7-e9bd-42c6-9682-c8830595df43
# ╠═81ef5fcb-d14f-4071-956b-c421783140c2
# ╠═273d3f59-ed58-4bbb-82a4-4a716ec15a15
# ╠═464909cb-d8cc-4c11-bd50-36dcd536ea17
