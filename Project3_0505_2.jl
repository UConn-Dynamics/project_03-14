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

# ╔═╡ 7d034eaa-47dc-4e7e-977b-edc9cab409d5
md"""
# Project_03 MBD Compound Pendulum

This project involves a rigid bar which is connected to a sliding block along a horizontal tracks. The sliding block is connected to a spring that stretches and compresses. The rigid bar L=0.4 m acts as a compound pendulum.
"""

# ╔═╡ eb1f3898-7861-482f-bc20-0c8ff9cd403c
md"""
## System Description

The system consists of:
- A **sliding block** along a horizontal track
- A **linear spring** attached to the block that streches and compresses
- A **rigid bar (length \(L = 0.4\) m)** pinned to the block which acts as a compound pendulum

There are two degrees of freedom in this example:
- \(x(t)\): horizontal displacement of the block
- \(\theta(t)\): angular position of the rigid bar

"""

# ╔═╡ c139bc1d-a312-4d39-8453-948b73dcae38
md"""
## 1. Determine Constraint Equations


"""

# ╔═╡ 63884579-742b-4a4f-ad36-c0127fca0d72
md"""
## System Parameters

The following physical parameters are used:

- Block mass: \(m_1 = 0.1\) kg  
- Bar mass: \(m_2 = 0.3\) kg  
- Spring stiffness: \(k = 10\) N/m  
- Gravity: \(g = 9.81\) m/s²  
- Bar length: \(L = 0.4\) m  

The bar is modeled as a slender rod with moment of inertia:

$J_2 = \frac{1}{12} m_2 L^2$
"""

# ╔═╡ d00e5659-8419-4cd1-a4b0-334de75beff3
md"""
## Equations of Motion

Two coupled nonlinear equations exist in this problem, both of which are derived from Newton–Euler mechanics:

"""

# ╔═╡ 1865b0c2-ef31-4a85-8e1f-3b5e97f2cdf0
md"""
These equations are solved using an ODE solver.
"""

# ╔═╡ 003d6873-3443-4d50-b0bc-8479214299a6
md"""
## 2. Create an augmented solution method for the dynamic motion of these two moving parts

"""

# ╔═╡ 2b126485-feb5-4cd3-9a16-a51717d345b3
md"""
Generalized Coordinates

The two-body system is described using the generalized coordinate vector

$\mathbf{q} =
\begin{bmatrix}
x_1 \\
y_1 \\
\theta_1 \\
x_2 \\
y_2 \\
\theta_2
\end{bmatrix}$

where body 1 is the sliding block and body 2 is the rigid bar.
"""

# ╔═╡ a4b6e722-7a89-4a32-aab1-c5bc2e075f35
md"""
Constraint Equations

The block is constrained to move only in the horizontal direction, so

$C_1 = y_1 = 0$

$C_2 = \theta_1 = 0$

The rigid bar is pinned to the block. Therefore, the pin location on the bar must remain coincident with the block location:

$C_3 = x_1 - x_2 + \frac{L}{2}\cos\theta_2 = 0$

$C_4 = y_1 - y_2 + \frac{L}{2}\sin\theta_2 = 0$

The full constraint vector is

$\mathbf{C}(\mathbf{q}) =
\begin{bmatrix}
y_1 \\
\theta_1 \\
x_1 - x_2 + \frac{L}{2}\cos\theta_2 \\
y_1 - y_2 + \frac{L}{2}\sin\theta_2
\end{bmatrix}
=
\mathbf{0}$
"""

# ╔═╡ ebf99e0d-548b-4f1f-aeb8-d85aa246b9bb
md"""
Applied Forces

The spring force acts horizontally on the block:

$F_s = -kx_1$

Gravity acts downward on both bodies:

$F_{g1} = -m_1g$

$F_{g2} = -m_2g$

The generalized force vector is

$\mathbf{Q} =
\begin{bmatrix}
-kx_1 \\
-m_1g \\
0 \\
0 \\
-m_2g \\
0
\end{bmatrix}$
"""

# ╔═╡ 59efded2-c383-4c15-9361-eaf62a6a2959
md"""
Mass Matrix

The mass matrix for the two-body system is

$\mathbf{M} =
\begin{bmatrix}
m_1 & 0 & 0 & 0 & 0 & 0 \\
0 & m_1 & 0 & 0 & 0 & 0 \\
0 & 0 & J_1 & 0 & 0 & 0 \\
0 & 0 & 0 & m_2 & 0 & 0 \\
0 & 0 & 0 & 0 & m_2 & 0 \\
0 & 0 & 0 & 0 & 0 & J_2
\end{bmatrix}$

Since the bar is modeled as a slender rod,

$J_2 = \frac{1}{12}m_2L^2$

Using \(m_2 = 0.3\) kg and \(L = 0.4\) m,

$J_2 = \frac{1}{12}(0.3)(0.4)^2$

$J_2 = 0.004 \ \text{kg}\cdot\text{m}^2$
"""

# ╔═╡ 5727812e-3fe1-4ab2-b8dd-c18b9699929a
md"""
Augmented Equations of Motion

The constrained equations of motion are written in augmented matrix form:

$\begin{bmatrix}
\mathbf{M} & \mathbf{C_q}^T \\
\mathbf{C_q} & \mathbf{0}
\end{bmatrix}
\begin{bmatrix}
\ddot{\mathbf{q}} \\
\boldsymbol{\lambda}
\end{bmatrix}
=
\begin{bmatrix}
\mathbf{Q} \\
-\dot{\mathbf{C_q}}\dot{\mathbf{q}}
\end{bmatrix}$

where

$\mathbf{C_q} = \frac{\partial \mathbf{C}}{\partial \mathbf{q}}$

and \(\boldsymbol{\lambda}\) represents the Lagrange multipliers, which are related to the constraint forces.
"""

# ╔═╡ 79d2c46a-05b3-433c-8d27-723509f1b400
md"""
Reduced Coordinates

Because the block only slides horizontally and the bar rotates about the pin, the system can also be described using two coordinates:

$x = x_1$

$\theta = \theta_2$

The bar center of mass is located at

$x_2 = x + \frac{L}{2}\cos\theta$

$y_2 = \frac{L}{2}\sin\theta$
"""

# ╔═╡ cc85b068-0aa9-4088-8498-4c0e98137b85
md"""
Reduced Equations of Motion

The reduced equations of motion are

$(m_1+m_2)\ddot{x}
-\frac{m_2L}{2}
\left(
\cos\theta \dot{\theta}^2
+
\sin\theta \ddot{\theta}
\right)
+kx=0$

$\left(J_2+\frac{m_2L^2}{4}\right)\ddot{\theta}
-\frac{m_2L}{2}\sin\theta \ddot{x}
+\frac{m_2gL}{2}\cos\theta=0$
"""

# ╔═╡ 08f4e782-923f-40f6-ace8-2e0c23a08431
md"""
Matrix Form Used in Julia

The equations are solved in Julia using the matrix form

$\mathbf{A}
\begin{bmatrix}
\ddot{x} \\
\ddot{\theta}
\end{bmatrix}
=
\mathbf{b}$

where

$\mathbf{A} =
\begin{bmatrix}
m_1+m_2 & -\frac{m_2L}{2}\sin\theta \\
-\frac{m_2L}{2}\sin\theta & J_2+\frac{m_2L^2}{4}
\end{bmatrix}$

and

$\mathbf{b} =
\begin{bmatrix}
\frac{m_2L}{2}\cos\theta \dot{\theta}^2-kx \\
-\frac{m_2gL}{2}\cos\theta
\end{bmatrix}$

At each time step, Julia solves

\$\begin{bmatrix}
\ddot{x} \\
\ddot{\theta}
\end{bmatrix}
=
\mathbf{A}^{-1}\mathbf{b}$
"""

# ╔═╡ 240bba9f-927d-438c-a111-57901b7fce4a
md"""
## Numerical Solution Method

The system is converted into a first-order system of ordinary differential equations:

$\mathbf{u} = [x, \dot{x}, \theta, \dot{\theta}]$

A matrix system is formed:

$\mathbf{A} \begin{bmatrix} \ddot{x} \\ \ddot{\theta} \end{bmatrix} = \mathbf{b}$

which is solved at each time step.

The equations are integrated using the **Runge–Kutta method** 
"""

# ╔═╡ 49e42a42-7c68-4a87-ac49-e5abd63b529b
md"""
## Initial Conditions

The system is initialized as:

- Initial block displacement: $x(0) = 0.05~m$  
- Initial block velocity: $\dot{x}(0) = 0$
- Initial angle: $\theta(0) = 60^\circ$
- Initial angular velocity: $\dot{\theta}(0) = 0$

These conditions produce oscillatory motion of both the block and the pendulum.
"""

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

# ╔═╡ ba203ffa-22ac-4bc1-9e44-28fadef2894f
md"""
## 3.  Visualize the motion of the system as the two parts complete at leastone oscillation

The numerical solution provides:

- Block displacement vs time
- Bar Angle Motion vs time
- Trajectory of the bar center of mass
- Animated motion of all componenents in the system

"""

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

# ╔═╡ a15f2328-aada-4be6-9988-970d70b08e50
md"""
## 4. Calculate and show (graph or vectors) the constraint forces acting on the 2-body system


"""

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

# ╔═╡ 3d2a30b3-d262-4aad-9f50-c95167c5097c
begin
	x  = [u[1] for u in sol.u]
	xd = [u[2] for u in sol.u]
	θ  = [u[3] for u in sol.u]
	θd = [u[4] for u in sol.u]
	t  = sol.t

	xbar = x .+ (L/2).*cos.(θ)
	ybar = (L/2).*sin.(θ)

	xtip = x .+ L.*cos.(θ)
	ytip = L.*sin.(θ)
end

# ╔═╡ c0039817-3304-499d-9dc6-392d961d619a
plot(
	t, xd,
	xlabel = "Time [s]",
	ylabel = "Block velocity ẋ [m/s]",
	label = "ẋ(t)",
	linewidth = 2,
	grid = true
)

# ╔═╡ af5d2a2f-e6b6-457e-a488-27e2b48fd07e
plot(
	t, θd,
	xlabel = "Time [s]",
	ylabel = "Angular velocity θ̇ [rad/s]",
	label = "θ̇(t)",
	linewidth = 2,
	grid = true
)

# ╔═╡ 385f2f3a-5d5a-4ba8-ae95-9ed55e09cd03
plot(
	x, xd,
	xlabel = "Block displacement x [m]",
	ylabel = "Block velocity ẋ [m/s]",
	label = "Block phase plot",
	linewidth = 2,
	grid = true
)

# ╔═╡ f3615ba8-7087-4395-a792-de70031c2aa1
plot(
	rad2deg.(θ), θd,
	xlabel = "Bar angle θ [deg]",
	ylabel = "Angular velocity θ̇ [rad/s]",
	label = "Bar phase plot",
	linewidth = 2,
	grid = true
)

# ╔═╡ 195d381d-6995-4deb-af37-d3a36a0eceee
plot(
	xtip, ytip,
	xlabel = "Tip x-position [m]",
	ylabel = "Tip y-position [m]",
	label = "Bar tip path",
	linewidth = 2,
	aspect_ratio = :equal,
	grid = true
)

# ╔═╡ 054e7496-5b28-4250-851a-40c6a4770d39
begin
	spring_PE = 0.5 .* k .* x.^2
	bar_PE = m2 .* g .* ybar

	block_KE = 0.5 .* m1 .* xd.^2
	bar_trans_KE = 0.5 .* m2 .* ((xd .- (L/2).*sin.(θ).*θd).^2 .+ ((L/2).*cos.(θ).*θd).^2)
	bar_rot_KE = 0.5 .* J2 .* θd.^2

	total_KE = block_KE .+ bar_trans_KE .+ bar_rot_KE
	total_PE = spring_PE .+ bar_PE
	total_E = total_KE .+ total_PE

	plot(
		t, total_KE,
		label = "Kinetic Energy",
		linewidth = 2,
		xlabel = "Time [s]",
		ylabel = "Energy [J]",
		title = "Energy of the System",
		grid = true
	)

	plot!(t, total_PE, label = "Potential Energy", linewidth = 2)
	plot!(t, total_E, label = "Total Energy", linewidth = 2)
end

# ╔═╡ 63109933-b837-4e4b-9380-b01ede925e56
begin
	anim_energy = @animate for i in 2:5:length(t)
		plot(
			t[1:i], total_KE[1:i],
			label = "Kinetic Energy",
			linewidth = 2,
			xlims = (minimum(t), maximum(t)),
			ylims = (0, maximum(total_E)*1.1),
			xlabel = "Time [s]",
			ylabel = "Energy [J]",
			title = "Energy Exchange",
			grid = true
		)

		plot!(t[1:i], total_PE[1:i], label = "Potential Energy", linewidth = 2)
		plot!(t[1:i], total_E[1:i], label = "Total Energy", linewidth = 2)
	end

	gif(anim_energy, "energy_exchange.gif", fps = 20)
end

# ╔═╡ Cell order:
# ╠═7d034eaa-47dc-4e7e-977b-edc9cab409d5
# ╠═eb1f3898-7861-482f-bc20-0c8ff9cd403c
# ╠═c139bc1d-a312-4d39-8453-948b73dcae38
# ╠═63884579-742b-4a4f-ad36-c0127fca0d72
# ╠═d00e5659-8419-4cd1-a4b0-334de75beff3
# ╠═1865b0c2-ef31-4a85-8e1f-3b5e97f2cdf0
# ╠═003d6873-3443-4d50-b0bc-8479214299a6
# ╠═2b126485-feb5-4cd3-9a16-a51717d345b3
# ╠═a4b6e722-7a89-4a32-aab1-c5bc2e075f35
# ╠═ebf99e0d-548b-4f1f-aeb8-d85aa246b9bb
# ╠═59efded2-c383-4c15-9361-eaf62a6a2959
# ╠═5727812e-3fe1-4ab2-b8dd-c18b9699929a
# ╠═79d2c46a-05b3-433c-8d27-723509f1b400
# ╠═cc85b068-0aa9-4088-8498-4c0e98137b85
# ╠═08f4e782-923f-40f6-ace8-2e0c23a08431
# ╠═240bba9f-927d-438c-a111-57901b7fce4a
# ╠═49e42a42-7c68-4a87-ac49-e5abd63b529b
# ╠═0d0636e0-3225-4e3d-ba54-180e3f726f14
# ╠═984ec51b-75e6-4c0e-b8f5-78610983da1b
# ╠═c6ec4d9c-cb63-4c25-a5ac-5ae3af3ceb88
# ╠═d0c26027-876d-4426-8250-e3def024046d
# ╠═d25ee215-43c7-4b16-95eb-25703ae2aa59
# ╠═ba203ffa-22ac-4bc1-9e44-28fadef2894f
# ╠═66f7cbe7-e9bd-42c6-9682-c8830595df43
# ╠═81ef5fcb-d14f-4071-956b-c421783140c2
# ╠═273d3f59-ed58-4bbb-82a4-4a716ec15a15
# ╠═a15f2328-aada-4be6-9988-970d70b08e50
# ╠═464909cb-d8cc-4c11-bd50-36dcd536ea17
# ╠═3d2a30b3-d262-4aad-9f50-c95167c5097c
# ╠═c0039817-3304-499d-9dc6-392d961d619a
# ╠═af5d2a2f-e6b6-457e-a488-27e2b48fd07e
# ╠═385f2f3a-5d5a-4ba8-ae95-9ed55e09cd03
# ╠═f3615ba8-7087-4395-a792-de70031c2aa1
# ╠═195d381d-6995-4deb-af37-d3a36a0eceee
# ╠═054e7496-5b28-4250-851a-40c6a4770d39
# ╠═63109933-b837-4e4b-9380-b01ede925e56
