### A Pluto.jl notebook ###
# v0.20.24

using Markdown
using InteractiveUtils

# ╔═╡ d0c26027-876d-4426-8250-e3def024046d
begin
	import Pkg   

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
## System Setup

**The system consists of:**
- A **sliding block** along a horizontal track
- A **linear spring** attached to the block that streches and compresses
- A **rigid bar (length = 0.4 m)** pinned to the block which acts as a compound pendulum

**2 DOF:**
- x(t): horizontal displacement of the block
- theta(t): angular position of the rigid bar

**Given Parameters:**

- Block mass: \(m_1 = 0.1\) kg  
- Bar mass: \(m_2 = 0.3\) kg  
- Spring stiffness: \(k = 10\) N/m  
- Gravity: \(g = 9.81\) m/s²  
- Bar length: \(L = 0.4\) m  
- Bar Moment of Inertia: $J_2 = \frac{1}{12} m_2 L^2$

"""

# ╔═╡ c139bc1d-a312-4d39-8453-948b73dcae38
md"""
## 1. Determine Constraint Equations


"""

# ╔═╡ d00e5659-8419-4cd1-a4b0-334de75beff3
md"""
**Equations of Motion:**

Two coupled nonlinear equations exist in this problem, both of which are derived from Newton–Euler mechanics. They can be solved using ODEs.

**Generalized Coordinates:**


$\mathbf{q} =
\begin{bmatrix}
x_1 \\
y_1 \\
\theta_1 \\
x_2 \\
y_2 \\
\theta_2
\end{bmatrix}$

The generalized coordinate matrix has two bodies (x_1 y_1 and theta_1) which is the sliding block and (x_2 y_2 theta_2) which is the rigid bar.

**Constraint Equations:**

Block is only constrained horizontally:

$C_1 = y_1 = 0$

$C_2 = \theta_1 = 0$

Bar is pinned to the block meaning that the block becomes a constraint on the bar.

$C_3 = x_1 - x_2 + \frac{L}{2}\cos\theta_2 = 0$

$C_4 = y_1 - y_2 + \frac{L}{2}\sin\theta_2 = 0$

The above constraint equations C_1 and C_2 from the block and C_3 and C_4 from the bar can be combined into the below matrix:

$\mathbf{C}(\mathbf{q}) =
\begin{bmatrix}
y_1 \\
\theta_1 \\
x_1 - x_2 + \frac{L}{2}\cos\theta_2 \\
y_1 - y_2 + \frac{L}{2}\sin\theta_2
\end{bmatrix}
=
\mathbf{0}$


**Note Reduced Coordinates:**

Due to the individual DOF of the block and the bar some coordinates can be eliminated.

$x = x_1$

$\theta = \theta_2$


"""

# ╔═╡ 75d68833-6313-4f71-b415-5da61c3b21ca
md"""
## 2. Create an augmented solution method for the dynamic motion of these two moving parts

"""

# ╔═╡ 5727812e-3fe1-4ab2-b8dd-c18b9699929a
md"""
**Constrained equations of motion in augmented form:**

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


$\mathbf{C_q} = \frac{\partial \mathbf{C}}{\partial \mathbf{q}}$


**Force introduced via $\boldsymbol{\lambda}$. Applied Forces (spring and gravity):**

Spring force on the block:

$F_s = -kx_1$

Gravity on block and bar:

$F_{g1} = -m_1g$

$F_{g2} = -m_2g$

**Force vector:**

$\mathbf{Q} =
\begin{bmatrix}
-kx_1 \\
-m_1g \\
0 \\
0 \\
-m_2g \\
0
\end{bmatrix}$



**Bar center of mass (COM) coordinate location:**

$x_2 = x + \frac{L}{2}\cos\theta$

$y_2 = \frac{L}{2}\sin\theta$





**Mass Matrix of the System:**


$\mathbf{M} =
\begin{bmatrix}
m_1 & 0 & 0 & 0 & 0 & 0 \\
0 & m_1 & 0 & 0 & 0 & 0 \\
0 & 0 & J_1 & 0 & 0 & 0 \\
0 & 0 & 0 & m_2 & 0 & 0 \\
0 & 0 & 0 & 0 & m_2 & 0 \\
0 & 0 & 0 & 0 & 0 & J_2
\end{bmatrix}$



**Bar Moment of Inertia:**

$J_2 = \frac{1}{12}m_2L^2$

$J_2 = \frac{1}{12}(0.3)(0.4)^2$

$J_2 = 0.004 \ \text{kg}\cdot\text{m}^2$




**Reduced Equations of Motion:**


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


**Convert into ODE:**

$\mathbf{u} = [x, \dot{x}, \theta, \dot{\theta}]$


$\mathbf{A} \begin{bmatrix} \ddot{x} \\ \ddot{\theta} \end{bmatrix} = \mathbf{b}$

"""

# ╔═╡ 301e0150-6013-4cfb-98f0-9b7169957eb9
md"""
## 3.  Visualize the motion of the system as the two parts complete at leastone oscillation

The numerical solution provides:

- Block displacement vs time
- Bar Angle Motion vs time
- Trajectory of the bar center of mass
- Animated motion of all componenents in the system


Initial Inputs:

- Initial block displacement: $x(0) = 0.05~m$  
- Initial block velocity: $\dot{x}(0) = 0$
- Initial angle: $\theta(0) = 60^\circ$
- Initial angular velocity: $\dot{\theta}(0) = 0$

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

# ╔═╡ 73de3209-ee3f-4a03-bc5f-12f4b03e0242
md"""
**Time vs. Block Displacement**

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

# ╔═╡ 0648c70f-b6d4-4194-9147-b5864a5acb39
md"""
**Time vs. Block Velocity**

"""

# ╔═╡ dba9e467-bc36-4abd-a904-c60ccbb73c1f
md"""
**Block Displacement vs. Block Velocity**

"""

# ╔═╡ 62df6878-575a-436f-b011-dbe77af2412d
md"""
**Time vs. Bar Angle**

"""

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

# ╔═╡ 7fb6d384-0266-493e-ad0f-8c8a3bd15044
md"""
**Time vs. Bar Angular Velocity**

"""

# ╔═╡ 51494a51-dc2d-41cd-bc90-cc97255e0506
md"""
**Bar Angle vs. Bar Angular Velocity**

"""

# ╔═╡ 50c16cb9-f124-46fb-8802-256121151ae0
md"""
**Bar Center of Mass Path**

"""

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

# ╔═╡ ffbedc3c-9651-4321-ac96-9690423f8106
md"""
**Bar Tip x vs. y**

"""

# ╔═╡ f7f0fd14-7dc0-4663-9d94-5f93f464459b
md"""
**Motion of the System**

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

# ╔═╡ b8ea6a32-dd64-473f-83d0-8e43f274e05c
plot(
	t, xd,
	xlabel = "Time [s]",
	ylabel = "Block velocity ẋ [m/s]",
	label = "ẋ(t)",
	linewidth = 2,
	grid = true
)

# ╔═╡ b074f0fe-82b1-4586-a678-eeecd3ef8520
plot(
	x, xd,
	xlabel = "Block displacement x [m]",
	ylabel = "Block velocity ẋ [m/s]",
	label = "Block phase plot",
	linewidth = 2,
	grid = true
)

# ╔═╡ 05d543ec-64dd-4dc2-847d-39faeb47f44f
plot(
	t, θd,
	xlabel = "Time [s]",
	ylabel = "Bar Angular velocity θ̇ [rad/s]",
	label = "θ̇(t)",
	linewidth = 2,
	grid = true
)

# ╔═╡ e8219de1-7d8f-422c-91e7-fc39d74cd7f5
plot(
	rad2deg.(θ), θd,
	xlabel = "Bar angle θ [deg]",
	ylabel = "Bar Angular velocity θ̇ [rad/s]",
	label = "Bar phase plot",
	linewidth = 2,
	grid = true
)

# ╔═╡ 93e31807-c170-473a-93d4-79c0268c25f6
plot(
	xtip, ytip,
	xlabel = "Bar Tip x-position [m]",
	ylabel = "Bar Tip y-position [m]",
	label = "Bar tip path",
	linewidth = 2,
	aspect_ratio = :equal,
	grid = true
)

# ╔═╡ 5bade53f-5cbb-479c-b2f5-b23811238217
md"""
**Energy of the System**

"""

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

# ╔═╡ 83d7baa0-ddda-4f0d-b89f-a5149dc03630
md"""
**Energy Exchange of the System**

"""

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

# ╔═╡ 5eb05d4b-de9b-4ed6-a96b-9eae9c8ac347
md"""
## 4. Calculate and show (graph or vectors) the constraint forces acting on the 2-body system


"""

# ╔═╡ 687c92e2-ae7a-4468-a81b-2b2c5aa0c108
md"""
**Force of Pin on Bar:**

$\mathbf{P} =
\begin{bmatrix}
F_{px} \\
F_{py}
\end{bmatrix}$

**Force Constraint Equations will will keep the bar and  block pinned together:**

$F_{px} = m_2 \ddot{x}_2$

$F_{py} - m_2g = m_2 \ddot{y}_2$

$F_{py} = m_2 \ddot{y}_2 + m_2g$

"""

# ╔═╡ 11c192e1-1819-4469-abd4-1c16d02061cc
begin
	# Recalculate accelerations from the solved motion
	xdd = similar(x)
	θdd = similar(θ)

	for i in eachindex(t)
		A = [
			m1 + m2                 -(m2*L/2)*sin(θ[i])
			-(m2*L/2)*sin(θ[i])    J2 + m2*L^2/4
		]

		b = [
			(m2*L/2)*cos(θ[i])*θd[i]^2 - k*x[i]
			-(m2*g*L/2)*cos(θ[i])
		]

		accel = A \ b

		xdd[i] = accel[1]
		θdd[i] = accel[2]
	end

	# Bar center of mass accelerations
	xbar_dd = xdd .- (L/2).*cos.(θ).*θd.^2 .- (L/2).*sin.(θ).*θdd
	ybar_dd = .-(L/2).*sin.(θ).*θd.^2 .+ (L/2).*cos.(θ).*θdd

	# Constraint forces at the pin acting on the bar
	F_px = m2 .* xbar_dd
	F_py = m2 .* ybar_dd .+ m2*g
end

# ╔═╡ 466a71ab-2372-4f25-957d-7c50aec8f6e3
md"""
**Time vs. Constraint Equation of Pin Force x and Pin Force y on the Bar**

"""

# ╔═╡ 5fddea89-6b79-4edc-b7db-7b0a62ac4cc9
begin
	plot(
		t, F_px,
		xlabel = "Time [s]",
		ylabel = "Constraint force [N]",
		label = "Pin force F_px on bar",
		linewidth = 2,
		grid = true
	)
	
	plot!(t, F_py, label = "Pin force F_py on bar", linewidth = 2)
end

# ╔═╡ 07d1ee00-41f6-463c-8e40-d42142956db0
md"""
**Bar Center of Mass Path x vs. y**

"""

# ╔═╡ bd49a38d-33d9-43b3-b68b-86dc2bf33354
begin
	skip = 1:15:length(t)

	plot(
		xbar, ybar,
		xlabel = "x-position [m]",
		ylabel = "y-position [m]",
		label = "bar center path",
		linewidth = 2,
		aspect_ratio = :equal,
		grid = true
	)

	quiver!(
		xbar[skip],
		ybar[skip],
		quiver = (0.03 .* F_px[skip], 0.03 .* F_py[skip]),
		label = "constraint force vectors"
	)
end

# ╔═╡ Cell order:
# ╠═7d034eaa-47dc-4e7e-977b-edc9cab409d5
# ╠═eb1f3898-7861-482f-bc20-0c8ff9cd403c
# ╠═c139bc1d-a312-4d39-8453-948b73dcae38
# ╠═d00e5659-8419-4cd1-a4b0-334de75beff3
# ╠═75d68833-6313-4f71-b415-5da61c3b21ca
# ╠═5727812e-3fe1-4ab2-b8dd-c18b9699929a
# ╠═301e0150-6013-4cfb-98f0-9b7169957eb9
# ╠═0d0636e0-3225-4e3d-ba54-180e3f726f14
# ╠═984ec51b-75e6-4c0e-b8f5-78610983da1b
# ╠═c6ec4d9c-cb63-4c25-a5ac-5ae3af3ceb88
# ╠═d0c26027-876d-4426-8250-e3def024046d
# ╠═d25ee215-43c7-4b16-95eb-25703ae2aa59
# ╠═73de3209-ee3f-4a03-bc5f-12f4b03e0242
# ╠═66f7cbe7-e9bd-42c6-9682-c8830595df43
# ╠═0648c70f-b6d4-4194-9147-b5864a5acb39
# ╠═b8ea6a32-dd64-473f-83d0-8e43f274e05c
# ╠═dba9e467-bc36-4abd-a904-c60ccbb73c1f
# ╠═b074f0fe-82b1-4586-a678-eeecd3ef8520
# ╠═62df6878-575a-436f-b011-dbe77af2412d
# ╠═81ef5fcb-d14f-4071-956b-c421783140c2
# ╠═7fb6d384-0266-493e-ad0f-8c8a3bd15044
# ╠═05d543ec-64dd-4dc2-847d-39faeb47f44f
# ╠═51494a51-dc2d-41cd-bc90-cc97255e0506
# ╠═e8219de1-7d8f-422c-91e7-fc39d74cd7f5
# ╠═50c16cb9-f124-46fb-8802-256121151ae0
# ╠═273d3f59-ed58-4bbb-82a4-4a716ec15a15
# ╠═ffbedc3c-9651-4321-ac96-9690423f8106
# ╠═93e31807-c170-473a-93d4-79c0268c25f6
# ╠═f7f0fd14-7dc0-4663-9d94-5f93f464459b
# ╠═464909cb-d8cc-4c11-bd50-36dcd536ea17
# ╠═3d2a30b3-d262-4aad-9f50-c95167c5097c
# ╠═5bade53f-5cbb-479c-b2f5-b23811238217
# ╠═054e7496-5b28-4250-851a-40c6a4770d39
# ╠═83d7baa0-ddda-4f0d-b89f-a5149dc03630
# ╠═63109933-b837-4e4b-9380-b01ede925e56
# ╠═5eb05d4b-de9b-4ed6-a96b-9eae9c8ac347
# ╠═687c92e2-ae7a-4468-a81b-2b2c5aa0c108
# ╠═11c192e1-1819-4469-abd4-1c16d02061cc
# ╠═466a71ab-2372-4f25-957d-7c50aec8f6e3
# ╠═5fddea89-6b79-4edc-b7db-7b0a62ac4cc9
# ╠═07d1ee00-41f6-463c-8e40-d42142956db0
# ╠═bd49a38d-33d9-43b3-b68b-86dc2bf33354
