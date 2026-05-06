### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ 7e6ddfa3-2d84-4c34-81bd-f0ff0e561001
begin
    using LinearAlgebra
    using Printf
end


# ╔═╡ 7e6ddfa3-2d84-4c34-81bd-f0ff0e561003
begin
    struct SVGBlock
        data::String
    end

    struct HTMLBlock
        data::String
    end

    Base.show(io::IO, ::MIME"image/svg+xml", svg::SVGBlock) = print(io, svg.data)
    Base.show(io::IO, ::MIME"text/html", html::HTMLBlock) = print(io, html.data)

    load_svg(path::AbstractString) = SVGBlock(read(path, String))
    load_html(path::AbstractString) = HTMLBlock(read(path, String))

    struct ModelParameters
        m1::Float64
        m2::Float64
        L::Float64
        k::Float64
        g::Float64
        J1::Float64
        J2::Float64
        alpha::Float64
        beta::Float64
        dt::Float64
        t_end::Float64
        block_width::Float64
        block_height::Float64
        wheel_radius::Float64
    end

    struct SimulationResult
        t::Vector{Float64}
        q::Matrix{Float64}
        v::Matrix{Float64}
        a::Matrix{Float64}
        lambda::Matrix{Float64}
        constraint_error::Vector{Float64}
        energy::Vector{Float64}
    end

    function default_parameters(; dt::Float64 = 1e-3, t_end::Float64 = 8.0)
        block_width = 0.10
        block_height = 0.06
        wheel_radius = 0.014
        m1 = 0.1
        m2 = 0.3
        L = 0.4
        k = 10.0
        g = 9.81
        J1 = m1 * (block_width^2 + block_height^2) / 12
        J2 = m2 * L^2 / 12
        alpha = 14.0
        beta = 14.0
        return ModelParameters(
            m1,
            m2,
            L,
            k,
            g,
            J1,
            J2,
            alpha,
            beta,
            dt,
            t_end,
            block_width,
            block_height,
            wheel_radius,
        )
    end

    mass_matrix(p::ModelParameters) = Diagonal([p.m1, p.m1, p.J1, p.m2, p.m2, p.J2])

    function make_initial_state(
        p::ModelParameters;
        x1::Float64 = 0.08,
        theta2::Float64 = deg2rad(55.0),
        x1dot::Float64 = 0.0,
        theta2dot::Float64 = 0.0,
    )
        r = 0.5 * p.L
        q = zeros(6)
        q[1] = x1
        q[2] = 0.0
        q[3] = 0.0
        q[4] = x1 + r * cos(theta2)
        q[5] = r * sin(theta2)
        q[6] = theta2

        v = zeros(6)
        v[1] = x1dot
        v[4] = x1dot - r * sin(theta2) * theta2dot
        v[5] = r * cos(theta2) * theta2dot
        v[6] = theta2dot
        return q, v
    end

    function constraint_vector(q::AbstractVector, p::ModelParameters)
        r = 0.5 * p.L
        theta2 = q[6]
        return [
            q[2]
            q[3]
            q[4] - q[1] - r * cos(theta2)
            q[5] - q[2] - r * sin(theta2)
        ]
    end

    function constraint_jacobian(q::AbstractVector, p::ModelParameters)
        r = 0.5 * p.L
        theta2 = q[6]
        return [
             0.0  1.0  0.0  0.0  0.0  0.0
             0.0  0.0  1.0  0.0  0.0  0.0
            -1.0  0.0  0.0  1.0  0.0  r * sin(theta2)
             0.0 -1.0  0.0  0.0  1.0 -r * cos(theta2)
        ]
    end

    function constraint_rhs(q::AbstractVector, v::AbstractVector, p::ModelParameters)
        r = 0.5 * p.L
        theta2 = q[6]
        omega2 = v[6]
        gamma = [
            0.0
            0.0
            -r * cos(theta2) * omega2^2
            -r * sin(theta2) * omega2^2
        ]
        c = constraint_vector(q, p)
        cdot = constraint_jacobian(q, p) * v
        return gamma - 2.0 * p.alpha * cdot - p.beta^2 * c
    end

    function generalized_forces(q::AbstractVector, p::ModelParameters)
        return [
            -p.k * q[1]
            -p.m1 * p.g
            0.0
            0.0
            -p.m2 * p.g
            0.0
        ]
    end

    function solve_augmented_system(q::AbstractVector, v::AbstractVector, p::ModelParameters)
        M = mass_matrix(p)
        Cq = constraint_jacobian(q, p)
        Z = zeros(size(Cq, 1), size(Cq, 1))
        A = [
            Matrix(M) Matrix(Cq')
            Cq        Z
        ]
        rhs = [
            generalized_forces(q, p)
            constraint_rhs(q, v, p)
        ]
        solution = A \ rhs
        return solution[1:6], solution[7:end]
    end

    function project_coordinates!(q::AbstractVector, p::ModelParameters)
        r = 0.5 * p.L
        q[2] = 0.0
        q[3] = 0.0
        q[4] = q[1] + r * cos(q[6])
        q[5] = q[2] + r * sin(q[6])
        return q
    end

    function project_velocities!(q::AbstractVector, v::AbstractVector, p::ModelParameters)
        r = 0.5 * p.L
        theta2 = q[6]
        omega2 = v[6]
        v[2] = 0.0
        v[3] = 0.0
        v[4] = v[1] - r * sin(theta2) * omega2
        v[5] = v[2] + r * cos(theta2) * omega2
        return v
    end

    function constrained_dynamics(q::AbstractVector, v::AbstractVector, p::ModelParameters)
        q_proj = copy(q)
        v_proj = copy(v)
        project_coordinates!(q_proj, p)
        project_velocities!(q_proj, v_proj, p)
        a, lambda = solve_augmented_system(q_proj, v_proj, p)
        return v_proj, a, lambda
    end

    function rk4_step(q::AbstractVector, v::AbstractVector, p::ModelParameters, dt::Float64)
        k1q, k1v, _ = constrained_dynamics(q, v, p)

        q2 = q .+ 0.5 * dt .* k1q
        v2 = v .+ 0.5 * dt .* k1v
        k2q, k2v, _ = constrained_dynamics(q2, v2, p)

        q3 = q .+ 0.5 * dt .* k2q
        v3 = v .+ 0.5 * dt .* k2v
        k3q, k3v, _ = constrained_dynamics(q3, v3, p)

        q4 = q .+ dt .* k3q
        v4 = v .+ dt .* k3v
        k4q, k4v, _ = constrained_dynamics(q4, v4, p)

        q_next = q .+ (dt / 6.0) .* (k1q .+ 2.0 .* k2q .+ 2.0 .* k3q .+ k4q)
        v_next = v .+ (dt / 6.0) .* (k1v .+ 2.0 .* k2v .+ 2.0 .* k3v .+ k4v)
        project_coordinates!(q_next, p)
        project_velocities!(q_next, v_next, p)
        return q_next, v_next
    end

    function total_energy(q::AbstractVector, v::AbstractVector, p::ModelParameters)
        M = mass_matrix(p)
        kinetic = 0.5 * dot(v, M * v)
        spring = 0.5 * p.k * q[1]^2
        gravity = p.m1 * p.g * q[2] + p.m2 * p.g * q[5]
        return kinetic + spring + gravity
    end

    function simulate_system(
        p::ModelParameters = default_parameters();
        x1::Float64 = 0.08,
        theta2::Float64 = deg2rad(55.0),
        x1dot::Float64 = 0.0,
        theta2dot::Float64 = 0.0,
    )
        q0, v0 = make_initial_state(p; x1 = x1, theta2 = theta2, x1dot = x1dot, theta2dot = theta2dot)
        n_steps = Int(round(p.t_end / p.dt))
        t = collect(range(0.0, step = p.dt, length = n_steps + 1))
        q_hist = zeros(6, n_steps + 1)
        v_hist = zeros(6, n_steps + 1)
        a_hist = zeros(6, n_steps + 1)
        lambda_hist = zeros(4, n_steps + 1)
        c_hist = zeros(n_steps + 1)
        e_hist = zeros(n_steps + 1)

        q = copy(q0)
        v = copy(v0)

        for i in eachindex(t)
            a, lambda = solve_augmented_system(q, v, p)
            q_hist[:, i] = q
            v_hist[:, i] = v
            a_hist[:, i] = a
            lambda_hist[:, i] = lambda
            c_hist[i] = norm(constraint_vector(q, p))
            e_hist[i] = total_energy(q, v, p)

            if i < length(t)
                q, v = rk4_step(q, v, p, p.dt)
            end
        end

        return SimulationResult(t, q_hist, v_hist, a_hist, lambda_hist, c_hist, e_hist)
    end

    function generalized_constraint_forces(result::SimulationResult, p::ModelParameters)
        n = length(result.t)
        qc = zeros(6, n)
        for i in 1:n
            Cq = constraint_jacobian(view(result.q, :, i), p)
            qc[:, i] = Cq' * view(result.lambda, :, i)
        end
        return qc
    end

    function padded_limits(values::AbstractVector; frac::Float64 = 0.08)
        vmin = minimum(values)
        vmax = maximum(values)
        if isapprox(vmin, vmax; atol = 1e-12, rtol = 1e-12)
            pad = max(1.0, abs(vmin) * frac + 1e-3)
        else
            pad = frac * (vmax - vmin)
        end
        return vmin - pad, vmax + pad
    end

    function svg_polyline(x::AbstractVector, y::AbstractVector, xmap, ymap)
        points = String[]
        for (xv, yv) in zip(x, y)
            push!(points, @sprintf("%.2f,%.2f", xmap(xv), ymap(yv)))
        end
        return join(points, " ")
    end


    function write_line_plot_svg(path::AbstractString, x::AbstractVector, series; title, xlabel, ylabel)
        width = 940
        height = 560
        margin_left = 90
        margin_right = 210
        margin_top = 70
        margin_bottom = 80
        plot_width = width - margin_left - margin_right
        plot_height = height - margin_top - margin_bottom

        all_y = reduce(vcat, [collect(item.values) for item in series])
        xmin, xmax = padded_limits(x; frac = 0.02)
        ymin, ymax = padded_limits(all_y)

        xmap(v) = margin_left + (v - xmin) / (xmax - xmin) * plot_width
        ymap(v) = margin_top + plot_height - (v - ymin) / (ymax - ymin) * plot_height

        open(path, "w") do io
            println(io, """<svg xmlns="http://www.w3.org/2000/svg" width="$width" height="$height" viewBox="0 0 $width $height">""")
            println(io, """<rect x="0" y="0" width="$width" height="$height" fill="white"/>""")
            println(io, """<text x="$(width / 2)" y="34" text-anchor="middle" font-size="26" font-family="Times New Roman">$title</text>""")
            println(io, """<line x1="$margin_left" y1="$(height - margin_bottom)" x2="$(width - margin_right)" y2="$(height - margin_bottom)" stroke="black" stroke-width="2"/>""")
            println(io, """<line x1="$margin_left" y1="$margin_top" x2="$margin_left" y2="$(height - margin_bottom)" stroke="black" stroke-width="2"/>""")

            for tick in 0:5
                xv = xmin + tick * (xmax - xmin) / 5
                px = xmap(xv)
                println(io, """<line x1="$px" y1="$(height - margin_bottom)" x2="$px" y2="$(height - margin_bottom + 8)" stroke="black" stroke-width="1.5"/>""")
                println(io, """<text x="$px" y="$(height - margin_bottom + 30)" text-anchor="middle" font-size="15" font-family="Arial">$(round(xv, digits = 2))</text>""")
                println(io, """<line x1="$px" y1="$margin_top" x2="$px" y2="$(height - margin_bottom)" stroke="#d8d8d8" stroke-width="1"/>""")
            end

            for tick in 0:5
                yv = ymin + tick * (ymax - ymin) / 5
                py = ymap(yv)
                println(io, """<line x1="$(margin_left - 8)" y1="$py" x2="$margin_left" y2="$py" stroke="black" stroke-width="1.5"/>""")
                println(io, """<text x="$(margin_left - 14)" y="$(py + 5)" text-anchor="end" font-size="15" font-family="Arial">$(round(yv, digits = 3))</text>""")
                println(io, """<line x1="$margin_left" y1="$py" x2="$(width - margin_right)" y2="$py" stroke="#d8d8d8" stroke-width="1"/>""")
            end

            for item in series
                points = svg_polyline(x, item.values, xmap, ymap)
                println(io, """<polyline fill="none" stroke="$(item.color)" stroke-width="2.5" points="$points"/>""")
            end

            println(io, """<text x="$(margin_left + plot_width / 2)" y="$(height - 20)" text-anchor="middle" font-size="20" font-family="Arial">$xlabel</text>""")
            println(io, """<text x="24" y="$(margin_top + plot_height / 2)" text-anchor="middle" font-size="20" font-family="Arial" transform="rotate(-90 24 $(margin_top + plot_height / 2))">$ylabel</text>""")

            legend_x = width - margin_right + 28
            legend_y = margin_top + 20
            for (index, item) in enumerate(series)
                y_line = legend_y + 28 * (index - 1)
                println(io, """<line x1="$legend_x" y1="$y_line" x2="$(legend_x + 32)" y2="$y_line" stroke="$(item.color)" stroke-width="3"/>""")
                println(io, """<text x="$(legend_x + 42)" y="$(y_line + 5)" font-size="16" font-family="Arial">$(item.label)</text>""")
            end

            println(io, "</svg>")
        end
    end

    function js_array(values::AbstractVector)
        return "[" * join((@sprintf("%.10f", value) for value in values), ",") * "]"
    end

    function write_animation_html(path::AbstractString, result::SimulationResult, params::ModelParameters)
        stride = 10
        indices = 1:stride:length(result.t)
        times = result.t[indices]
        x1 = result.q[1, indices]
        theta2 = result.q[6, indices]
        x_left = minimum(result.q[1, :]) - 0.32
        x_right = maximum(result.q[1, :]) + params.L + 0.15
        y_bottom = -0.13
        y_top = params.L + 0.12

        html = """
<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8" />
  <title>Project 03 Motion</title>
  <style>
    body { font-family: Arial, sans-serif; margin: 12px; background: #fafafa; color: #222; }
    #wrap { max-width: 980px; }
    canvas { background: white; border: 1px solid #ccc; display: block; }
    .controls { margin-top: 14px; display: flex; gap: 12px; align-items: center; flex-wrap: wrap; }
    button { padding: 8px 14px; font-size: 14px; }
    input[type="range"] { width: 420px; }
    .note { margin-top: 10px; color: #555; line-height: 1.5; }
  </style>
</head>
<body>
<div id="wrap">
  <canvas id="canvas" width="920" height="560"></canvas>
  <div class="controls">
    <button id="toggle">Pause</button>
    <label>Frame <input id="frame" type="range" min="0" max="$(length(times) - 1)" value="0" /></label>
    <span id="stamp"></span>
  </div>
  <div class="note">
    Slider motion, spring extension, and rigid-bar rotation are drawn directly from the constrained dynamics solution.
  </div>
</div>
<script>
const times = $(js_array(times));
const x1 = $(js_array(x1));
const theta = $(js_array(theta2));
const params = { L: $(params.L), bw: $(params.block_width), bh: $(params.block_height), wr: $(params.wheel_radius) };
const world = { xmin: $(x_left), xmax: $(x_right), ymin: $(y_bottom), ymax: $(y_top) };
const canvas = document.getElementById("canvas");
const ctx = canvas.getContext("2d");
const slider = document.getElementById("frame");
const stamp = document.getElementById("stamp");
const button = document.getElementById("toggle");
let playing = true;
let frame = 0;

function mapX(x) {
  const pad = 50;
  return pad + (x - world.xmin) / (world.xmax - world.xmin) * (canvas.width - 2 * pad);
}
function mapY(y) {
  const padTop = 40;
  const padBottom = 50;
  return padTop + (world.ymax - y) / (world.ymax - world.ymin) * (canvas.height - padTop - padBottom);
}
function line(x1, y1, x2, y2, width = 2) {
  ctx.lineWidth = width;
  ctx.beginPath();
  ctx.moveTo(x1, y1);
  ctx.lineTo(x2, y2);
  ctx.stroke();
}
function circle(x, y, r, fill = false) {
  ctx.beginPath();
  ctx.arc(x, y, r, 0, Math.PI * 2);
  if (fill) {
    ctx.fill();
  } else {
    ctx.stroke();
  }
}
function drawSpring(xStart, xEnd, y) {
  const n = 13;
  ctx.beginPath();
  ctx.moveTo(mapX(xStart), mapY(y));
  for (let i = 1; i <= 2 * n; i++) {
    const x = xStart + (xEnd - xStart) * i / (2 * n + 1);
    const yy = y + (i % 2 === 0 ? 0.02 : -0.02);
    ctx.lineTo(mapX(x), mapY(yy));
  }
  ctx.lineTo(mapX(xEnd), mapY(y));
  ctx.stroke();
}
function drawFrame(idx) {
  const xb = x1[idx];
  const th = theta[idx];
  const xPin = xb;
  const yPin = 0.0;
  const xTip = xPin + params.L * Math.cos(th);
  const yTip = yPin + params.L * Math.sin(th);
  ctx.clearRect(0, 0, canvas.width, canvas.height);
  ctx.strokeStyle = "#000";
  ctx.fillStyle = "#000";
  ctx.lineWidth = 2;

  line(mapX(world.xmin), mapY(-0.08), mapX(world.xmax), mapY(-0.08), 3);
  for (let i = 0; i < 18; i++) {
    const x = world.xmin + (world.xmax - world.xmin) * i / 17;
    line(mapX(x - 0.02), mapY(-0.10), mapX(x), mapY(-0.08), 1);
  }

  const wallX = world.xmin + 0.05;
  line(mapX(wallX), mapY(-0.02), mapX(wallX), mapY(0.12), 3);
  for (let i = 0; i < 5; i++) {
    const y = -0.01 + 0.03 * i;
    line(mapX(wallX - 0.035), mapY(y - 0.015), mapX(wallX), mapY(y), 1);
  }

  drawSpring(wallX, xb - 0.5 * params.bw, 0.0);

  const left = xb - 0.5 * params.bw;
  const right = xb + 0.5 * params.bw;
  const top = 0.5 * params.bh;
  const bottom = -0.5 * params.bh;
  ctx.strokeRect(mapX(left), mapY(top), mapX(right) - mapX(left), mapY(bottom) - mapY(top));

  circle(mapX(xb - 0.22 * params.bw), mapY(bottom - params.wr), mapX(xb + params.wr) - mapX(xb), false);
  circle(mapX(xb + 0.22 * params.bw), mapY(bottom - params.wr), mapX(xb + params.wr) - mapX(xb), false);

  line(mapX(xPin), mapY(yPin), mapX(xTip), mapY(yTip), 5);
  circle(mapX(xPin), mapY(yPin), 5, true);
  circle(mapX(xTip), mapY(yTip), 4, true);

  ctx.font = "18px Arial";
  ctx.fillText(`t = \${times[idx].toFixed(2)} s`, 60, 28);
  stamp.textContent = `time = \${times[idx].toFixed(2)} s`;
  slider.value = idx;
}

slider.addEventListener("input", (event) => {
  frame = Number(event.target.value);
  drawFrame(frame);
});
button.addEventListener("click", () => {
  playing = !playing;
  button.textContent = playing ? "Pause" : "Play";
});

function animate() {
  if (playing) {
    frame = (frame + 1) % times.length;
    drawFrame(frame);
  }
  requestAnimationFrame(animate);
}

drawFrame(0);
requestAnimationFrame(animate);
</script>
</body>
</html>
"""

        write(path, html)
    end


    function generate_outputs(
        base_dir::AbstractString;
        params::ModelParameters,
        result::SimulationResult,
        constraint_forces::AbstractMatrix,
    )
        output_dir = joinpath(base_dir, "outputs")
        mkpath(output_dir)

        write_line_plot_svg(
            joinpath(output_dir, "block_position.svg"),
            result.t,
            [(values = result.q[1, :], label = "x1 (m)", color = "#2458A6")];
            title = "Block translation",
            xlabel = "Time (s)",
            ylabel = "x1 (m)",
        )

        write_line_plot_svg(
            joinpath(output_dir, "bar_angle.svg"),
            result.t,
            [(values = rad2deg.(result.q[6, :]), label = "theta2 (deg)", color = "#B84732")];
            title = "Bar rotation",
            xlabel = "Time (s)",
            ylabel = "theta2 (deg)",
        )

        write_line_plot_svg(
            joinpath(output_dir, "constraint_forces.svg"),
            result.t,
            [
                (values = constraint_forces[1, :], label = "Block Cx", color = "#2458A6"),
                (values = constraint_forces[2, :], label = "Block Cy", color = "#2F855A"),
                (values = constraint_forces[4, :], label = "Bar Cx", color = "#B84732"),
                (values = constraint_forces[5, :], label = "Bar Cy", color = "#7D4BA0"),
            ];
            title = "Generalized constraint forces",
            xlabel = "Time (s)",
            ylabel = "Force / moment",
        )

        write_line_plot_svg(
            joinpath(output_dir, "constraint_error.svg"),
            result.t,
            [(values = log10.(result.constraint_error .+ 1e-18), label = "log10 ||C(q)||", color = "#444444")];
            title = "Constraint stability",
            xlabel = "Time (s)",
            ylabel = "log10 ||C(q)||",
        )

        write_line_plot_svg(
            joinpath(output_dir, "energy.svg"),
            result.t,
            [(values = result.energy, label = "E(t)", color = "#C77700")];
            title = "Energy history",
            xlabel = "Time (s)",
            ylabel = "Energy (J)",
        )

        write_animation_html(joinpath(output_dir, "system_motion.html"), result, params)
        return output_dir
    end
end

# ╔═╡ 7e6ddfa3-2d84-4c34-81bd-f0ff0e561004
md"""
## Model Description

The generalized coordinate vector is defined as

```math
\mathbf{q} = [x_1,\ y_1,\ \theta_1,\ x_2,\ y_2,\ \theta_2]^T
```

The holonomic constraint equations are

```math
\mathbf{C}(\mathbf{q}) =
\begin{bmatrix}
y_1 \\
\theta_1 \\
x_2 - x_1 - \frac{L}{2}\cos\theta_2 \\
y_2 - y_1 - \frac{L}{2}\sin\theta_2
\end{bmatrix}
= \mathbf{0}
```

These four constraints reduce the system to two physical degrees of freedom: the horizontal displacement of the slider and the rotation angle of the bar.
"""

# ╔═╡ 7e6ddfa3-2d84-4c34-81bd-f0ff0e561005
NOTEBOOK_DIR = @__DIR__

# ╔═╡ 7e6ddfa3-2d84-4c34-81bd-f0ff0e561006
params = default_parameters(dt = 1e-3, t_end = 8.0)

# ╔═╡ 7e6ddfa3-2d84-4c34-81bd-f0ff0e561007
simulation_setup = (
    x1 = 0.08,
    theta2 = deg2rad(55.0),
    x1dot = 0.0,
    theta2dot = 0.0,
)

# ╔═╡ 7e6ddfa3-2d84-4c34-81bd-f0ff0e561008
result = simulate_system(params; simulation_setup...)

# ╔═╡ 7e6ddfa3-2d84-4c34-81bd-f0ff0e561009
constraint_forces = generalized_constraint_forces(result, params)

# ╔═╡ 7e6ddfa3-2d84-4c34-81bd-f0ff0e561010
output_dir = generate_outputs(
    NOTEBOOK_DIR;
    params = params,
    result = result,
    constraint_forces = constraint_forces,
)

# ╔═╡ 7e6ddfa3-2d84-4c34-81bd-f0ff0e561011
md"""
## Augmented Equations of Motion

At every time step, the accelerations and Lagrange multipliers are obtained by solving

```math
\begin{bmatrix}
\mathbf{M} & \mathbf{C_q}^T \\
\mathbf{C_q} & \mathbf{0}
\end{bmatrix}
\begin{bmatrix}
\ddot{\mathbf{q}} \\
\boldsymbol{\lambda}
\end{bmatrix}
=
\begin{bmatrix}
\mathbf{Q}_e \\
\boldsymbol{\gamma} - 2\alpha \dot{\mathbf{C}} - \beta^2 \mathbf{C}
\end{bmatrix}
```

The last two terms on the right-hand side are Baumgarte stabilization terms that suppress constraint drift. The integrator uses a fixed-step RK4 scheme and projects the state back onto the constraint manifold after each step.
"""

# ╔═╡ 7e6ddfa3-2d84-4c34-81bd-f0ff0e561012
md"""
## Numerical Summary

- Maximum block displacement: **$(round(maximum(abs.(result.q[1, :])), digits = 4)) m**
- Maximum bar angle magnitude: **$(round(maximum(abs.(rad2deg.(result.q[6, :]))), digits = 2)) deg**
- Maximum constraint error norm: **$(round(maximum(result.constraint_error), sigdigits = 3))**
- Simulated duration: **$(round(last(result.t), digits = 2)) s**
- Output directory: `$(output_dir)`
"""

# ╔═╡ 7e6ddfa3-2d84-4c34-81bd-f0ff0e561013
md"## Block translation"

# ╔═╡ 7e6ddfa3-2d84-4c34-81bd-f0ff0e561014
load_svg(joinpath(output_dir, "block_position.svg"))

# ╔═╡ 7e6ddfa3-2d84-4c34-81bd-f0ff0e561015
md"## Bar rotation"

# ╔═╡ 7e6ddfa3-2d84-4c34-81bd-f0ff0e561016
load_svg(joinpath(output_dir, "bar_angle.svg"))

# ╔═╡ 7e6ddfa3-2d84-4c34-81bd-f0ff0e561017
md"## Constraint forces"

# ╔═╡ 7e6ddfa3-2d84-4c34-81bd-f0ff0e561018
load_svg(joinpath(output_dir, "constraint_forces.svg"))

# ╔═╡ 7e6ddfa3-2d84-4c34-81bd-f0ff0e561019
md"## Constraint stability"

# ╔═╡ 7e6ddfa3-2d84-4c34-81bd-f0ff0e561020
load_svg(joinpath(output_dir, "constraint_error.svg"))

# ╔═╡ 7e6ddfa3-2d84-4c34-81bd-f0ff0e561021
md"## Energy history"

# ╔═╡ 7e6ddfa3-2d84-4c34-81bd-f0ff0e561022
load_svg(joinpath(output_dir, "energy.svg"))

# ╔═╡ 7e6ddfa3-2d84-4c34-81bd-f0ff0e561024
load_html(joinpath(output_dir, "system_motion.html"))


# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This notebook uses only Julia standard libraries.
"""

# ╔═╡ Cell order:
# ╠═7e6ddfa3-2d84-4c34-81bd-f0ff0e561001
# ╠═7e6ddfa3-2d84-4c34-81bd-f0ff0e561003
# ╟─7e6ddfa3-2d84-4c34-81bd-f0ff0e561004
# ╟─7e6ddfa3-2d84-4c34-81bd-f0ff0e561011
# ╠═7e6ddfa3-2d84-4c34-81bd-f0ff0e561005
# ╠═7e6ddfa3-2d84-4c34-81bd-f0ff0e561006
# ╠═7e6ddfa3-2d84-4c34-81bd-f0ff0e561007
# ╠═7e6ddfa3-2d84-4c34-81bd-f0ff0e561008
# ╠═7e6ddfa3-2d84-4c34-81bd-f0ff0e561009
# ╠═7e6ddfa3-2d84-4c34-81bd-f0ff0e561010
# ╟─7e6ddfa3-2d84-4c34-81bd-f0ff0e561012
# ╟─7e6ddfa3-2d84-4c34-81bd-f0ff0e561013
# ╠═7e6ddfa3-2d84-4c34-81bd-f0ff0e561014
# ╟─7e6ddfa3-2d84-4c34-81bd-f0ff0e561015
# ╠═7e6ddfa3-2d84-4c34-81bd-f0ff0e561016
# ╟─7e6ddfa3-2d84-4c34-81bd-f0ff0e561017
# ╠═7e6ddfa3-2d84-4c34-81bd-f0ff0e561018
# ╟─7e6ddfa3-2d84-4c34-81bd-f0ff0e561019
# ╠═7e6ddfa3-2d84-4c34-81bd-f0ff0e561020
# ╟─7e6ddfa3-2d84-4c34-81bd-f0ff0e561021
# ╠═7e6ddfa3-2d84-4c34-81bd-f0ff0e561022
# ╠═7e6ddfa3-2d84-4c34-81bd-f0ff0e561024
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
