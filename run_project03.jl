using Printf

if !isdefined(@__MODULE__, :Project03Core)
    include(joinpath(@__DIR__, "project03_core.jl"))
end

using .Project03Core

default_output_dir(base_dir::AbstractString = @__DIR__) = joinpath(base_dir, "outputs")

function write_csv(path::AbstractString, result::SimulationResult, constraint_forces::AbstractMatrix)
    headers = [
        "t", "x1", "y1", "theta1", "x2", "y2", "theta2",
        "vx1", "vy1", "omega1", "vx2", "vy2", "omega2",
        "ax1", "ay1", "alpha1", "ax2", "ay2", "alpha2",
        "lambda1", "lambda2", "lambda3", "lambda4",
        "constraint_error", "energy",
        "Qc1", "Qc2", "Qc3", "Qc4", "Qc5", "Qc6",
    ]
    open(path, "w") do io
        println(io, join(headers, ","))
        for i in eachindex(result.t)
            row = vcat(
                result.t[i],
                result.q[:, i],
                result.v[:, i],
                result.a[:, i],
                result.lambda[:, i],
                result.constraint_error[i],
                result.energy[i],
                constraint_forces[:, i],
            )
            println(io, join((@sprintf("%.10f", value) for value in row), ","))
        end
    end
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
    body { font-family: Arial, sans-serif; margin: 24px; background: #fafafa; color: #222; }
    h1 { font-family: "Times New Roman", serif; }
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
  <h1>Spring-block compound pendulum animation</h1>
  <canvas id="canvas" width="920" height="560"></canvas>
  <div class="controls">
    <button id="toggle">Pause</button>
    <label>Frame <input id="frame" type="range" min="0" max="$(length(times) - 1)" value="0" /></label>
    <span id="stamp"></span>
  </div>
  <div class="note">
    The animation is exported directly from the constrained dynamics solution. The slider moves horizontally, the spring compresses and stretches, and the rigid bar rotates about the pin on the block.
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

function write_summary_readme(base_dir::AbstractString = @__DIR__)
    summary_text = """
# Project 03 Solution

This folder contains a complete Julia solution for the spring-block compound pendulum assignment.

## Files

- `project03_core.jl`: constrained multibody model and solver
- `run_project03.jl`: batch script that runs the simulation and exports plots
- `project_03_pluto.jl`: Pluto notebook version with derivation and result links
- `outputs/`: generated figures, browser animation, and CSV time history

## Modeling assumptions

1. The block reference point is taken at the pin, so the slider constraints are `y1 = 0` and `theta1 = 0`.
2. The bar is modeled as a uniform rigid rod of length `L = 0.4 m` with center of mass at `L / 2`.
3. Since the assignment figure does not specify block dimensions, the block inertia is estimated from a `0.10 m x 0.06 m` rectangle. This affects only the constrained rotational reaction.
4. Constraint stabilization uses Baumgarte terms plus a projection back to the constraint manifold after each RK4 step.

## Constraints used

`C(q) = [y1, theta1, x2 - x1 - L/2 cos(theta2), y2 - y1 - L/2 sin(theta2)]^T = 0`

## Run

```julia
julia run_project03.jl
```

Open `outputs/system_motion.html` in a browser to watch the animation.
"""

    write(joinpath(base_dir, "README.md"), summary_text)
end

function generate_outputs(
    base_dir::AbstractString = @__DIR__;
    params::ModelParameters = default_parameters(dt = 1e-3, t_end = 8.0),
    result::SimulationResult = simulate_system(params; x1 = 0.08, theta2 = deg2rad(55.0), x1dot = 0.0, theta2dot = 0.0),
    constraint_forces::AbstractMatrix = generalized_constraint_forces(result, params),
)
    output_dir = default_output_dir(base_dir)
    mkpath(output_dir)

    write_csv(joinpath(output_dir, "simulation_history.csv"), result, constraint_forces)

    write_line_plot_svg(
        joinpath(output_dir, "block_position.svg"),
        result.t,
        [
            (values = result.q[1, :], label = "x1 (m)", color = "#2458A6"),
        ];
        title = "Block translation",
        xlabel = "Time (s)",
        ylabel = "x1 (m)",
    )

    write_line_plot_svg(
        joinpath(output_dir, "bar_angle.svg"),
        result.t,
        [
            (values = rad2deg.(result.q[6, :]), label = "theta2 (deg)", color = "#B84732"),
        ];
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
        [
            (values = log10.(result.constraint_error .+ 1e-18), label = "log10 ||C(q)||", color = "#444444"),
        ];
        title = "Constraint stability",
        xlabel = "Time (s)",
        ylabel = "log10 ||C(q)||",
    )

    write_line_plot_svg(
        joinpath(output_dir, "energy.svg"),
        result.t,
        [
            (values = result.energy, label = "E(t)", color = "#C77700"),
        ];
        title = "Energy history",
        xlabel = "Time (s)",
        ylabel = "Energy (J)",
    )

    write_animation_html(joinpath(output_dir, "system_motion.html"), result, params)
    return output_dir
end

function run_batch_export(base_dir::AbstractString = @__DIR__)
    params = default_parameters(dt = 1e-3, t_end = 8.0)
    result = simulate_system(params; x1 = 0.08, theta2 = deg2rad(55.0), x1dot = 0.0, theta2dot = 0.0)
    constraint_forces = generalized_constraint_forces(result, params)
    write_summary_readme(base_dir)
    return generate_outputs(base_dir; params = params, result = result, constraint_forces = constraint_forces)
end

if abspath(PROGRAM_FILE) == @__FILE__
    output_dir = run_batch_export()
    println("Outputs written to: ", output_dir)
end
