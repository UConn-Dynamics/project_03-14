### A Pluto.jl notebook ###
# v0.20.21

using Markdown
using InteractiveUtils

# ╔═╡ 9ba23e1e-cf31-4697-a343-f1ec2d60b2e2
begin
    function locate_project03_dir()
        candidates = unique(filter(isdir, [
            @__DIR__,
            pwd(),
            dirname(@__DIR__),
            dirname(pwd()),
        ]))
        for dir in candidates
            if isfile(joinpath(dir, "project03_core.jl")) && isfile(joinpath(dir, "run_project03.jl"))
                return dir
            end
        end
        error("Could not locate project03_core.jl and run_project03.jl. Please open the delivered notebook file directly from its folder.")
    end

    const PROJECT03_DIR = locate_project03_dir()
    include(joinpath(PROJECT03_DIR, "run_project03.jl"))
    using .Project03Core
end

# ╔═╡ a16327fc-f552-42eb-8d3f-d09e63662d17
md"""
# Project_03 - Multibody Dynamic Modeling

This notebook solves the spring-block plus compound pendulum system with an augmented constrained dynamics formulation and links to the exported results in `outputs/`.

The generalized coordinate vector is

```math
\mathbf{q} = [x_1,\ y_1,\ \theta_1,\ x_2,\ y_2,\ \theta_2]^T
```

and the holonomic constraints are

```math
\mathbf{C}(\mathbf{q}) =
\begin{bmatrix}
y_1 \\
\theta_1 \\
x_2 - x_1 - \frac{L}{2}\cos\theta_2 \\
y_2 - y_1 - \frac{L}{2}\sin\theta_2
\end{bmatrix}
= \mathbf{0}.
```

That leaves two physical degrees of freedom: slider translation and bar rotation.
"""

# ╔═╡ 3139ee1f-f4bf-48b5-b104-868633fdf347
params = default_parameters(dt = 1e-3, t_end = 8.0)

# ╔═╡ 0c7b4733-b9b1-417a-8090-c88e1a57d506
result = simulate_system(params; x1 = 0.08, theta2 = deg2rad(55.0), x1dot = 0.0, theta2dot = 0.0)

# ╔═╡ 1dca6b35-47c8-441c-b478-bca52aab6804
constraint_forces = generalized_constraint_forces(result, params)

# ╔═╡ 2e531bc4-6716-42fd-9efd-2dd6082448ea
output_dir = begin
    write_summary_readme(PROJECT03_DIR)
    generate_outputs(PROJECT03_DIR; params = params, result = result, constraint_forces = constraint_forces)
end

# ╔═╡ 79e88d93-7687-4d6b-b6db-c1acb8179d86
md"""
## Augmented equations

At every time step, the acceleration and Lagrange multipliers are computed from

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
\end{bmatrix},
```

where the last two terms are Baumgarte stabilization terms. The notebook uses a fixed-step RK4 integrator and projects the state back to the constraint manifold after each step.
"""

# ╔═╡ c0f14d88-bd93-4724-97f4-40e08417e915
md"""
## Quick numerical summary

- Maximum block displacement: **$(round(maximum(abs.(result.q[1, :])), digits = 4)) m**
- Maximum bar angle magnitude: **$(round(maximum(abs.(rad2deg.(result.q[6, :]))), digits = 2)) deg**
- Maximum constraint error norm: **$(round(maximum(result.constraint_error), sigdigits = 3))**
- Simulated duration: **$(round(last(result.t), digits = 2)) s**
"""

# ╔═╡ cfa3282d-af33-4d20-b23f-0d6899c7cda6
md"""
## Exported result figures

### Block translation
![](outputs/block_position.svg)

### Bar rotation
![](outputs/bar_angle.svg)

### Constraint forces
![](outputs/constraint_forces.svg)

### Constraint stability
![](outputs/constraint_error.svg)

### Energy history
![](outputs/energy.svg)
"""

# ╔═╡ c2ce4c24-2e38-4c0f-8f05-8fe7b0bdcf7d
md"""
## Notes

- The bar is treated as a uniform rigid rod, so `J_2 = m_2 L^2/12`.
- The slider block inertia is estimated from a `0.10 m x 0.06 m` rectangle because the prompt does not specify its geometry.
- The generated batch outputs are written by `run_project03.jl` into the `outputs/` folder.
- Open `outputs/system_motion.html` in a browser to view the animation.
- Output directory: `$(output_dir)`
"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This notebook uses only the Julia standard library and local files.
"""

# ╔═╡ Cell order:
# ╠═9ba23e1e-cf31-4697-a343-f1ec2d60b2e2
# ╠═a16327fc-f552-42eb-8d3f-d09e63662d17
# ╠═3139ee1f-f4bf-48b5-b104-868633fdf347
# ╠═0c7b4733-b9b1-417a-8090-c88e1a57d506
# ╠═1dca6b35-47c8-441c-b478-bca52aab6804
# ╠═2e531bc4-6716-42fd-9efd-2dd6082448ea
# ╠═79e88d93-7687-4d6b-b6db-c1acb8179d86
# ╠═c0f14d88-bd93-4724-97f4-40e08417e915
# ╠═cfa3282d-af33-4d20-b23f-0d6899c7cda6
# ╠═c2ce4c24-2e38-4c0f-8f05-8fe7b0bdcf7d
# ╠═00000000-0000-0000-0000-000000000001
# ╠═00000000-0000-0000-0000-000000000002
