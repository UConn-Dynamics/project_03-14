module Project03Core

using LinearAlgebra

export ModelParameters,
       SimulationResult,
       default_parameters,
       make_initial_state,
       constraint_vector,
       constraint_jacobian,
       generalized_forces,
       solve_augmented_system,
       simulate_system,
       generalized_constraint_forces,
       block_polygon,
       wheel_centers,
       bar_endpoints

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

function block_polygon(x::Float64, p::ModelParameters)
    w = p.block_width
    h = p.block_height
    xs = x .+ 0.5 .* [-w, w, w, -w, -w]
    ys = [0.0, 0.0, h, h, 0.0] .- 0.5 * h
    return xs, ys
end

function wheel_centers(x::Float64, p::ModelParameters)
    w = p.block_width
    r = p.wheel_radius
    y = -0.5 * p.block_height - r
    return [(x - 0.22 * w, y), (x + 0.22 * w, y)]
end

function bar_endpoints(q::AbstractVector, p::ModelParameters)
    x_pin = q[1]
    y_pin = q[2]
    theta2 = q[6]
    x_tip = x_pin + p.L * cos(theta2)
    y_tip = y_pin + p.L * sin(theta2)
    return (x_pin, y_pin), (x_tip, y_tip)
end

end
