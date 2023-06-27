using Chemfiles
using ProgressBars
using CUDA

function calculate_distance_matrix!(distance_matrix::CuMatrix{Float64}, frame::Frame)::CuMatrix{Float64}
    n_atoms = size(distance_matrix, 1)

    for i in 1:n_atoms
        for j in i+1:n_atoms
            dist = distance(frame, i - 1, j - 1)::Float64
            distance_matrix[i, j] = dist
            distance_matrix[j, i] = dist
        end
    end

    return distance_matrix
end

function distance_cutoff(distance::Float64, rcutoff::Float64=6.0)::Float64
    if distance > rcutoff
        return 0.0
    else
        return 0.5 * (cos(π * distance / rcutoff) + 1.0)
    end
end

function compute_g2_element(distance::Float64, eta::Float64, rcutoff::Float64, rshift::Float64)::Float64
    if distance > 0.0
        return exp(-eta * (distance - rshift)^2) * distance_cutoff(distance, rcutoff)
    else
        return 0.0
    end
end

function compute_g2(distances::CuVector{Float64}, eta::Float64, rcutoff::Float64, rshift::Float64)::Float64
    return sum(compute_g2_element(d, eta, rcutoff, rshift) for d in CuArray(distances))
end

function build_g2_matrix(distance_matrix::CuMatrix{Float64}, eta, rcutoff, rshift)
    N = size(distance_matrix, 1)
    g2_matrix = CUDA.zeros(Float64, N)

    for i in 1:N
        distance_vector = distance_matrix[i, :]
        g2_matrix[i] = compute_g2(distance_vector, eta, rcutoff, rshift)
    end

    return g2_matrix
end

function compute_mean_g2(trajectory, eta, rcutoff, rshift, n_atoms)
    distance_matrix = CUDA.zeros(Float64, n_atoms, n_atoms)

    mean_g2_value = 0.0
    for frame in trajectory
        calculate_distance_matrix!(distance_matrix, frame)
        g2_matrix = build_g2_matrix(distance_matrix, eta, rcutoff, rshift)
        mean_g2_value += sum(g2_matrix) / length(g2_matrix)
    end

    return mean_g2_value
end

function main()
    CUDA.allowscalar(false)

    filename = "100CH3OH-CG-200.xtc"
    trajectory = Trajectory(filename)
    n_atoms = size(read(trajectory))::Int

    rcutoff = 10.0
    eta_range = collect(0.0:1.0:5.0)
    rshift_range = collect(0.0:1.0:5.0)
    n_eta = length(eta_range)
    n_rshift = length(rshift_range)
    g2_values_matrix = CUDA.zeros(Float64, n_eta, n_rshift)

    println("Iterations: $(length(eta_range) * length(rshift_range))")

    for (i, eta) in tqdm(enumerate(eta_range))
        for (j, rshift) in enumerate(rshift_range)
            trajectory = Trajectory(filename)
            g2_values_matrix[i, j] = compute_mean_g2(trajectory, eta, rcutoff, rshift, n_atoms)
        end
    end

    output_file = "julia-g2_values.txt"
    f = open(output_file, "w")

    for row in eachrow(g2_values_matrix)
        line = join(CUDA.readd(row), "    ")
        write(f, line * "\n")
    end

    close(f)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
