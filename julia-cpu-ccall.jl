# Julia code with C function calls

# Load the necessary packages
using Chemfiles
using ProgressBars

# Function to calculate distance matrix
function calculate_distance_matrix!(distance_matrix::Matrix{Float64}, frame::Frame)::Matrix{Float64}
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

# Function to compute g2 matrix element
function compute_g2(distances::Vector{Float64}, eta::Float64, rcutoff::Float64, rshift::Float64)::Float64
    result = 0.0
    for d in distances
        result += @ccall "./g2.so".compute_g2_element(d::Cdouble, eta::Cdouble, rcutoff::Cdouble, rshift::Cdouble)::Cdouble
        # result += compute_g2_element(d, eta, rcutoff, rshift)
    end
    return result
end


# Function to build g2 matrix
function build_g2_matrix(distance_matrix::Matrix{Float64}, eta, rcutoff, rshift)
    N = size(distance_matrix, 1)
    g2_matrix = zeros(Float64, N)

    for i in 1:N
        distance_vector = distance_matrix[i, :]
        g2_matrix[i] = compute_g2(distance_vector, eta, rcutoff, rshift)
    end

    return g2_matrix
end

# Function to compute mean g2
function compute_mean_g2(trajectory, eta, rcutoff, rshift, n_atoms)
    distance_matrix = zeros(Float64, n_atoms, n_atoms)

    mean_g2_value = 0.0
    for frame in trajectory
        calculate_distance_matrix!(distance_matrix, frame)
        g2_matrix = build_g2_matrix(distance_matrix, eta, rcutoff, rshift)
        mean_g2_value += sum(g2_matrix) / length(g2_matrix)
    end

    return mean_g2_value
end

# Main function
function main()
    filename = "100CH3OH-CG-200.xtc"
    trajectory = Trajectory(filename)
    n_atoms = size(read(trajectory))::Int

    rcutoff = 10.0
    eta_range = collect(0.0:1.0:10.0)
    rshift_range = collect(0.0:1.0:10.0)
    n_eta = length(eta_range)
    n_rshift = length(rshift_range)
    g2_values_matrix = zeros(Float64, n_eta, n_rshift)

    println("Iterations: $(length(eta_range) * length(rshift_range))")

    for (i, eta) in tqdm(enumerate(eta_range))
        for (j, rshift) in enumerate(rshift_range)
            trajectory = Trajectory(filename)
            g2_values_matrix[i, j] = compute_mean_g2(trajectory, eta, rcutoff, rshift, n_atoms)
        end
    end

    output_file = "julia-ccall-g2_values.txt"
    f = open(output_file, "w")

    for row in eachrow(g2_values_matrix)
        line = join(row, "    ")
        write(f, line * "\n")
    end

    close(f)
end

if abspath(PROGRAM_FILE) == @__FILE__
    main()
end
