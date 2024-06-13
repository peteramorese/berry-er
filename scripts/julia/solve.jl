#### Created with the help of ChatGPT ####

using CSV
using DataFrames
using JuMP
using GLPK
using MosekTools

# Function to read matrix A and vector b from CSV files
function read_matrix(file::String)
    open(file, "r") do file
        rows, cols = split(readline(file)) |> x -> (parse(Int, x[1]), parse(Int, x[2]))
        mat = Array{Float64}(undef, rows, cols)
        for i in 1:rows
            mat[i, :] = split(readline(file)) |> x -> parse.(Float64, x)
        end
        return mat
    end
end

# Function to solve Ax >= b using JuMP and GLPK solver
function solve_linear_program(A::Matrix{Float64}, b::Vector{Float64}, N::Int)
    rows, cols = size(A)
    #model = Model(GLPK.Optimizer)
    model = Model(optimizer_with_attributes(Mosek.Optimizer))

    tolerance = 1e-4
    set_optimizer_attribute(model, "MSK_DPAR_INTPNT_TOL_PFEAS", tolerance)  # Primal feasibility tolerance
    set_optimizer_attribute(model, "MSK_DPAR_INTPNT_TOL_DFEAS", tolerance)  # Dual feasibility tolerance
    set_optimizer_attribute(model, "MSK_DPAR_INTPNT_TOL_INFEAS", tolerance)  # Infeasibility
    set_optimizer_attribute(model, "MSK_DPAR_INTPNT_TOL_REL_GAP", tolerance)  # Relative gap tolerance

    @variable(model, x[1:cols] >= 0)
    for i in 1:rows
        @constraint(model, sum(A[i,j] * x[j] for j in 1:cols) >= b[i])
    end
    objective_vec = Vector{Float64}(undef, cols)
    objective_vec[end] = N
    objective_vec[end - 1] = 1.0
    @objective(model, Min, sum(objective_vec[j] * x[j] for j in 1:cols)) 
    optimize!(model)
    status = termination_status(model) 
    if status == MOI.OPTIMAL
        return value.(x)
    elseif status == MOI.INFEASIBLE_OR_UNBOUNDED
        println("The problem is infeasible or unbounded")
    else
        println("The problem terminated with status ", status)
        return value.(x)
    end
end

# Main function to run the script
function main()
    # Replace "A.csv" and "b.csv" with your actual file names
    N = 10
    file_directory = "../../build/bin/";
    A = read_matrix(file_directory * "A.txt")
    b = vec(read_matrix(file_directory * "b.txt"))

    gamma_cond = zeros(1, size(A, 2))
    gamma_cond[end] = 1
    eta_cond = zeros(1, size(A, 2))
    eta_cond[end - 1] = 1
    A = vcat(A, gamma_cond)
    A = vcat(A, eta_cond)
    push!(b, 0)
    push!(b, 0)

    #println("A: ", A)

    vars = solve_linear_program(A, b, N)
    coeffs = vars[1:end - 2]
    η = vars[end - 1];
    γ = vars[end];
    println("Eta: ", η)
    println("Gamma: ", γ)
    println("Psafe: ", 1 - (η + N * γ))
    println("Coeffs:")
    for c in coeffs
        print(c, " ")
    end
    print("\n")
end

main()
