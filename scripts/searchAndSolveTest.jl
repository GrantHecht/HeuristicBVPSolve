using MPI
using HeuristicNLSolve

function main()
    # ====== PERFORM SETUP ON ALL PROCS
    MPI.Init() 

    # Construct communicators
    color   = (MPI.Comm_rank(MPI.COMM_WORLD) == 
               MPI.Comm_size(MPI.COMM_WORLD) - 1 ? 2 : 1)
    comm    = MPI.Comm_split(MPI.COMM_WORLD, color, 0)

    # Define function
    function nleqs!(F, J, x)
        if !(F === nothing)
            F[1] = x[2] - x[1] - 1
            F[2] = x[1]^2 - x[2] + 1
        end
        if !(J === nothing)
            J[1,1] = -1
            J[1,2] = 1
            J[2,1] = 2*x[1]
            J[2,2] = -1
        end
        return nothing
    end

    function psoCost(x)
        c1 = x[2] - x[1] - 1
        c2 = x[1]^2 - x[2] + 1
        return c1^2 + c2^2
    end

    # Create problem
    prob = Problem(psoCost, nleqs!, [-50, -50], [50, 50])

    # ===== EXECUTE THE FOLLOWING CODE DEPENDING ON COLOR
    if color == 1
        # Set swarm options
        options = SwarmOptions(;display=false, resetDistance=0.1, maxResets=10000,
            maxTotalTime    = 10.0,
            solverComm      = true,
            solOutFile      = "/Users/granthec/.julia/dev/HeuristicNLSolve/scripts/test.txt")

        # Instantiate DMSPSO
        mspso   = DMSPSO(prob; numParticlesPerSwarm = 50, commIterationBuffer = 5, comm = comm) 

        # Perform optimization
        optimize!(mspso, options)
    else
        # Set solver options
        options = SolverOptions(;fTol = 1e-8, factor = 2.0,
            solOutFile = "/Users/granthec/.julia/dev/HeuristicNLSolve/scripts/sols.txt")

        # Create solver
        solver  = SolverServer(prob)

        # Solve 
        start!(solver, options)
    end

    # Finalize MPI
    MPI.Finalize()
end

main()