
using MPI
using HeuristicNLSolve

function main()
    MPI.Init() 

    # Define function
    function psoCost(x)
        c1 = x[2] - x[1] - 1
        c2 = x[1]^2 - x[2] + 1
        return c1^2 + c2^2
    end

    # Create problem
    prob = Problem(psoCost, [-50, -50], [50, 50])

    # Set swarm options
    options = SwarmOptions(;display=true, resetDistance=0.1, maxResets=100000,
        maxTotalTime    = 10.0,
        solOutFile      = "/Users/granthec/.julia/dev/HeuristicNLSolve/scripts/test.txt")

    # Instantiate DMSPSO
    mspso   = DMSPSO(prob; numParticlesPerSwarm = 50, commIterationBuffer = 5) 

    # Perform optimization
    optimize!(mspso, options)

    # Finalize MPI
    MPI.Finalize()
end

main()
