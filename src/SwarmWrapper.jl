
# A struct to wrap a swarm which lives on a unique process.
# DMSPSO manages each local swarm wrapper intance to perform 
# distributed multi-swarm PSO
struct SwarmWrapper{T,S,fType}
    # Optimization problem
    prob::Problem{fType, S}

    # Swarm 
    swarm::Swarm{T}

    # PSO specific parameters
    inertiaRange::Tuple{T,T}
    minNeighborFrac::T
    selfAdjustWeight::T
    socialAdjustWeight::T
end