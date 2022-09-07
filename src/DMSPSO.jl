# Abstract MSPSO
abstract type MSPSO end

# DMSPSO - {D}istributed {M}ulti-{S}warm {P}article {S}warm {O}ptimization
mutable struct DMSPSO{T,S,fType} <: MSPSO
    # Process IDs with SwarmWrapper
    procIDs::Vector{Int}

    # Swarm reset counter
    resetCount::Int
end