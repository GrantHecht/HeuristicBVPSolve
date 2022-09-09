mutable struct Swarm{T<:AbstractFloat}
    # Vector of particles
    particles::Vector{Particle{T}}

    # Preallocated vector of integers for neighborhood selection
    nVec::Vector{Int}

    # Global best objective function value
    b::T

    # Location of global best objective function value
    d::Vector{T}

    # Neighborhood size
    n::Int

    # Inertia
    w::T

    # Adaptive inertia counter
    c::Int

    # Adjustment weights
    y₁::T
    y₂::T
end

# UndefInitializer constructor
function Swarm{T}(::UndefInitializer) where {T}
    return Swarm{T}(Vector{Particle{T}}(undef, 0),
        Vector{Int}(undef, 0), T(0.0), Vector{T}(undef, 0),
        Int(0), T(0.0), Int(0), T(0.0), T(0.0))
end

# Standard constructor
function Swarm{T}(nDims::Integer, nParticles::Integer) where {T}
    if nDims < 0
        throw(ArgumentError("nDims must be greater than 0."))
    end
    if nParticles < 0
        throw(ArgumentError("nParticles must be greater than 0."))
    end
    return Swarm{T}([Particle{T}(nDims) for i in 1:nParticles],
        Vector{Int}(1:nParticles), T(0.0), Vector{T}(undef, nDims),
        Int(0), T(0.0), Int(0), T(0.0), T(0.0))
end

# Methods
Base.length(s::Swarm) = length(s.particles)
Base.eachindex(s::Swarm) = eachindex(s.particles)

function Base.getindex(s::Swarm, i::Int)
    1 <= i <= length(s) || throw(BoundsError(s, i))
    return s.particles[i]
end

function Base.setindex!(s::Swarm, v, i::Int)
    1 <= i <= length(s) || throw(BoundsError(s, i))
    s.particles[i] = v
    return nothing
end

function feval!(s::Swarm, f::Function, opts::SwarmOptions; init = false)
    # Evaluate objective function
    @inbounds begin
        if opts.useParallel
            Threads.@threads for i in eachindex(s)
                s[i].fx = f(s[i].x)
            end
        else
            for i in eachindex(s)
                s[i].fx = f(s[i].x)
            end
        end
    end

    # Check bjective function values if desired
    if opts.funValCheck 
        fValCheck(s)
    end

    # Update each particles best objective function value and it's location
    if !init
        @inbounds for i in eachindex(s)
            if s[i].fx < s[i].fp
                s[i].p .= s[i].x
                s[i].fp = s[i].fx
            end
        end
    else
        @inbounds begin
            @inbounds for i in eachindex(s)
                s[i].fp = s[i].fx
            end
        end
    end
    return nothing
end

function fValCheck(s::Swarm)
    @inbounds begin
        for i in eachindex(s)
            if isinf(s[i].fx) || isnan(s[i].fx)
                throw(ErrorException("Objective function return Inf or NaN!"))
            end
        end
    end
end

function setGlobalBest!(s::Swarm)
    updated = false
    @inbounds begin
        for i in eachindex(s)
            if s[i].fp < s.b 
                s.b = s[i].fp
                s.d .= s[i].p 
                updated = true
            end
        end
    end
    return updated
end

function updateVelocities!(s::Swarm)
    n = length(s.d)
    @inbounds begin
        for i in eachindex(s)
            # Shuffle vector containing integers 1:n
            # first m != i will be neighbors
            shuffle!(s.nVec)

            # Determine fbest(S)
            fbest = Inf
            best = 0
            incr = 0
            for j in 1:s.n 
                s.nVec[j] == i ? incr = 1 : ()
                k = s.nVec[j + incr]
                if s[k].fp < fbest
                    fbest = s[k].fp
                    best = k
                end
            end

            # Update i's velocity
            for j in 1:n
                s[i].v[j] = s.w*s[i].v[j] + s.y₁*rand()*(s[i].p[j] - s[i].x[j]) +
                    s.y₂*rand()*(s[best].p[j] - s[i].x[j])
            end
        end
    end
    return nothing
end

function step!(s::Swarm)
    @inbounds for i in eachindex(s)
        s[i].x .= s[i].x .+ s[i].v
    end
end

function enforceBounds!(s::Swarm, LB, UB)
    @inbounds begin 
        for i in eachindex(s)
            for j in 1:length(s.d)
                if s[i].x[j] > UB[j] 
                    s[i].x[j] = UB[j]
                    s[i].v[j] = 0.0
                elseif s[i].x[j] < LB[j]
                    s[i].x[j] = LB[j] 
                    s[i].v[j] = 0.0
                end
            end
        end
    end
    return nothing
end

# Initializes position and velocities of particles by sampling from 
# a uniform distribution
function uniformInitialization!(swarm::Swarm, prob::Problem, opts::SwarmOptions)
    # Get N: Number of diamensions
    N = length(swarm[1])

    # Get Boundary Constraints
    LB  = prob.LB
    UB  = prob.UB

    # Check if initial bounds on positions have been set
    useInitBnds = false
    if length(opts.iLB) == N && length(opts.iUB) == N
        useInitBnds = true
        iLB = opts.iLB            
        iUB = opts.iUB
    end

    # Initialize particle positions and velocities
    @inbounds begin
        for d in 1:N
            # Get local bounds for d-axis
            lLB = useInitBnds ? (LB[d] < iLB[d] ? iLB[d] : LB[d]) : LB[d]
            lUB = useInitBnds ? (UB[d] > iUB[d] ? iUB[d] : UB[d]) : UB[d]
            for p in eachindex(swarm)
                # Position information
                swarm[p].x[d] = lLB + (lUB - lLB)*rand()
                swarm[p].p[d] = swarm[p].x[d]

                # Velocity 
                r = useInitBnds ? min(lUB-lLB,UB[d]-LB[d]) : lUB - lLB 
                swarm[p].v[d] = -r + 2*r*rand()
            end
        end
    end
    
    return nothing
end