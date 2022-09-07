mutable struct Particle{T<:AbstractFloat}
    x::Vector{T}    # Particle position
    v::Vector{T}    # Particle velocity
    p::Vector{T}    # Particle best position

    fx::T           # Particle current function value
    fp::T           # Particle best function value

    # Inner constructor
    function Particle{T}(nDims::Integer) where {T}
        if nDims < 0
            throw(ArgumentError("nDims cannot be less than 0."))
        end
        return new{T}(Vector{T}(undef, nDims),
                      Vector{T}(undef, nDims),
                      Vector{T}(undef, nDims),
                      T(0.0), T(0.0))
    end
end

# UndefInitializer Particle constructor
function Particle{T}(::UndefInitializer) where {T}
    return Particle{T}(0)
end

# Particle methods
Base.length(p::Particle) = length(p.x)