struct SwarmOptions{T<:AbstractFloat, U<:AbstractVector, CF<:Union{Function,Nothing}}
    # Display options
    display::Bool
    displayInterval::Int

    # Function tolerance
    funcTol::T

    # Check function val for NaN and Inf
    funValCheck::Bool

    # Initial bounds
    iUB::U
    iLB::U

    # Max iterations 
    maxIters::Int

    # Max stall Iterations
    maxStallIters::Int

    # Max stall time
    maxStallTime::T

    # Max time 
    maxTime::T

    # Objective limit
    objLimit::T

    # Use multithreading
    useParallel::Bool

    # Callback function 
    callback::CF
end

function SwarmOptions(;display=false, displayInterval=1, funcTol::T=1e6,
    funValCheck=true, iUB::Uu=nothing, iLB::Ul=nothing, maxIters=1000,
    maxStallIters=25, maxStallTime=500, maxTime=1800, objLimit=-Inf, 
    useParallel=false, callback::CF=nothing) where 
    {T<:AbstractFloat, Uu<:Union{Nothing,Vector}, Ul<:Union{Nothing,Vector},
     CF<:Union{Nothing, Function}}

    U = Vector{T} 
    if iUB === nothing
        iUB = U([])
    else
        iUB = U(iUB)
    end
    if iLB === nothing
        iLB = U([])
    else
        iLB = U(iLB)
    end
        
    return Options{T,U,CF}(display, displayInterval, funcTOl,
        funValCheck, iUB, iLB, maxIters, maxStallIters, masStallTime,
        maxTime, objLimit, useParallel, callback)
end