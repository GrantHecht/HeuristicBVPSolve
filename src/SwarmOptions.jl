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

    # Max total time for MSPSO solve
    maxTotalTime::T

    # Objective limit
    objLimit::T

    # Use multithreading
    useParallel::Bool

    # Callback function 
    callback::CF

    # Reset distance
    resetDistance::T

    # Maximum number of resets
    maxResets::Int

    # File output flag
    fileOutput::Bool

    # Solution output file
    solOutFile::String
end

function SwarmOptions(;display=false, displayInterval=1, funcTol::T=1e6,
    funValCheck=true, iUB::Uu=nothing, iLB::Ul=nothing, maxIters=1000,
    maxStallIters=25, maxStallTime::T=500.0, maxTime::T=1800.0, objLimit::T=-Inf, 
    useParallel=false, callback::CF=nothing, resetDistance::T=1.0, maxResets=25,
    solOutFile="NoFileOutput", maxTotalTime::T=3600.0) where 
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

    # Check if outputing to file
    fileOutput = true
    if solOutFile == "NoFileOutput"
        fileOutput = false
    end
    if fileOutput == true && MPI.Comm_rank(MPI.COMM_WORLD) == 0
        if isfile(solOutFile)
            @warn "Solution output files already exists and will be overwritten!"
            rm(solOutFile)
        end
        touch(solOutFile)
    end

    # If multithreading is on, check that we have acces to multiple threads
    if useParallel == true
        if Threads.nthreads == 1
            @warn "Multithreading was turned on but process only has access to 1 thread. Turning multithreading off."
            useParallel = false
        end
    end
        
    return SwarmOptions{T,U,CF}(display, displayInterval, funcTol,
        funValCheck, iUB, iLB, maxIters, maxStallIters, maxStallTime,
        maxTime, maxTotalTime, objLimit, useParallel, callback, 
        resetDistance, maxResets, fileOutput, solOutFile)
end
