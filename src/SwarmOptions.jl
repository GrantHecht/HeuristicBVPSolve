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

    # Communicate with solver proc
    solverComm::Bool

    # Solver rank
    solverRank::Int

    # Solution output file
    solOutFile::String

    # Print debug info
    printDebugInfo::Bool
end

function SwarmOptions(;display=false, displayInterval=1, funcTol::T=1e6,
    funValCheck=true, iUB::Uu=nothing, iLB::Ul=nothing, maxIters=1000,
    maxStallIters=25, maxStallTime::T=500.0, maxTime::T=1800.0, objLimit::T=-Inf, 
    useParallel=false, callback::CF=nothing, resetDistance::T=1.0, maxResets=25,
    solOutFile="NoFileOutput", solverComm=false, solverRank=-1, maxTotalTime::T=3600.0,
    printDebugInfo=false) where {T<:AbstractFloat, Uu<:Union{Nothing,Vector}, 
                                 Ul<:Union{Nothing,Vector}, CF<:Union{Nothing, Function}}

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
        close(open(solOutFile, "w"))
    end

    # If multithreading is on, check that we have acces to multiple threads
    if useParallel == true
        if Threads.nthreads == 1
            @warn "Multithreading was turned on but process only has access to 1 thread. Turning multithreading off."
            useParallel = false
        end
    end

    # If communicating with solver and solver rank unset, default to largest rank in MPI.COMM_WORLD
    if solverComm == true && solverRank == -1
        solverRank = MPI.Comm_size(MPI.COMM_WORLD) - 1
    end
        
    return SwarmOptions{T,U,CF}(display, displayInterval, funcTol,
        funValCheck, iUB, iLB, maxIters, maxStallIters, maxStallTime,
        maxTime, maxTotalTime, objLimit, useParallel, callback, 
        resetDistance, maxResets, fileOutput, solverComm, solverRank, 
        solOutFile, printDebugInfo)
end
