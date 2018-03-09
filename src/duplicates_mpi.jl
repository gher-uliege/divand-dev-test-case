
import MPI
import Quadtrees

include("check_duplicates.jl")

function main()
    MPI.Init()

    comm = MPI.COMM_WORLD

    rank = MPI.Comm_rank(comm)
    world_size = MPI.Comm_size(comm)

    x = ([1,2,4,1,2.46,2.54],[2,3,5,2,1,1])
    Nobs = length(x[1])

    delta = [0.2,0.2]


    n = length(x)
    #dup = check_duplicates(x,delta; maxcap = 100)

    x1min = minimum(x[1])
    x1max = maximum(x[1])



    label = collect(1:Nobs)

    x1min_loc = x1min +    rank  * (x1max-x1min)/world_size 
    x1max_loc = x1min + (rank+1) * (x1max-x1min)/world_size

    # make domain overlap

    x1min_loc -= delta[1]
    x1max_loc += delta[1]

    sel = x1min_loc .<= x[1] .<= x1max_loc
    xloc = ntuple(i -> x[i][sel],n)
    labelloc = label[sel]

    @show xloc

    duploc = check_duplicates(
        xloc,delta;
        maxcap = 100,
        label = labelloc)

    MPI.send(duploc, 0, rank+32, comm)

    @show duploc


    if rank == 0        
        duplocall = Vector{Vector{Set{Int}}}(world_size)
        for i = 1:world_size
            duplocall[i],status = MPI.recv(i-1,  i+32-1, comm)
        end

        @show duplocall
        
        # merge all duploc
        
        dup = dupset(cat(1,duploc...))

        @show dup
    end

    MPI.Finalize()

end

main()
