
import MPI
import Quadtrees

include("check_duplicates.jl")

function dup_mpi(comm,x,delta)
    
rank = MPI.Comm_rank(comm)
world_size = MPI.Comm_size(comm)


if rank == 0
    Nobs = length(x[1])

    n = length(x)

    x1min = minimum(x[1])
    x1max = maximum(x[1])
    label = collect(1:Nobs)

    for r = 0:world_size-1


        x1min_loc = x1min +    r  * (x1max-x1min)/world_size 
        x1max_loc = x1min + (r+1) * (x1max-x1min)/world_size

        # make domain overlap
        
        x1min_loc -= delta[1]
        x1max_loc += delta[1]
        
        sel = x1min_loc .<= x[1] .<= x1max_loc
        xloc = ntuple(i -> x[i][sel],n)
        labelloc = label[sel]

        #@show "send to",r, r+30
        MPI.send(xloc, r, r+30, comm)
        #@show "send to",r, r+40
        MPI.send(labelloc, r, r+40, comm)
    end
end

#@show "get xloc",rank, rank+30
xloc,status = MPI.recv(0,  rank+30, comm)
#@show "get label",rank, rank+40
labelloc,status = MPI.recv(0,  rank+40, comm)


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
    
    return dup
end

end

#function main()
MPI.Init()

comm = MPI.COMM_WORLD
delta = [0.2,0.2]
x = ([1,2,4,1,2.46,2.54],[2,3,5,2,1,1])

dup_mpi(comm,x,delta)

MPI.Finalize()

#end

#main()
