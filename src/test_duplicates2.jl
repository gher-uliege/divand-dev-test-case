

x = ([1,2,4,1,2.46,2.54],[2,3,5,2,1,1])
Nobs = length(x[1])

delta = [0.2,0.2]


n = length(x)
#dup = check_duplicates(x,delta; maxcap = 100)

x1min = minimum(x[1])
x1max = maximum(x[1])

rank = 0
world_size = 2

duploc = Vector{Vector{Set{Int}}}(world_size)
label = collect(1:Nobs)

for rank = 0:world_size-1
    x1min_loc = x1min +    rank  * (x1max-x1min)/world_size 
    x1max_loc = x1min + (rank+1) * (x1max-x1min)/world_size
    
    # make domain overlap
    
    x1min_loc -= delta[1]
    x1max_loc += delta[1]
    
    sel = x1min_loc .<= x[1] .<= x1max_loc
    xloc = ntuple(i -> x[i][sel],n)
    labelloc = label[sel]

    @show xloc
    
    duploc[rank+1] = check_duplicates(
        xloc,delta;
        maxcap = 100,
        label = labelloc)

    @show duploc[rank+1]
end

@show duploc

# merge all duploc

dup = dupset(cat(1,duploc...))

