module quadtree

using Base.Test
import Base.length
#using PyPlot

"""
quadtree (of the higher-dimensional equivalent)
T the type of the coordinates
TA the type of the attributes
N number of dimensions
"""

type QT{T,TA,N}
    children :: Vector{QT{T,TA,N}}  # vector of child nodes (empty if node is a leaf)
    points :: Array{T,2}            # list of coordinates (only non-empty if node is a leaf)
    min :: Vector{T}                # minimum of bounding box
    max :: Vector{T}                # maximim of bounding box
    attribs :: Vector{TA}           # additional attributes (only non-empty if node is a leaf)w
end

"""create empty quadtree"""
QT{T}(TA::DataType,min::Vector{T}, max::Vector{T}) = QT(QT{T,TA,size(min,1)}[],Matrix{T}(0,size(min,1)),min,max,TA[])

"""create a quadtree"""

QT{T,TA}(points::Array{T,2},attribs::Vector{TA}) = QT(QT{T,TA,size(points,2)}[],points,minimum(points,1)[:],maximum(points,1)[:],attribs)

QT{T,TA}(points::Array{T,2}, min::Vector{T}, max::Vector{T}, attribs::Vector{TA}) = QT(QT{T,TA,size(points,2)}[],points,min,max,attribs)
#QT{T,TA,N}(children::Vector{QT{T,TA,N}}, min::Vector{T}, max::Vector{T}) = QT(children,Array{T,2}(0,N),min,max,TA[])



"""       
             x1
  +----------+ 
  |          |
  |   +      |
  |   y      |
  +----------+
 x0

"""

function inside(x0,x1,y)
    insd = true

    for i = 1:length(y)
        insd = insd & (x0[i] <= y[i] <= x1[i])
    end
    return insd
end

"Test if the n-th bit in a is set. The least significant bit is n = 1."
bitget(a,n) = Bool((a & (1 << (n-1))) >> (n-1))





function intersect_(x0,x1,y0,y1)
    # number of dimensions
    n = size(x0,1)

    intrsct = false
    
    for i = 1:2^n
        inside = true
        for j = 1:n
            if bitget(i-1, j)
                inside = inside & (x0[j] <= y0[j] <= x1[j])
            else
                inside = inside & (x0[j] <= y1[j] <= x1[j])
            end
        end

        intrsct = intrsct | inside
    end

    return intrsct
end

"""    
Test of the rectanges defined by x0,x1  and y0,y1 intersects
             x1
  +----------+ 
  |          |
  |   +----------+ y1
  |   |      |   |
  +----------+   |
 x0   |          |
      |          |
      +----------+
     y0
"""

intersect(x0,x1,y0,y1) = intersect_(x0,x1,y0,y1) || intersect_(y0,y1,x0,x1)

#qt = QT_(X)



isleaf(qt) = length(qt.children) == 0

"""
number of points per node
it is always zero for non-leaf nodes
"""
Base.length(qt::QT) = size(qt.points,1)

inside(qt::QT,y) = inside(qt.min,qt.max,y)
Base.intersect(qt::QT,y0,y1) = intersect(qt.min,qt.max,y0,y1)
Base.ndims{T,TA,N}(qt::QT{T,TA,N}) = N

function count(qt::QT)
    if isleaf(qt)
        return length(qt)
    else
        c = 0
        for child in qt.children
            c = c + count(child)
        end

        return c
    end
end    

"""
sucess = add!(qt,x,attrib,max_cap = 10)
Add point `x` with the attribute `attrib` to the quadtree `qt`.
`sucess` is true if `x`is within the bounds of the quadtree node `qt` (otherwise 
false and the point has not been added)
"""
function add!(qt::QT,x,attrib,max_cap = 10)

    if !inside(qt,x)
        return false
    else                
        if isleaf(qt)
            qt.points  = cat(1,qt.points,x')
            push!(qt.attribs,attrib)
            
            # split if necessary
            rsplit!(qt, max_cap)
            return true
        else
            # try to add to all children and returns on first sucessful
            for child in qt.children
                if add!(child,x,attrib,max_cap)
                    return true
                end
            end

            # should never happen
            return false
        end
    end
end

"""
split a single node
"""

function split!{T,TA,N}(qt::QT{T,TA,N})
    # N: number of dimenions
    
    if isleaf(qt)
        xcenter = (qt.max + qt.min)/2

        nchildren = 2^N
        children = Vector{QT{T,TA,N}}(nchildren)

        # bounds of child
        cmin = Vector{T}(N)
        cmax = Vector{T}(N)        
        sel = trues(size(qt.points,1))
        
        for i = 1:nchildren
            sel[:] = true

            for j = 1:N
                # all corners of a hypercube
                if bitget(i-1, j)
                    sel = sel .& (qt.points[:,j] .<= xcenter[j])
                    cmin[j] = qt.min[j]
                    cmax[j] = xcenter[j]                    
                else
                    sel = sel .& (qt.points[:,j] .> xcenter[j])
                    cmin[j] = xcenter[j]                    
                    cmax[j] = qt.max[j]
                end
            end
            
            children[i] = QT(qt.points[sel,:],copy(cmin),copy(cmax),qt.attribs[sel])
        end
        
        
        # # in 2D
        
        # children = Vector{QT}(4)
        
        # sel1 = (qt.points[:,1] .<= xcenter[1]) .& (qt.points[:,2] .<= xcenter[2]);
        # children[1] = QT(qt.points[sel1,:],qt.min,xcenter)
        
        # sel2 = (qt.points[:,1] .<= xcenter[1]) .& (qt.points[:,2] .> xcenter[2]);
        # children[2] = QT(qt.points[sel2,:],[qt.min[1],xcenter[2]],[xcenter[1],qt.max[2]])
        
        # sel3 = (qt.points[:,1]  .> xcenter[1]) .& (qt.points[:,2] .<= xcenter[2]);
        # children[3] = QT(qt.points[sel3,:],[xcenter[1],qt.min[2]],[qt.max[1],xcenter[2]])

        # sel4 = (qt.points[:,1]  .> xcenter[1]) .& (qt.points[:,2] .> xcenter[2]);
        # children[4] = QT(qt.points[sel4,:],[xcenter[1],xcenter[2]],[qt.max[1],qt.max[2]])

        qt.children = children
        
    end
end

"""
recursive split
"""

function rsplit!{T,TA,N}(qt::QT{T,TA,N}, max_cap = 10)


    
    if isleaf(qt)
        if length(qt) < max_cap
            # no enougth points, nothing to do
            return
        end
        
        # all points are equal, stop recursion
        if all(qt.points .== qt.points[1,:]')
            return
        end

        split!(qt)
    end

    for child in qt.children
        rsplit!(child,max_cap)
    end

end


"""
points,attribs = within(qt,min,max)
Search all points within a bounding box defined by the vectors `min` and `max`.
"""

function within{T,TA,N}(qt::QT{T,TA,N}, min, max)
    if !Base.intersect(qt, min, max)
        # nothing to do
        return Array{T,2}(0,N), TA[]
    end

    if isleaf(qt)
        sel = Vector{Bool}(length(qt))
        for i = 1:length(qt)
            sel[i] =  inside(min,max,qt.points[i,:])
        end

        return qt.points[sel,:], qt.attribs[sel]
    end
    
    points = Array{T,2}(0,N)
    attribs = TA[]
    
    for child in qt.children
        cpoints,cattribs = within(child, min, max)
        
        if size(cpoints,1) != 0
            points = vcat(points,cpoints)
            attribs = vcat(attribs,cattribs)
        end
    end
    
    #    if !Base.intersect(qt, min, max) & (size(points,1) > 0)
    #        @show points,min,max
    #        @show qt
    #    end
    return points,attribs
    
end

function qplot(qt::QT)
    plot([qt.min[1], qt.max[1], qt.max[1], qt.min[1], qt.min[1]],
         [qt.min[2], qt.min[2], qt.max[2], qt.max[2], qt.min[2]])
end

function rplot(qt::QT)
    qplot(qt)
    for child in qt.children
        #@show child
        rplot(child)
    end
end

function Base.show(io::IO,qt::QT; indent = "  ")
    if isleaf(qt)
        print_with_color(:green, io, indent,"Leaf $(length(qt))")
    else
        print_with_color(:blue, io, indent,"Node ")
    end
    print(io,"  from $(qt.min) to $(qt.max)\n")

    if !isleaf(qt)
        for child in qt.children
            show(io,child; indent = indent * "  ")
        end
    end
end


function test()
    @testset "quadtree" begin
        @test [bitget(42,i) for i = 6:-1:1] == [true,false,true,false,true,false]
        
        @test inside([0,0],[1,1],[0.5,0.5]) == true
        @test inside([0,0],[1,1],[1.5,1.5]) == false
        

        @test intersect([0,0],[1,1],[0.5,0.5],[2,2]) == true
        @test intersect([0,0],[1,1],[1.5,1.5],[2,2]) == false
        # one rectange contains the other
        @test intersect([0,0],[1,1],[-1,-1],[2,2]) == true
        
        qt = QT(Int,[0.,0.],[1.,1.])
        add!(qt,[0.1,0.1],1)
        add!(qt,[0.2,0.2],2)
        add!(qt,[0.7,0.7],3)
        add!(qt,[0.9,0.1],4)

        split!(qt)


        @test isleaf(qt) == false
        @test isleaf(qt.children[1]) == true


        #X = rand(1000000,2)
        X = rand(10000,2)
        attribs = collect(1:size(X,1))

        @time begin
            qt2 = QT(X,attribs)
            rsplit!(qt2,5)
        end

        @time begin
            qt2 = QT(X,attribs)
            rsplit!(qt2,5)
        end



        #rplot(qt2)
        #plot(X[:,1],X[:,2],"b.")

        xmin = [0.3,0.3]
        xmax = [0.51,0.51]

        function simplesearch(X,attribs,xmin,xmax)
            sel = trues(size(X,1))
            for j = 1:size(X,2)
                sel = sel .& (xmin[j] .<= X[:,j] .<= xmax[j])
            end
            ind = find(sel)
            return X[ind,:],attribs[ind]
        end


        @time xref,attribs_ref = simplesearch(X,attribs,xmin,xmax)

        @time xin,attribs_res = within(qt2,xmin,xmax)
        @time xin,attribs_res = within(qt2,xmin,xmax)


        #@show xin
        #@show xref

        @test sortrows(xref) == sortrows(xin)
        @test sort(attribs_ref) == sort(attribs_res)


        # progressively add all points
        qt3 = QT(Int,[0.,0.],[1.,1.])

        for i = 1:size(X,1)
            add!(qt3,X[i,:],i)
        end

        @time xin,attribs_res = within(qt3,xmin,xmax)
        @test sortrows(xref) == sortrows(xin)
        @test sort(attribs_ref) == sort(attribs_res)


        # Test in 1D - 4D

        for n = 1:4
            @show n
            X = rand(100,n)
            attribs = collect(1:size(X,1))

            qtND = QT(X,attribs)
            rsplit!(qtND)

            @test ndims(qtND) == n
            @test count(qtND) == size(X,1)
            
            xmin = fill(0.0,(n,))
            xmax = fill(0.5,(n,))
            
            @time xref,attribs_ref = simplesearch(X,attribs,xmin,xmax)
            @time xin,attribs_res = within(qtND,xmin,xmax)

            @test sortrows(xref) == sortrows(xin)
            @test sort(attribs_ref) == sort(attribs_res)

            s = IOBuffer()
            show(s,qtND)
            @test contains(String(take!(s)),"Node")

        end

    end
end # function test()

export  QT, rsplit!, add!, show, ndims, count
    
end # module
