
"""
    dupl = check_duplicates(x,delta)

Based the coordinates `x` (a tuple of longitude `lons`, latitudes `lats`, depth (`zs`) 
and time (`times`)) check of points who are in the same spatio-temporal bounding
 box of a length `delta`. `delta` is a vector with 4 elements corresponding to 
longitude, latitude, depth and time
(in days). `mult` is a vector of the same length than `lons` with the number of
time an observation is present within this bounding box and `duplicates` a list
of duplicates and each element of `duplicates` corresponds to the index in
`lons`, `lats`, `zs` and times`.
"""

function check_duplicates(x,delta; maxcap = 100, label = collect(1:size(x[1],1)))
    n = length(x)
    Nobs = length(x[1])

    X = Array{Float64,2}(n,Nobs)
    for i = 1:n
        if eltype(x[i]) <: DateTime
            for j = 1:Nobs
                X[i,j] = Dates.Millisecond(x[i][j] - DateTime(1900,1,1)).value/24/60/60/1000
            end
        else
            X[i,:] = x[i]
        end
    end

    @show size(X)
    qt = Quadtrees.QTnew(X,label)
    @time Quadtrees.rsplit!(qt, maxcap)
    
    #mult = Vector{Int}(size(X,1))
    duplicates = Vector{Set{Int}}(0)
    delta2 = delta/2

    xmin = zeros(n)
    xmax = zeros(n)
    
    @time @fastmath @inbounds for i = 1:Nobs
        for j = 1:n
            xmin[j] = X[j,i] - delta2[j]
            xmax[j] = X[j,i] + delta2[j]
        end

        index = Quadtrees.within(qt,xmin,xmax)

        #@show index
        #mult = length(index)
        if length(index) > 1
            push!(duplicates,Set(index))
            #@show index
        end
    end

    #return mult,duplicates
    return collect(Set(duplicates))
end

function dupset(duplicates)
    d = Vector{Set{Int}}()
    sizehint!(d,length(duplicates) ÷ 2)
    
    for i = 1:length(duplicates)
        if !(duplicates[i] in d)
            push!(d,duplicates[i])
        end
    end

    return d
end





"""
    dupl = check_duplicates(x,delta)

Based the coordinates `x` (a tuple of longitude `lons`, latitudes `lats`, depth (`zs`) 
and time (`times`)) check of points who are in the same spatio-temporal bounding
 box of a length `delta`. `delta` is a vector with 4 elements corresponding to 
longitude, latitude, depth and time
(in days). `mult` is a vector of the same length than `lons` with the number of
time an observation is present within this bounding box and `duplicates` a list
of duplicates and each element of `duplicates` corresponds to the index in
`lons`, `lats`, `zs` and times`.
"""

function check_duplicatesv(x,value,delta,deltavalue; maxcap = 100, label = collect(1:size(x[1],1)))
    n = length(x)
    Nobs = length(x[1])

    X = Array{Float64,2}(n,Nobs)
    for i = 1:n
        if eltype(x[i]) <: DateTime
            for j = 1:Nobs
                X[i,j] = Dates.Millisecond(x[i][j] - DateTime(1900,1,1)).value/24/60/60/1000
            end
        else
            X[i,:] = x[i]
        end
    end

    @show size(X)
    qt = Quadtrees.QTnew(X,label)
    @time Quadtrees.rsplit!(qt, maxcap)
    
    #mult = Vector{Int}(size(X,1))
    duplicates = Vector{Set{Int}}(0)
    delta2 = delta/2

    xmin = zeros(n)
    xmax = zeros(n)
    
    @time @fastmath @inbounds for i = 1:Nobs
        for j = 1:n
            xmin[j] = X[j,i] - delta2[j]
            xmax[j] = X[j,i] + delta2[j]
        end

        index = Quadtrees.within(qt,xmin,xmax)

        #@show index

        if length(index) > 1
            # check for values
            vv = value[index]
            ii = sortperm(vv)
            
            istart = 1
            for i=1:length(vv)-1
                if vv[ii[i+1]] - vv[ii[i]] > deltavalue
                    #@show istart:i;

                    if i > istart
                        push!(duplicates,Set(index[ii[istart:i]]))
                    end
                    
                    istart=i+1
                end
            end
            i = length(vv)
            if i > istart
                push!(duplicates,Set(index[ii[istart:i]]))
            end
            
            #push!(duplicates,Set(index))
            #@show index
        end
    end

    #return mult,duplicates
    #return dupset(duplicates)
    return duplicates
end



"""
    dupl = checkduplicates(x1,value1,v2,value2,delta,deltavalue)

Based the coordinates `x` (a tuple of longitude `lons`, latitudes `lats`, depth (`zs`) 
and time (`times`)) check of points who are in the same spatio-temporal bounding
 box of a length `delta`. `delta` is a vector with 4 elements corresponding to 
longitude, latitude, depth and time
(in days). `mult` is a vector of the same length than `lons` with the number of
time an observation is present within this bounding box and `duplicates` a list
of duplicates and each element of `duplicates` corresponds to the index in
`lons`, `lats`, `zs` and times`.
"""

function catx(x)
    Nobs = length(x[1])

    X = Array{Float64,2}(n,Nobs)
    for i = 1:n
        if eltype(x[i]) <: DateTime
            for j = 1:Nobs
                X[i,j] = Dates.Millisecond(x[i][j] - DateTime(1900,1,1)).value/24/60/60/1000
            end
        else
            X[i,:] = x[i]
        end
    end
    return X
end

function checkduplicates(x1,value1,x2,value2,delta,deltavalue;
                         maxcap = 100, label = collect(1:size(x2[1],1)))
    n = length(x1)
    Nobs1 = length(x1[1])
    Nobs2 = length(x2[1])

    X1 = catx(x1)
    X2 = catx(x2)

    @show size(X1)
    qt = Quadtrees.QTnew(X1,label)
    @time Quadtrees.rsplit!(qt, maxcap)
    
    #mult = Vector{Int}(size(X,1))
    duplicates = Vector{Set{Int}}(0)
    delta2 = delta/2

    xmin = zeros(n)
    xmax = zeros(n)
    
    @time @fastmath @inbounds for i = 1:Nobs2
        for j = 1:n
            xmin[j] = X2[j,i] - delta2[j]
            xmax[j] = X2[j,i] + delta2[j]
        end

        index = Quadtrees.within(qt,xmin,xmax)

        #@show index

        if length(index) > 1
            # check for values
            vv = value[index]
            ii = sortperm(vv)
            
            istart = 1
            for i=1:length(vv)-1
                if vv[ii[i+1]] - vv[ii[i]] > deltavalue
                    #@show istart:i;

                    if i > istart
                        push!(duplicates,Set(index[ii[istart:i]]))
                    end
                    
                    istart=i+1
                end
            end
            i = length(vv)
            if i > istart
                push!(duplicates,Set(index[ii[istart:i]]))
            end
            
            #push!(duplicates,Set(index))
            #@show index
        end
    end

    #return mult,duplicates
    #return dupset(duplicates)
    return duplicates
end





function check_duplicates2(x,delta; maxcap = 100)
    function work(n,X,qt,delta2)
        function check(i)
            xmin = zeros(n)
            xmax = zeros(n)
            
            for j = 1:n
                xmin[j] = X[j,i] - delta2[j]
                xmax[j] = X[j,i] + delta2[j]
            end
            
            index = Quadtrees.within(qt,xmin,xmax)
            
            #mult = length(index)
            if length(index) > 1
                return index
            else
                return Int[]
            end
        end
        
        return check
    end
    
    n = length(x)
    Nobs = length(x[1])

    X = Array{Float64,2}(n,Nobs)
    for i = 1:n
        if eltype(x[i]) <: DateTime
            for j = 1:Nobs
                X[i,j] = Dates.Millisecond(x[i][j] - DateTime(1900,1,1)).value/24/60/60/1000
            end
        else
            X[i,:] = x[i]
        end
    end

    @show size(X)
    qt = Quadtrees.QTnew(X,collect(1:size(lon,1)))
    @time Quadtrees.rsplit!(qt, maxcap)
    
    #mult = Vector{Int}(size(X,1))
    duplicates = Vector{Vector{Int}}(0)
    delta2 = delta/2

    xmin = zeros(n)
    xmax = zeros(n)

    c = work(n,X,qt,delta2)
    dup = pmap(c,1:Nobs; batch_size = 1000)
#    @time @fastmath @inbounds for i = 1:Nobs
#        check(i)
#    end

    #return mult,duplicates
    return filter(x -> !isempty(x),dup)
end





function check_duplicates3(x,delta; maxcap = 100, label = collect(1:size(x[1],1)))
    n = length(x)
    Nobs = length(x[1])

    X = Array{Float64,2}(n,Nobs)
    for i = 1:n
        if eltype(x[i]) <: DateTime
            for j = 1:Nobs
                X[i,j] = Dates.Millisecond(x[i][j] - DateTime(1900,1,1)).value/24/60/60/1000
            end
        else
            X[i,:] = x[i]
        end
    end

    @show size(X)
    rtree = LibSpatialIndex.RTree(4)

    for i = 1:Nobs
        LibSpatialIndex.insert!(rtree, i, X[:,i])
    end
    
    qt = Quadtrees.QTnew(X,label)
    @time Quadtrees.rsplit!(qt, maxcap)
    
    #mult = Vector{Int}(size(X,1))
    duplicates = Vector{Set{Int}}(0)
    delta2 = delta/2

    xmin = zeros(n)
    xmax = zeros(n)
    
    @time @fastmath @inbounds for i = 1:Nobs
        for j = 1:n
            xmin[j] = X[j,i] - delta2[j]
            xmax[j] = X[j,i] + delta2[j]
        end

        index = Quadtrees.within(qt,xmin,xmax)

        #@show index
        #mult = length(index)
        if length(index) > 1
            push!(duplicates,Set(index))
            #@show index
        end
    end

    #return mult,duplicates
    return dupset(duplicates)
end
