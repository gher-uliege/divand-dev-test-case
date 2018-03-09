
"""
    mult,duplicates = check_duplicates(lons,lats,zs,times,delta)

Based on longitude `lons`, latitudes `lats`, depth (`zs`) and time (`times`) check of points
who are in the same spatio-temporal bounding box of a length `delta`. `delta` is
a vector with 4 elements corresponding to longitude, latitude, depth and time
(in days). `mult` is a vector of the same length than `lons` with the number of
time an observation is present within this bounding box and `duplicates` a list
of duplicates and each element of `duplicates` corresponds to the index in
`lons`, `lats`, `zs` and times`.
"""

function check_duplicates(x,delta; maxcap = 100)
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
    
    @time @fastmath @inbounds for i = 1:Nobs
        for j = 1:n
            xmin[j] = X[j,i] - delta2[j]
            xmax[j] = X[j,i] + delta2[j]
        end

        index = Quadtrees.within(qt,xmin,xmax)

        #mult = length(index)
        if length(index) > 1
            push!(duplicates,index)
            #@show index
        end
    end

    #return mult,duplicates
    return duplicates
end




function check_duplicates2(x,delta)
    lon,lat,z,time = x
    time2 = [Dates.Millisecond(t - DateTime(1900,1,1)).value/24/60/60/1000 for t in time]
    X = [lon lat z time2]
    n = size(X,2)

    @time qt = Quadtrees.QT(X,collect(1:size(lon,1)))
    @time Quadtrees.rsplit!(qt, 100)
    
    count = zeros(Int,size(X,1))
    delta2 = delta/2

    xmin = zeros(n)
    xmax = zeros(n)
    
    @time     @fastmath @inbounds for i = 1:size(X,1)
        for j = 1:n
            xmin[j] = X[i,j] - delta2[j]
            xmax[j] = X[i,j] + delta2[j]
        end
        #p,index = Quadtrees.within(qt,X[i,:] - delta2,X[i,:] + delta2)
        
        Quadtrees.withincount!(qt,xmin,xmax,count)
    end

    return count
end
