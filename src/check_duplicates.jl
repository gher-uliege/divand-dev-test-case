
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
function check_duplicates(lon,lat,z,time,delta)
    time2 = [Dates.Millisecond(t - DateTime(1900,1,1)).value/24/60/60/1000 for t in time]
    X = [lon lat z time2]
    qt = Quadtrees.QT(X,collect(1:size(lon,1)))
    mult = Vector{Int}(size(X,1))
    duplicates = Vector{Vector{Int}}(0)
    delta2 = delta/2

    @fastmath @inbounds for i = 1:size(X,1)
    #for i = 1:size(X,1)
        p,index = Quadtrees.within(qt,X[i,:] - delta2,X[i,:] + delta2)
        mult = length(index)
        if length(index) > 1
            push!(duplicates,index)
            #@show index
        end
    end

    return mult,duplicates
end
