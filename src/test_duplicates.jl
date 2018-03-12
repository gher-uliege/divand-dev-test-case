import divand
push!(LOAD_PATH,"/home/abarth/src/Mustache.jl/src");
import divand
import Quadtrees

include("check_duplicates.jl")

obsname = "/home/abarth/projects/Julia/divand.jl/examples/WOD-Salinity.nc"
value,lon,lat,depth,obstime,ids = divand.loadobs(Float64,obsname,"Salinity")


sel = [1:10_000; 1]

sel = [1:100; 1]
#@time count = check_duplicates2((lon[sel],lat[sel],depth[sel],obstime[sel]),[0.001,0.001,1.,1./24]);

@time dup = check_duplicates((lon[sel],lat[sel],depth[sel],obstime[sel]),[0.001,0.001,1.,1./24]);
#Profile.clear_malloc_data()
#@time dup = check_duplicates(lon[sel],lat[sel],depth[sel],obstime[sel],[0.001,0.001,1.,1./24]);

#mult,duplicates = check_duplicates(lon[sel],lat[sel],depth[sel],time[sel],[0.1,0.1,1.,0.1])


#mult,duplicates = check_duplicates([0,1,1,0,0],[0,0,lat[sel],depth[sel],obstime[sel],[0.1,0.1,1.,0.1])


# sel = 1:10000

# @time mult,duplicates = check_duplicates(lon[sel],lat[sel],depth[sel],obstime[sel],[0.1,0.1,1.,0.1])

# lont,latt,zt,obstimet,deltat = lon[sel],lat[sel],depth[sel],obstime[sel],[0.1,0.1,1.,0.1]

# obstime2 = [Dates.Millisecond(t - DateTime(1900,1,1)).value/24/60/60/1000 for t in obstimet]
#     X = [lont latt zt obstime2]
#     @time qt = Quadtrees.QT(X,collect(1:size(lont,1)))
#     mult = Vector{Int}(size(X,1))
#     duplicates = Vector{Vector{Int}}(0)
#     delta2 = delta/2


#9.412261 
# 9.412261 seconds (295.07 M allocations: 13.439 GiB, 21.03% gc time)
