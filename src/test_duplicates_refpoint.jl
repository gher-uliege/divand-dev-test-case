import divand

function duplicates_refpoint(lon,lat,depth,time,value)


dl=0.0005    # latitude (degree)
# 0.0005    # longitude
dz= 0.5          # metres
dt= 0.0208333 # jours (1/2 heurs)
dv=0.01
dtms=Dates.DateTime("1971-10-23T08:45:00.001")-Dates.DateTime("1971-10-23T08:15:00.001")

imax=size(value)[1]

distjm=zeros(Float64,imax)

iref=0

if iref==0
    mlon=mean(lon)
    mlat=mean(lat)
    mdepth=mean(depth)
    mvalue=mean(value)
else
    mlon=lon[iref]
    mlat=lat[iref]
    mdepth=depth[iref]
    mvalue=value[iref]
end


# critical distance, here for 4 dimensions pre-sorted via distjm
dellim=4.1

# distances with respect to a reference point
distjm=((lon.-mlon)/dl).^2+((lat.-mlat)/dl).^2+((depth.-mdepth)/dz).^2+((value.-mvalue)/dv).^2

# Sort them
jm=sortperm(distjm)



isduplicate = SharedArray{Int}(imax)
isduplicate[:]=0
#isduplicate=zeros(Int,imax);

# calculate upper bound for difference in distances beyond which distance is certain to be larger then dellim
jbound=dellim+2*sqrt.(dellim*distjm[jm[:]])


# Try out parallel run as whoever writes in a duplicate into the array marks it as duplicate.
# I have used addprocs(2) to get two parallel workers for the loop
# As a better balance is needed because work is not uniform , I used random indexes
# Seems to work very nicely

indexes=randperm(imax-1)
@sync @parallel for iii=1:imax-1
    i=indexes[iii]
    # instead of standard loop which would lead to unbalanced load
    #@sync @parallel for i=1:imax-1
    # iii=i

    ip=jm[i]
    j=i+1

    #if you want to check and make a total exhaustive research, comment the second part of the test
    while ((j<=imax) && ((distjm[jm[j]]-distjm[ip]) < jbound[j]  ))
        jp=jm[j]

        # Final test can be replaced by distance if considered better than in box decision
#dd=((lon[ip]-lon[jp])/dl)^2+((lat[ip]-lat[jp])/dl)^2+((depth[ip]-depth[jp])/dz)^2+((value[ip]-value[jp])/dv)^2
        #if (dd<4.01 && abs(time[ip]-time[jp])<dtms)


        # Room for optimization here: use first the test most likely to fail as then the rest is not done
        # So probably time and any other parameter which was not used in the preordering
        if ((abs(time[ip]-time[jp])<dtms) && (abs(value[ip]-value[jp])<dv) && (abs((lon[ip]-lon[jp]))<dl) && (abs(lat[ip]-lat[jp]) < dl) && (abs(depth[ip]-depth[jp])<dz) )

            # storing of duplicates in the following way makes sur
            # each data point which has a duplicate somewhere is marked and ONE of the duplicates is indicated
            # triplet should be easy to identify but quadruplets might be seen as two pairs of duplicates
            # either run again check after elimination or encode here another way to mark the duplicates
            # in this case, be carefull when thinking about parallel execution
            isduplicate[jp]=ip
            isduplicate[ip]=jp
            #@show j,i,ip,jp,distjm[jp],distjm[ip]
        end
        j=j+1


    end
    if mod(iii,10000)==0
        @show j-i,iii,distjm[jm[j]]
    end
end

  

    return isduplicate
end


obsname = joinpath(dirname(@__FILE__),"..","..","divand-example-data","Provencal","WOD-Salinity.nc"))
value,lon,lat,depth,obstime,ids = divand.loadobs(Float64,obsname,"Salinity")


# warm-up
sel = [1:1_000; 1];
dep = duplicates_refpoint(lon[sel],lat[sel],depth[sel],obstime[sel],value[sel]);

tic()
isduplicate  = duplicates_refpoint(lon,lat,depth,obstime,value)
toc() 
ii = isduplicate.>0
@show size(value[ii])
