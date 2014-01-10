import tracpy
import tracpy.calcs
import numpy
import glob
import netCDF4 as netCDF
import numpy as np
import pdb
import op

seed = 'C'

loccalcs = 'calcs/'
loctracks = 'tracks/'

# Loop through tracks files and run calculations
Files = glob.glob(loctracks + seed + '*.nc')

# Highest res case for comparison
dc = netCDF.Dataset(loctracks + seed + '_dx250_V25.nc') # control case
Uc = dc.variables['U'][:]; Vc = dc.variables['V'][:]
Sc = np.sqrt(op.resize(Uc,1)**2 + op.resize(Vc,0)**2)
dc.close()

# lonpc = dc.variables['lonp'][:]; latpc = dc.variables['latp'][:]; tpc = dc.variables['tp'][:];

loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'
grid = tracpy.inout.readgrid(loc)  # grid file, Gulf model domain
# grid = tracpy.inout.readgrid(loc[1], vert_filename=loc[0], llcrnrlon=llcrnrlon, urcrnrlon=urcrnrlon,
#                                   llcrnrlat=llcrnrlat, urcrnrlat=urcrnrlat)  # grid file, Gulf model domain

for File in Files:

    print File

    # Read in tracks
    d = netCDF.Dataset(File)
    lonp = d.variables['lonp'][:]; latp = d.variables['latp'][:]; tp = d.variables['tp'][:];

    # Project into x and y coords
    xp, yp = grid['basemap'](lonp,latp)

    # 1st moments
    M1x, nnans = tracpy.calcs.moment1(xp)
    np.savez(loccalcs + File.split('/')[-1][:-3] + 'M1x.npz', M=M1x, t=tp, nnans=nnans)
    M1y, nnans = tracpy.calcs.moment1(yp)
    np.savez(loccalcs + File.split('/')[-1][:-3] + 'M1y.npz', M=M1y, t=tp, nnans=nnans)

    # 2nd moments
    M2x, nnans = tracpy.calcs.moment2(xp, M1x)
    np.savez(loccalcs + File.split('/')[-1][:-3] + 'M2x.npz', M=M2x, t=tp, nnans=nnans)
    M2y, nnans = tracpy.calcs.moment2(yp, M1y)
    np.savez(loccalcs + File.split('/')[-1][:-3] + 'M2y.npz', M=M2y, t=tp, nnans=nnans)

    # 3rd moments
    M3x, nnans = tracpy.calcs.moment3(xp, M1x)
    np.savez(loccalcs + File.split('/')[-1][:-3] + 'M3x.npz', M=M3x, t=tp, nnans=nnans)
    M3y, nnans = tracpy.calcs.moment3(yp, M1y)
    np.savez(loccalcs + File.split('/')[-1][:-3] + 'M3y.npz', M=M3y, t=tp, nnans=nnans)
 
    # # 4th moments
    # M4x, nnans = tracpy.calcs.moment4(xp, M1x)
    # np.savez(loccalcs + File.split('/')[-1][:-3] + 'M4x.npz', M=M4x, t=tp, nnans=nnans)
    # M4y, nnans = tracpy.calcs.moment4(yp, M1y)
    # np.savez(loccalcs + File.split('/')[-1][:-3] + 'M4y.npz', M=M4y, t=tp, nnans=nnans)


    ## Transport comparisons
    U = d.variables['U'][:]; V = d.variables['V'][:]
    S = np.sqrt(op.resize(U,1)**2 + op.resize(V,0)**2)
    # pdb.set_trace()
    diffnorm = np.linalg.norm(S-Sc, 'fro')
    np.savez(loccalcs + File.split('/')[-1][:-3] + 'diffnorm.npz', diffnorm=diffnorm)


    d.close()
