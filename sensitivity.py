'''
Script to run drifters backward from Galveston Bay to examine the Bay's 
connectivity with the shelf region.
'''

import matplotlib as mpl
mpl.use("Agg") # set matplotlib to use the backend that does not require a windowing system
import numpy as np
import os
import netCDF4 as netCDF
import pdb
import matplotlib.pyplot as plt
import tracpy
import init
from datetime import datetime, timedelta
import glob
from matplotlib.mlab import find

# def distance(pt, pts):
#     '''
#     Return distance between a point, pt, and an array of points, pts.
#     pt      [x, y]
#     pts     [xs, ys]
#     '''

#     # return np.sqrt((pt[0,:]-pts[0,:])**2 + (pt[1,:]-pts[1,:])**2)
#     # Leave as squared since that is what the dispersion calculation uses
#     # and it doesn't change the index order
#     pt = pt#/1000.
#     pts = pts#/1000.
#     return ((pt[0]-pts[0,:])**2 + (pt[1]-pts[1,:])**2) # in kilometers

# function to compute great circle distance between point lat1 and lon1 
# and arrays of points given by lons, lats or both same length arrays
# Haversine formula
def get_dist(lon1,lons,lat1,lats): 
    lon1 = lon1*np.pi/180.
    lons = lons*np.pi/180.
    lat1 = lat1*np.pi/180.
    lats = lats*np.pi/180.

    earth_radius = 6373.
    distance = earth_radius*2.0*np.arcsin(np.sqrt(np.sin(0.50*(lat1-lats))**2 \
                                       + np.cos(lat1)*np.cos(lats) \
                                       * np.sin(0.50*(lon1-lons))**2))
    return distance

def calc_dispersion(name, grid=None, r=1.):


    if grid is None:
        loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'
        grid = tracpy.inout.readgrid(loc)
    else:
        grid = grid

    # Read in tracks
    d = netCDF.Dataset(name)
    # d = netCDF.Dataset('tracks/2006-02-01T00C_doturb2_ah20.nc')
    lonp = d.variables['lonp'][:]
    latp = d.variables['latp'][:]
    t = d.variables['tp'][:]
    d.close()

    dist = np.zeros((lonp.shape[0],lonp.shape[0]))
    for idrifter in xrange(lonp.shape[0]):
        # dist contains all of the distances from other drifters for each drifter
        dist[idrifter,:] = get_dist(lonp[idrifter,0], lonp[:,0], latp[idrifter,0], latp[:,0])

    # let the index in axis 0 be the drifter id
    ID = np.arange(lonp.shape[0])

    pairs = []
    for idrifter in xrange(lonp.shape[0]):
        ind = find(dist[idrifter,:]<=r)
        for i in ind:
            if ID[idrifter] != ID[i]:
                pairs.append([min(ID[idrifter], ID[i]), 
                                max(ID[idrifter], ID[i])])

    pairs_set = set(map(tuple,pairs))
    pairs = map(list,pairs_set)# now pairs has only unique pairs of drifters
    # pairs.sort() #unnecessary but handy for checking work

    D2 = np.ones(lonp.shape[1])*np.nan
    nnans = np.zeros(lonp.shape[1]) # to collect number of non-nans over all drifters for a time
    for ipair in xrange(len(pairs)):
        dist = get_dist(lonp[pairs[ipair][0],:], lonp[pairs[ipair][1],:], 
                    latp[pairs[ipair][0],:], latp[pairs[ipair][1],:])
        D2 = np.nansum(np.vstack([D2, dist**2]), axis=0)
        nnans = nnans + ~np.isnan(dist)
    D2 = D2.squeeze()/nnans #len(pairs) # average over all pairs

    return D2, t

def run():

    # Location of TXLA model output
    loc = 'http://barataria.tamu.edu:8080/thredds/dodsC/NcML/txla_nesting6.nc'

    # Make sure necessary directories exist
    if not os.path.exists('tracks'):
        os.makedirs('tracks')
    if not os.path.exists('figures'):
        os.makedirs('figures')

    grid = tracpy.inout.readgrid(loc)

    # startdates = np.array([datetime(2006, 7, 1, 0, 1)])
    date = np.array([datetime(2006, 2, 1, 0, 1)])#, datetime(2006, 7, 1, 0, 1)])
    # pdb.set_trace()

    # drifter initial separation distance in meters
    dxs = np.array([20000,15000,10000,7500,5000,2500,1000,750,500,250])
    # Max volume a drifter can represent
    Vs = np.array([100,75,50,25])
    # Drifter array test areas
    Seeds = ['A','B','C']


    # loop through run parameters
    for dx in dxs:
        for Vmax in Vs:
            for Seed in Seeds:

                # Read in simulation initialization
                nstep, ndays, ff, tseas, ah, av, lon0, lat0, z0, zpar, do3d, doturb, \
                        grid, dostream, N, T0, U, V, name = init.disp(Seed, dx, Vmax, 
                                                            date, loc, grid=grid)

                # If the particle trajectories have not been run, run them
                if not os.path.exists('tracks/' + name + '.nc'):

                    # Run tracpy
                    lonp, latp, zp, t, grid, T0, U, V \
                        = tracpy.run.run(loc, nstep, ndays, ff, date, tseas, ah, av, \
                                            lon0, lat0, z0, zpar, do3d, doturb, name, \
                                            grid=grid, dostream=dostream, N=N, T0=T0, U=U, V=V)

                # If basic figures don't exist, make them
                if not os.path.exists('figures/' + name + '*.png'):

                    # Read in and plot tracks
                    d = netCDF.Dataset('tracks/' + name + '.nc')
                    lonp = d.variables['lonp'][:]
                    latp = d.variables['latp'][:]
                    tracpy.plotting.tracks(lonp, latp, name, grid=grid)
                    # tracpy.plotting.hist(lonp, latp, name, grid=grid, which='hexbin')
                    d.close()
                    # Do transport plot
                    tracpy.plotting.transport(name='', fmod=name, extraname=name,
                        Title='Transport on Shelf', dmax=1.0)
                    plt.close("all")
   
        # # Do transport plot
        # tracpy.plotting.transport(name='', fmod=startdate.isoformat()[0:7] + '*', 
        #     extraname=startdate.isoformat()[0:7], Title='Transport on Shelf', dmax=1.0)


if __name__ == "__main__":
    run()    