#!/usr/bin/python
import math
import numpy
from pytags_coord_tools import ll2utm, utm2ll
from scipy.spatial import Delaunay
import matplotlib.pyplot as plt




def grid_sr_dist(sta_lat, sta_lon, props):
        maxlon = float(props.getmaxlon())
        minlon = float(props.getminlon())
        maxlat = float(props.getmaxlat())
        minlat = float(props.getminlat())
        nlon = int(props.getnlon())
        nlat = int(props.getnlat())
        sta_length = props.getstationlength()

        lons = numpy.linspace(minlon,maxlon,num=nlon)
        lats = numpy.linspace(minlat,maxlat,num=nlat)
        ngp = nlon*nlat #total number of grid points

        xrs = numpy.zeros([sta_length,ngp])
        yrs = numpy.zeros([sta_length,ngp])


        LONS = numpy.zeros([ngp,1])
        LATS = numpy.zeros([ngp,1])

        n = 0
        for i in range (0, nlon):
                for j in range (0, nlat):
                        LONS[n,0] = lons[i]
                        LATS[n,0] = lats[j]
                        n=n+1

        for i in range (0, sta_length):
                for j in range (0, ngp):
                        (x1,y1) = ll2utm(sta_lon[i],sta_lat[i], numpy.mean(lons), numpy.mean(lats))
                        (x2,y2) = ll2utm(LONS[j],LATS[j], numpy.mean(lons), numpy.mean(lats))


                        xrs[i,j] = (x1-x2)
                        yrs[i,j] = (y1-y2)

        return (xrs,yrs,LONS,LATS,ngp)

def triangulate(sta_lon,sta_lat,props):
        sta_length = props.getstationlength()
        pts = numpy.zeros([sta_length,2])
        maxtrilength = float(props.getmaxtrilength())

        for i in range (0, sta_length):
                (x1,y1) = ll2utm(sta_lon[i],sta_lat[i], numpy.mean(sta_lon), numpy.mean(sta_lat))
                pts[i,0] = x1
                pts[i,1] = y1

        tri = Delaunay(pts)
        triangles = tri.simplices.copy()
        [tri_length,tri_length2] = numpy.shape(triangles)

        triindexlarge = numpy.zeros([tri_length,1])
        for i in range (0, tri_length):
                tridim1 = math.sqrt(math.pow(pts[triangles[i,0],0]-pts[triangles[i,1],0],2)+ math.pow(pts[triangles[i,0],1]-pts[triangles[i,1],1],2))/1000       
                tridim2 = math.sqrt(math.pow(pts[triangles[i,0],0]-pts[triangles[i,2],0],2)+ math.pow(pts[triangles[i,0],1]-pts[triangles[i,2],1],2))/1000 
                tridim3 = math.sqrt(math.pow(pts[triangles[i,1],0]-pts[triangles[i,2],0],2)+ math.pow(pts[triangles[i,1],1]-pts[triangles[i,2],1],2))/1000


                if (tridim1 >= maxtrilength) or (tridim2 >= maxtrilength) or (tridim3 >= maxtrilength):
                        triindexlarge[i,0] = 1
                        
        a1 = numpy.where(triindexlarge == 0)[0]
        a1 = numpy.array(a1)
        triangles_culled = triangles[a1,:]

        [tri_length,tri_length2] = numpy.shape(triangles_culled)

        return (triangles_culled, tri_length)

