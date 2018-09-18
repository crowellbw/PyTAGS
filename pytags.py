#!/usr/bin/python
import math
import numpy
import csv
from pytags_coord_tools import ll2utm
from pytags_spatial_tools import grid_sr_dist, triangulate
from pytags_straincalcs import strain, strain_tri
from pytags_paraminit import Properties
from pathlib2 import Path
import calendar
import matplotlib.pyplot as plt
from datetime import datetime, timedelta
#######################################################################
#PyTAGS - A Python package for the Temporal Analysis of GNSS Strains
#Inversion of interseismic GPS displacement or velocities for strain 
#and translation in grid or in triangulated subnetworks
#Written by Brendan Crowell, University of Washington
#v1.0, last modified November 7, 2017
#######################################################################
#Parameters
#######################################################################
props=Properties('pytags.props')

strainmode = props.getstrainmode()
temporalmode = props.gettemporalmode()

if int(strainmode) == 0:
    print ("Using triangulated strain mode")
elif int(strainmode) == 1:
    print ("Using gridded strain mode")
else:
    print ("Incorrect strain mode - killing run - please enter 0 or 1 in pytags.props, you entered",int(strainmode))
    quit()

if int(temporalmode) == 0:
    print ("Using time series temporal mode")
elif int(temporalmode) == 1:
    print ("Using static velocity temporal mode")
else:
    print ("Incorrect temporal mode - killing run - please enter 0 or 1 in pytags.props, you entered",int(temporalmode))
    quit()


sitefile = props.getsitefile() #stations to use

startyear = props.getstartyear()
endyear = props.getendyear()
datadir = props.getdatadir()
datasource = props.getdatasource()
datasuffix = props.getdatasuffix()

#########################################################################
###Create time matrices
#########################################################################

if int(temporalmode) == 0:
    ndays = 0
    y365 = numpy.ones([365,1])
    y366 = numpy.ones([366,1])
    d365 = numpy.linspace(0,364,365)[numpy.newaxis]
    d366 = numpy.linspace(0,365,366)[numpy.newaxis]



    for i in range (int(startyear), int(endyear)+1):
        isleap = calendar.isleap(i)
        if str(isleap) == 'True':
            ndays = ndays+366
        else:
            ndays = ndays+365
            

    isleap = calendar.isleap(int(startyear))
    if str(isleap) == 'True':
        YMAT = int(startyear)*y366
        DMAT = numpy.transpose(d366)
    else:
        YMAT = int(startyear)*y365
        DMAT = numpy.transpose(d365)

    for i in range (int(startyear)+1, int(endyear)+1):
        isleap = calendar.isleap(i)
        if str(isleap) == 'True':
            YMAT = numpy.vstack([YMAT,i*y366])
            DMAT = numpy.vstack([DMAT,numpy.transpose(d366)])
        else:
            YMAT = numpy.vstack([YMAT,i*y365])
            DMAT = numpy.vstack([DMAT,numpy.transpose(d365)])

#######################################################################
#Load Station information and data into matrix
#######################################################################
sta_length = props.getstationlength()

sta_lat = numpy.zeros([sta_length,1])
sta_lon = numpy.zeros([sta_length,1])
sites=[]
if int(temporalmode) == 0:
    N = numpy.nan*numpy.ones([ndays,sta_length])
    E = numpy.nan*numpy.ones([ndays,sta_length])

elif int(temporalmode) == 1:
    N = numpy.nan*numpy.ones([sta_length,1])
    E = numpy.nan*numpy.ones([sta_length,1])


k=0
with open(sitefile, 'rt') as f:
    reader = csv.reader(f, delimiter=',')
    for row in reader:
        if not row[0].startswith("#"):
            sta_lat[k,0] = row[1]
            sta_lon[k,0] = row[0]
            site = row[7]
            sites.append(str.lower(site))
            if int(temporalmode) == 0:
            
                if datasource == 'sopac':
                    GPSfile = datadir + str.lower(site) + datasuffix
                    filetest = Path(GPSfile)
                    if filetest.is_file():
                        print ("loading", GPSfile)

                        with open(GPSfile, 'rt') as g:
                            rows = (line.split() for line in g)
                            for grow in rows:
                                if not grow[0].startswith("#"):
                                    ygps = grow[1]
                                    dgps = grow[2]
                                    north = float(grow[3])
                                    east = float(grow[4])

                                    a1 = numpy.where((int(ygps) == YMAT) & (int(dgps) == DMAT))[0]
                                    a1 = numpy.array(a1)
                                    if len(a1) > 0:
                                        N[a1[0],k] = north/1000
                                        E[a1[0],k] = east/1000
                    else:
                        print ("file", GPSfile, "does not exist")

                elif datasource == 'pbo':

                    GPSfile = datadir + str.upper(site) + datasuffix
                    filetest = Path(GPSfile)
                    if filetest.is_file():
                        print ("loading", GPSfile)

                        with open(GPSfile, 'rt') as g:
                            rows = (line.split() for line in g)
                            for grow in rows:
                                if grow[0].isdigit():
                                    
                                    mjd = float(grow[2])
                                    mjd_1 = datetime(1858,11,17,0,0,0)+timedelta(days=mjd)
                                    tt = mjd_1.timetuple()
                                    doy = tt.tm_yday-1
                                    yr = tt.tm_year
                                    north = float(grow[15])
                                    east = float(grow[16])


                                    a1 = numpy.where((int(yr) == YMAT) & (int(doy) == DMAT))[0]
                                    a1 = numpy.array(a1)
                                    if len(a1) > 0:
                                        N[a1[0],k] = north
                                        E[a1[0],k] = east

                    elif datasource == 'ngl':

                        GPSfile = datadir + str.upper(site) + datasuffix
                        filetest = Path(GPSfile)
                        if filetest.is_file():
                            print ("loading", GPSfile)

                            with open(GPSfile, 'rt') as g:
                                rows = (line.split() for line in g)
                                for grow in rows:
                                    if grow[0].isdigit():
                                        
                                        mjd = float(grow[2])
                                        mjd_1 = datetime(1858,11,17,0,0,0)+timedelta(days=mjd)
                                        tt = mjd_1.timetuple()
                                        doy = tt.tm_yday
                                        yr = tt.tm_year
                                        north = float(grow[15])
                                        east = float(grow[16])

                                        a1 = numpy.where((int(yr) == YMAT) & (int(doy) == DMAT))[0]
                                        a1 = numpy.array(a1)
                                        if len(a1) > 0:
                                            N[a1[0],k] = north
                                            E[a1[0],k] = east
                    else:
                        print ("file", GPSfile, "does not exist")
                else:
                    print ("Incorrect datasource, currently only accepting sopac or pbo formatted data")
                    exit()
            elif int(temporalmode) == 1:
                E[k,0] = row[2]
                N[k,0] = row[3]
                
        
            k=k+1



#########################################################################
###Compute strains if triangulation mode
#########################################################################
if int(strainmode) == 0:
    print ("Delaunay triangulating stations")
    [triangles,tri_length] = triangulate(sta_lon,sta_lat,props)

    
    sta_lat1 = numpy.zeros([tri_length,1])
    sta_lon1 = numpy.zeros([tri_length,1])

    sta_lat2 = numpy.zeros([tri_length,1])
    sta_lon2 = numpy.zeros([tri_length,1])

    sta_lat3 = numpy.zeros([tri_length,1])
    sta_lon3 = numpy.zeros([tri_length,1])

    tri_lat = numpy.zeros([tri_length,1])
    tri_lon = numpy.zeros([tri_length,1])

    for i in range (0,tri_length):
        ind1 = triangles[i,1]
        
        sta_lat1[i,0] = sta_lat[triangles[i,0],0]
        sta_lat2[i,0] = sta_lat[triangles[i,1],0]
        sta_lat3[i,0] = sta_lat[triangles[i,2],0]
        sta_lon1[i,0] = sta_lon[triangles[i,0],0]
        sta_lon2[i,0] = sta_lon[triangles[i,1],0]
        sta_lon3[i,0] = sta_lon[triangles[i,2],0]

    sta1 = []
    sta2 = []
    sta3 = []

    if int(temporalmode) == 0:
        N1 = numpy.nan*numpy.ones([ndays,tri_length])
        E1 = numpy.nan*numpy.ones([ndays,tri_length])

        N2 = numpy.nan*numpy.ones([ndays,tri_length])
        E2 = numpy.nan*numpy.ones([ndays,tri_length])

        N3 = numpy.nan*numpy.ones([ndays,tri_length])
        E3 = numpy.nan*numpy.ones([ndays,tri_length])

        for i in range (0, tri_length):
            ind1 = triangles[i,0]
            ind2 = triangles[i,1]
            ind3 = triangles[i,2]

            sta1.append(str.lower(sites[ind1]))
            sta2.append(str.lower(sites[ind2]))
            sta3.append(str.lower(sites[ind3]))
            
            for j in range (0, ndays):
                N1[j,i]=N[j,ind1]
                N2[j,i]=N[j,ind2]
                N3[j,i]=N[j,ind3]

                E1[j,i]=E[j,ind1]
                E2[j,i]=E[j,ind2]
                E3[j,i]=E[j,ind3]
                
        print ("Beginning strain calculations for", tri_length, "trianglular subnetworks")
        subnet_outfile = 'output/subnetworks_triangulated_overview.txt'
        fstrain_subnetworks = open(subnet_outfile,'w')
        fstrain_subnetworks.write('#'+'Site1'+' '+'Lon1'+' '+'Lat1'+' '+'Site2'+' '+'Lon2'+' '+'Lat2'+' '+'Site3'+' '+'Lon3'+' '+'Lat3'+'\n')
        for i in range (0, tri_length):
            outfile = 'output/triangle'+str(i+1)+'_'+ sta1[i]+ '_' + sta2[i] + '_' + sta3[i] + '.txt'
            fstrain = open(outfile,'w')

            print ("Subnetwork", i+1, "of", tri_length,"; stations", sta1[i], sta2[i], sta3[i])

            lat0 = sta_lat1[i,0]
            lon0 = sta_lon1[i,0]

            lat1 = sta_lat2[i,0]
            lon1 = sta_lon2[i,0]

            lat2 = sta_lat3[i,0]
            lon2 = sta_lon3[i,0]

            latp0 = "{0:.4f}".format(float(lat0))
            latp1 = "{0:.4f}".format(float(lat1))
            latp2 = "{0:.4f}".format(float(lat2))

            lonp0 = "{0:.4f}".format(float(lon0))
            lonp1 = "{0:.4f}".format(float(lon1))
            lonp2 = "{0:.4f}".format(float(lon2))

            fstrain.write('##########################################################################'+'\n')
            fstrain.write('#'+'PyTAGS Triangulated Subnetwork Strain overview'+'\n')
            fstrain.write('##########################################################################'+'\n')
            fstrain.write('#'+'Data source'+' '+props.getdatasource()+'\n')
            fstrain.write('#'+'Data file suffix'+' '+props.getdatasuffix()+'\n')
            fstrain.write('#'+'Run between'+' '+props.getstartyear()+' and '+props.getendyear()+'\n')
            fstrain.write('#'+'Triangle'+' '+str(i+1)+'\n')
            fstrain.write('#'+'Station1'+' '+sta1[i]+'\n')
            fstrain.write('#'+'Lon1'+' '+lonp0+'\n')
            fstrain.write('#'+'Lat1'+' '+latp0+'\n')
            fstrain.write('#'+'Station2'+' '+sta2[i]+'\n')
            fstrain.write('#'+'Lon2'+' '+lonp1+'\n')
            fstrain.write('#'+'Lat2'+' '+latp1+'\n')
            fstrain.write('#'+'Station3'+' '+sta3[i]+'\n')
            fstrain.write('#'+'Lon3'+' '+lonp2+'\n')
            fstrain.write('#'+'Lat3'+' '+latp2+'\n')
            fstrain.write('#'+'Year'+' '+'DOY'+' '+'Decimal time'+' '+'E1'+' '+'E2'+' '+'Een'+' '+'w'+' '+'theta'+'\n')
            fstrain.write('##########################################################################'+'\n')

            fstrain_subnetworks.write(sta1[i]+' '+lonp0+' '+latp0+' '+sta2[i]+' '+lonp1+' '+latp1+' '+sta3[i]+' '+lonp2+' '+latp2+' '+'\n')
           
            for j in range (0, ndays):
                if (math.isnan(N1[j,i]) == False and math.isnan(N2[j,i]) == False and math.isnan(N3[j,i]) == False):
                    n0 = N1[j,i]
                    e0 = E1[j,i]

                    n1 = N2[j,i]
                    e1 = E2[j,i]
                    
                    n2 = N3[j,i]
                    e2 = E3[j,i]

                    L = strain_tri(lon0,lat0,lon1,lat1,lon2,lat2,n0,e0,n1,e1,n2,e2)
                   
                    isleap = calendar.isleap(int(YMAT[j,0]))
                    if str(isleap) == 'True':
                        DTIME = YMAT[j,0] + (DMAT[j,0]+0.5)/366
                    else:
                        DTIME = YMAT[j,0] + (DMAT[j,0]+0.5)/365
                    timestr = "{0:.4f}".format(float(DTIME))
                    yearstr = "{0:.0f}".format(float(YMAT[j,0]))
                    daystr = "{0:.0f}".format(float(DMAT[j,0]))
                    
                    e1str = "{0:.4e}".format(float(L[0]))
                    e2str = "{0:.4e}".format(float(L[1]))
                    eenstr = "{0:.4e}".format(float(L[2]))
                    wstr = "{0:.4e}".format(float(L[4]))
                    thetastr = "{0:.2f}".format(float(L[3]*180/math.pi))

                    fstrain.write(yearstr+' '+daystr+' '+timestr+' '+e1str+' '+e2str+' '+eenstr+' '+wstr+' '+thetastr+'\n')
            fstrain.close()
        fstrain_subnetworks.close()
        
    elif int(temporalmode) == 1:
        N1 = numpy.nan*numpy.ones([tri_length,1])
        E1 = numpy.nan*numpy.ones([tri_length,1])

        N2 = numpy.nan*numpy.ones([tri_length,1])
        E2 = numpy.nan*numpy.ones([tri_length,1])

        N3 = numpy.nan*numpy.ones([tri_length,1])
        E3 = numpy.nan*numpy.ones([tri_length,1])

        for i in range (0, tri_length):
            ind1 = triangles[i,0]
            ind2 = triangles[i,1]
            ind3 = triangles[i,2]

            sta1.append(str.lower(sites[ind1]))
            sta2.append(str.lower(sites[ind2]))
            sta3.append(str.lower(sites[ind3]))

            N1[i,0]=N[ind1,0]
            N2[i,0]=N[ind2,0]
            N3[i,0]=N[ind3,0]

            E1[i,0]=E[ind1,0]
            E2[i,0]=E[ind2,0]
            E3[i,0]=E[ind3,0]
            
        print ("Beginning strain calculations for", tri_length, "trianglular subnetworks")
        outfile = 'output/strain_triangulated_staticmode.txt'
        fstrain = open(outfile,'w')
        fstrain.write('#'+'Site1'+' '+'Lon1'+' '+'Lat1'+' '+'Site2'+' '+'Lon2'+' '+'Lat2'+' '+'Site3'+' '+'Lon3'+' '+'Lat3'+' '+'E1'+' '+'E2'+' '+'Een'+' '+'w'+' '+'theta'+'\n')
        for i in range (0, tri_length):
            print ("Subnetwork", i+1, "of", tri_length,"; stations", sta1[i], sta2[i], sta3[i])

            lat0 = sta_lat1[i,0]
            lon0 = sta_lon1[i,0]

            lat1 = sta_lat2[i,0]
            lon1 = sta_lon2[i,0]

            lat2 = sta_lat3[i,0]
            lon2 = sta_lon3[i,0]

            latp0 = "{0:.4f}".format(float(lat0))
            latp1 = "{0:.4f}".format(float(lat1))
            latp2 = "{0:.4f}".format(float(lat2))

            lonp0 = "{0:.4f}".format(float(lon0))
            lonp1 = "{0:.4f}".format(float(lon1))
            lonp2 = "{0:.4f}".format(float(lon2))


            if (math.isnan(N1[i,0]) == False and math.isnan(N2[i,0]) == False and math.isnan(N3[i,0]) == False):
                n0 = N1[i,0]
                e0 = E1[i,0]

                n1 = N2[i,0]
                e1 = E2[i,0]
                
                n2 = N3[i,0]
                e2 = E3[i,0]

                L = strain_tri(lon0,lat0,lon1,lat1,lon2,lat2,n0,e0,n1,e1,n2,e2)                  
            

                e1str = "{0:.4e}".format(float(L[0]))
                e2str = "{0:.4e}".format(float(L[1]))
                eenstr = "{0:.4e}".format(float(L[2]))
                wstr = "{0:.4e}".format(float(L[4]))
                thetastr = "{0:.2f}".format(float(L[3]*180/math.pi))


                fstrain.write(sta1[i]+' '+lonp0+' '+latp0+' '+sta2[i]+' '+lonp1+' '+latp1+' '+sta3[i]+' '+lonp2+' '+latp2+' '+e1str+' '+e2str+' '+eenstr+' '+wstr+' '+thetastr+'\n')
        fstrain.close()
    
#########################################################################
###Compute strains for grid mode
#########################################################################
if int(strainmode) == 1:
    print ("Computing source receiver distances in grid")
    [xrs,yrs,LONS,LATS,ngp] = grid_sr_dist(sta_lat,sta_lon,props)
    if int(temporalmode) == 0:
        for j in range (0, ngp):
            latstr2 = "{0:.0f}".format(float(LATS[j,0]*10000))
            lonstr2 = "{0:.0f}".format(float(LONS[j,0]*10000))
            latstr = "{0:.4f}".format(float(LATS[j,0]))
            lonstr = "{0:.4f}".format(float(LONS[j,0]))
            print ("Computing strain for grid", j+1, "of", ngp)
            outfile = 'output/grid'+ str(j+1)+'_'+latstr2+'_'+lonstr2+'.txt'
            fstrain = open(outfile,'w')
            fstrain.write('#'+'Year'+' '+'DOY'+' '+'Decimal time'+' '+'Lon,'+' '+'Lat,'+' '+'E_Translation,'+' '+'N_translation,'+' '+'E1,'+' '+'E2,'+' '+'Exy,'+' '+'w;'+' '+'(Units determined by input data - assumed to be m or m/yr)'+'\n')
            for i in range (0, ndays):
                a1 = numpy.where(numpy.isnan(N[i,:]) == False)[0]
                a1 = numpy.array(a1)
                if len(a1) > 0:
                    xs = xrs[a1,j]
                    ys = yrs[a1,j]
                    nd = N[i,a1]
                    ed = E[i,a1]

                    L = strain(xs,ys,ed,nd,props)

                    isleap = calendar.isleap(int(YMAT[i,0]))
                    if str(isleap) == 'True':
                        DTIME = YMAT[i,0] + (DMAT[i,0]+0.5)/366
                    else:
                        DTIME = YMAT[i,0] + (DMAT[i,0]+0.5)/365
                        
                    timestr = "{0:.4f}".format(float(DTIME))
                    
                    yearstr = "{0:.0f}".format(float(YMAT[i,0]))
                    daystr = "{0:.0f}".format(float(DMAT[i,0]))

                    dxstr = "{0:.4e}".format(float(L[0]))
                    dystr = "{0:.4e}".format(float(L[1]))
                    e1str = "{0:.4e}".format(float(L[2]))
                    e2str = "{0:.4e}".format(float(L[3]))
                    exystr = "{0:.4e}".format(float(L[4]))
                    wstr = "{0:.4e}".format(float(L[5]))
                    theta = "{0:.2f}".format(float(L[6]*180/math.pi))                              

                    fstrain.write(yearstr+' '+daystr+' '+timestr+' '+lonstr+' '+latstr+' '+dxstr+' '+dystr+' '+e1str+' '+e2str+' '+exystr+' '+wstr+' '+theta+'\n')
            fstrain.close()

        print ("Finishing writing strains to file")
    elif int(temporalmode) == 1:
        outfile = 'output/strain_gridded_staticmode.txt'
        fstrain = open(outfile,'w')
        fstrain.write('#'+'Lon,'+' '+'Lat,'+' '+'E_Translation,'+' '+'N_translation,'+' '+'E1,'+' '+'E2,'+' '+'Exy,'+' '+'w;'+' '+'(Units determined by input data - assumed to be m or m/yr)'+'\n')
        for j in range (0, ngp):
            a1 = numpy.where(numpy.isnan(N[:,0]) == False)[0]
            a1 = numpy.array(a1)
            if len(a1) > 0:
                xs = xrs[a1,j]
                ys = yrs[a1,j]
                nd = N[a1,0]
                ed = E[a1,0]

                L = strain(xs,ys,ed,nd,props)
                

                latstr = "{0:.4f}".format(float(LATS[j,0]))
                lonstr = "{0:.4f}".format(float(LONS[j,0]))
                dxstr = "{0:.4e}".format(float(L[0]))
                dystr = "{0:.4e}".format(float(L[1]))
                e1str = "{0:.4e}".format(float(L[2]))
                e2str = "{0:.4e}".format(float(L[3]))
                exystr = "{0:.4e}".format(float(L[4]))
                wstr = "{0:.4e}".format(float(L[5]))
                theta = "{0:.2f}".format(float(L[6]*180/math.pi))
                            

                fstrain.write(lonstr+' '+latstr+' '+dxstr+' '+dystr+' '+e1str+' '+e2str+' '+exystr+' '+wstr+'\n')
        fstrain.close()



