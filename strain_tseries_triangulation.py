#!/usr/bin/python
import math
import numpy
import csv
from coord_tools import ll2utm
from strain_calc import strain_tri
import calendar

#######################################################################
#Inversion of interseismic GPS velocities for strain and translation in grid
#Written by Brendan Crowell, University of Washington
#5/19/2017
#######################################################################
#Parameters
#######################################################################
sitefile = 'cascadia_triangles.txt' #stations to use

#outfile = 'cascadia_strain_output.txt
#fstrain = open(outfile,'w')
startyear = 2000
endyear = 2017
datadir = 'data/'

#########################################################################
###Create data matrices
#########################################################################
ndays = 0
y365 = numpy.ones([365,1])
y366 = numpy.ones([366,1])
d365 = numpy.linspace(0,364,365)[numpy.newaxis]
d366 = numpy.linspace(0,365,366)[numpy.newaxis]


for i in range (startyear, endyear+1):
    isleap = calendar.isleap(i)
    if isleap == 'True':
        ndays = ndays+366
    else:
        ndays = ndays+365

isleap = calendar.isleap(startyear)
if isleap == 'True':
    YMAT = startyear*y366
    DMAT = numpy.transpose(d366)
else:
    YMAT = startyear*y365
    DMAT = numpy.transpose(d365)

for i in range (startyear+1, endyear+1):
    isleap = calendar.isleap(i)
    if isleap == 'True':
        YMAT = numpy.vstack([YMAT,i*y366])
        DMAT = numpy.vstack([DMAT,numpy.transpose(d366)])
    else:
        YMAT = numpy.vstack([YMAT,i*y365])
        DMAT = numpy.vstack([DMAT,numpy.transpose(d365)])


#######################################################################
#Load Station information and data
#######################################################################
k=0
with open(sitefile) as f:
    for i, l in enumerate(f):
        if not l.startswith("#"):
            k = k+1

sta_length = k

sta_lat1 = numpy.zeros([sta_length,1])
sta_lon1 = numpy.zeros([sta_length,1])

sta_lat2 = numpy.zeros([sta_length,1])
sta_lon2 = numpy.zeros([sta_length,1])

sta_lat3 = numpy.zeros([sta_length,1])
sta_lon3 = numpy.zeros([sta_length,1])

tri_lat = numpy.zeros([sta_length,1])
tri_lon = numpy.zeros([sta_length,1])



T = numpy.nan*numpy.ones([ndays,sta_length])

sta1 = []
sta2 = []
sta3 = []


N1 = numpy.nan*numpy.ones([ndays,sta_length])
E1 = numpy.nan*numpy.ones([ndays,sta_length])

N2 = numpy.nan*numpy.ones([ndays,sta_length])
E2 = numpy.nan*numpy.ones([ndays,sta_length])

N3 = numpy.nan*numpy.ones([ndays,sta_length])
E3 = numpy.nan*numpy.ones([ndays,sta_length])


k=0
with open(sitefile, 'rt') as f:
    reader = csv.reader(f, delimiter=',')
    for row in reader:
        if not row[0].startswith("#"):
            site1 = row[0]
            site2 = row[1]
            site3 = row[2]

            sta1.append(str.lower(site1))
            sta2.append(str.lower(site2))
            sta3.append(str.lower(site3))

            tri_lat[k,0] = row[3]
            tri_lon[k,0] = row[4]

            sta_lat1[k,0] = row[5]
            sta_lon1[k,0] = row[6]

            sta_lat2[k,0] = row[7]
            sta_lon2[k,0] = row[8]
            
            sta_lat3[k,0] = row[9]
            sta_lon3[k,0] = row[10]


            GPSfile1 = datadir + str.lower(site1) + 'FilterResid.neu'
            with open(GPSfile1, 'rt') as g:
                rows = (line.split() for line in g)
                for grow in rows:
                    if not grow[0].startswith("#"):
                        t = grow[0]
                        ygps = grow[1]
                        dgps = grow[2]
                        north = float(grow[3])
                        east = float(grow[4])

                        a1 = numpy.where((int(ygps) == YMAT) & (int(dgps) == DMAT))[0]
                        a1 = numpy.array(a1)
                        if len(a1) > 0:
                            N1[a1[0],k] = north
                            E1[a1[0],k] = east
                            T[a1[0],k] = t

            GPSfile2 = datadir + str.lower(site2) + 'FilterResid.neu'
            with open(GPSfile2, 'rt') as g:
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
                            N2[a1[0],k] = north
                            E2[a1[0],k] = east

            GPSfile3 = datadir + str.lower(site3) + 'FilterResid.neu'
            with open(GPSfile3, 'rt') as g:
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
                            N3[a1[0],k] = north
                            E3[a1[0],k] = east
                
            
            k=k+1


print sta1

#########################################################################
###Compute Strain, write out to file
#########################################################################


for i in range (0, ndays):
    k=0
    for j in range (0, sta_length):
        if (math.isnan(N1[i,j]) == False and math.isnan(N2[i,j]) == False and math.isnan(N3[i,j]) == False):

            timer = T[i,j]
            
            n0 = N1[i,j]
            e0 = E1[i,j]

            n1 = N2[i,j]
            e1 = E2[i,j]
            
            n2 = N3[i,j]
            e2 = E3[i,j]

            lat0 = sta_lat1[j,0]
            lon0 = sta_lon1[j,0]

            lat1 = sta_lat2[j,0]
            lon1 = sta_lon2[j,0]

            lat2 = sta_lat3[j,0]
            lon2 = sta_lon3[j,0]

            L = strain_tri(lon0,lat0,lon1,lat1,lon2,lat2,n0,e0,n1,e1,n2,e2)
            
            yearstr = "{0:.0f}".format(float(YMAT[i,0]))
            daystr = "{0:.0f}".format(float(DMAT[i,0]))
            timestr = "{0:.4f}".format(float(timer))

            l1str = "{0:.6e}".format(float(L[0]))
            l2str = "{0:.6e}".format(float(L[1]))
            l3str = "{0:.6e}".format(float(L[2]))
            l4str = "{0:.6e}".format(float(L[3]))
            l5str = "{0:.6e}".format(float(L[4]*180/math.pi))
            l6str = "{0:.6e}".format(float(L[5]))
            l7str = "{0:.6e}".format(float(L[6]))
                        
            #outfile = 'output/' + str(k) + '.txt'
            outfile = 'output/' + sta1[k]+ '_' + sta2[k] + '_' + sta3[k] + '.txt'
            fstrain = open(outfile,'a')
            fstrain.write(timestr+' '+yearstr+' '+daystr+' '+l1str+' '+l2str+' '+l3str+' '+l4str+' '+l5str+' '+l6str+' '+l7str+'\n')
            fstrain.close()

        k=k+1

            


        






