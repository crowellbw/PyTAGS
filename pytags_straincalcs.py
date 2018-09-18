#!/usr/bin/python
import math
import numpy
import obspy.signal
from pytags_coord_tools import ll2utm




def strain(x, y, e, n, props):
        l1 = len(x)
        G = numpy.zeros([2*l1,6])
        U = numpy.zeros([2*l1,1])
        W = numpy.zeros([2*l1,1])
        D = numpy.sqrt(numpy.power(x,2)+numpy.power(y,2))/1000
        gridweight = float(props.getgriddistweight())
        for i in range (0, l1):
                U[2*i,0] = e[i]
                U[2*i+1,0] = n[i]
                w = numpy.exp(-numpy.power(D[i],2)/2/gridweight/gridweight)
                G[2*i,0] = 1*w #Translation terms, east
                G[2*i,1] = 0

                G[2*i+1,0] = 0 #Translation terms, north
                G[2*i+1,1] = 1*w

                G[2*i,2] = x[i]*w #strain terms, east
                G[2*i,3] = y[i]*w

                G[2*i+1,4] = x[i]*w #strain terms, north
                G[2*i+1,5] = y[i]*w

                W[2*i,0] = w

                W[2*i+1,0] = w

                
        L = numpy.linalg.lstsq(G,numpy.multiply(W,U))[0]
        dx = L[0]
        dy = L[1]
        exx = L[2]
        exy = L[3]
        eyx = L[4]
        eyy = L[5]

        Exy = 0.5*(exy+eyx) #shear strain rate
        e1 = (exx+eyy)/2+math.sqrt(math.pow(exx-eyy,2)/4+math.pow(Exy,2))
        e2 = (exx+eyy)/2-math.sqrt(math.pow(exx-eyy,2)/4+math.pow(Exy,2))

        if e1 > e2:
                E1 = e1
                E2 = e2
        else:
                E1 = e2
                E2 = e1
        w = -0.5*(exy-eyx)

        theta = 1/2.0*math.atan2(2.0*Exy,exx-eyy)


        return (dx,dy,E1,E2,Exy,w,theta)

def strain_tri(lon0,lat0,lon1,lat1,lon2,lat2,n0,e0,n1,e1,n2,e2):
        meanlat = (lat0+lat1+lat2)/3.0;
        meanlon = (lon0+lon1+lon2)/3.0;
        (x0,y0) = ll2utm(lon0,lat0, meanlon, meanlat)
        (x1,y1) = ll2utm(lon1,lat1, meanlon, meanlat)
        (x2,y2) = ll2utm(lon2,lat2, meanlon, meanlat)

        dx1 = x1-x0
        dx2 = x2-x0

        dy1 = y1-y0
        dy2 = y2-y0

        dn1 = (n1-n0)
        dn2 = (n2-n0)

        de1 = (e1-e0)
        de2 = (e2-e0)

        X = numpy.array([[dx1,dy1],[dx2,dy2]])
        detX = numpy.linalg.det(X)

        w = ((dy2*dn1-dy1*dn2)-(dx1*de2-dx2*de1))/2/detX

        Eee = (dy2*de1-dy1*de2)/detX

        Enn = (dx1*dn2-dx2*dn1)/detX

        Een = (dx1*de2-dx2*de1)/detX+w

        theta = 1/2.0*math.atan2(2.0*Een,Enn-Eee)


        e1 = Eee*math.cos(theta)*math.cos(theta)+Enn*math.sin(theta)*math.sin(theta)-2*Een*math.sin(theta)*math.cos(theta)
        e2 = Eee*math.sin(theta)*math.sin(theta)+Enn*math.cos(theta)*math.cos(theta)+2*Een*math.sin(theta)*math.cos(theta)

        if (e2 > e1):
                e1new = e2
                e2new = e1
                e1 = e1new
                e2 = e2new
        if (theta < 0):
                theta = math.pi+theta

        return (e1,e2,Een,theta,w)

