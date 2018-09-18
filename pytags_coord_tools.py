#!/usr/bin/python
import math
import numpy

#This converts WGS84 coordinates to ITRF xyz

def lla2ecef(lat,lon,alt):
	lat = lat*math.pi/180
	lon = lon*math.pi/180
	a = 6378137
	e = 8.1819190842622e-2

	N = a/numpy.sqrt(1-numpy.power(e,2)*numpy.power(numpy.sin(lat),2))

	x = (N+alt)*numpy.cos(lat)*numpy.cos(lon)
	y = (N+alt)*numpy.cos(lat)*numpy.sin(lon)
	z = ((1-numpy.power(e,2))*N+alt)*numpy.sin(lat)
	
	return (x, y, z)

#This converts ITRF xyz coordinates to WGS84

def ecef2lla(x,y,z):
	a = 6378137
	e = 8.1819190842622e-2

	b = numpy.sqrt(numpy.power(a,2)*(1-numpy.power(e,2)))
	ep = numpy.sqrt((numpy.power(a,2)-numpy.power(b,2))/numpy.power(b,2))
	p = numpy.sqrt(numpy.power(x,2)+numpy.power(y,2))
	th = numpy.atan2(a*z,b*p)
	lon = numpy.atan2(y,x)
	lat = numpy.atan2((z+numpy.power(ep,2)*b*numpy.power(numpy.sin(th),3)),(p-numpy.power(e,2)*a*numpy.power(numpy.cos(th),3)))
	N = a/numpy.sqrt(1-numpy.power(e,2)*numpy.power(numpy.sin(lat),2))
	alt = p/numpy.cos(lat)-N

	lon = lon*180/math.pi
	lat = lat*180/math.pi

	return (lat,lon,alt)

#This takes displacements in x, y, z and converts them to north, east up

def dxyz2dneu(dx,dy,dz,lat,lon):
	lat = lat*math.pi/180
	lon = lon*math.pi/180
	dn = -numpy.sin(lat)*numpy.cos(lon)*dx-numpy.sin(lat)*numpy.sin(lon)*dy+numpy.cos(lat)*dz
	de = -numpy.sin(lon)*dx+numpy.cos(lon)*dy
	du = numpy.cos(lat)*numpy.cos(lon)*dx+numpy.cos(lat)*numpy.sin(lon)*dy+numpy.sin(lat)*dz
	return (dn, de, du)

#This takes lat and lon values and converts to UTM northing and easting for the Pacific Northwest. 
#The central meridian is fixed by lon0 (no lookup table included), so dont use if you are covering a large geographic area (probably greater than 3 UTM zones)

def ll2utm(lon,lat,lon0,lat0):
	#WGS84 parameters
        a = 6378137.0000
        esq = 0.006694380069978522
        epsq = esq/(1-esq)
        k0 = 0.9996
        if lat0 < 0.0:
                Noff = 10000000.0
        else:
                Noff = 0.0



        lon = lon*math.pi/180
        lat = lat*math.pi/180
        lon0 = lon0*math.pi/180 

        A = (lon-lon0)*numpy.cos(lat)
        v = a/numpy.sqrt(1-esq*numpy.power(numpy.sin(lat),2))
        T = numpy.power(numpy.tan(lat),2)
        C = esq*numpy.power(numpy.cos(lat),2)/(1-esq)

        M = a*(
                (1-esq/4-3*numpy.power(esq,2)/64-5*numpy.power(esq,3)/256)*lat
                -(3*esq/8+3*numpy.power(esq,2)/32+45*numpy.power(esq,3)/1024)*numpy.sin(2*lat)
                +(15*numpy.power(esq,2)/256+45*numpy.power(esq,3)/1024)*numpy.sin(4*lat)
                -(35*numpy.power(esq,3)/3072)*numpy.sin(6*lat)
        )

        UTMNorthing = k0*( M + v*numpy.tan(lat)*( A*A/2 + (5-T+9*C+4*C*C)*numpy.power(A,4)/24 + (61-58*T+T*T+600*C-330*epsq)*numpy.power(A,6)/720 ) ) +Noff

        UTMEasting = k0*v*( A + (1-T+C)*numpy.power(A,3)/6 + (5-18*T+T*T+72*C-58*epsq)*numpy.power(A,5)/120 )  + 500000.0

        return (UTMEasting, UTMNorthing)


#Takes UTM easting and northing with a central meridian and converts to lat and lon. Obviously running ll2utm then putting output into utm2ll should result in the original value (off by a small amount due to truncation of the UTM parameters)
def utm2ll(UTMEasting,UTMNorthing,lon0,NorS):
        #WGS84 parameters
        a = 6378137.0000
        esq = 0.006694380069978522
        epsq = esq/(1-esq)
        k0 = 0.9996
        lon0 = lon0*math.pi/180


        if NorS > 0.0:
                Noff = 0.0
        else:
                Noff = 10000000.0

        M1 = (UTMNorthing-Noff)/k0
        mu1 = M1/(a * (1 - esq/4 - 3/64*numpy.power(esq,2) - 5/256*numpy.power(esq,3)) )
        e1 = ( 1 - numpy.sqrt(1-esq) ) / ( 1 + numpy.sqrt(1-esq) )
        lat1 = mu1 + (3*e1/2 - 27/32*numpy.power(e1,3))*numpy.sin(2*mu1) + (21/16*numpy.power(e1,2) - 55/32*numpy.power(e1,4))*numpy.sin(4*mu1) + (151/96*numpy.power(e1,3))*numpy.sin(6*mu1) + (1097/512*numpy.power(e1,4))*numpy.sin(8*mu1)
        T1 = numpy.power(numpy.tan(lat1),2)
        C1 = epsq*numpy.power(numpy.cos(lat1),2)

        v1 = a/numpy.sqrt(1 - esq*numpy.power(numpy.sin(lat1),2))
        p1 = a*(1-esq)/numpy.power(1 - esq*numpy.power(numpy.sin(lat1),2),1.5)
        D = (UTMEasting - 500000.0)/v1/k0

        lat = lat1 - (v1*numpy.tan(lat1)/p1) * ( numpy.power(D,2)/2 - (5 + 3*T1 + 10*C1 - 4*numpy.power(C1,2) - 9*epsq)*numpy.power(D,4)/24 + (61 + 90*T1 +298*C1 + 45*numpy.power(T1,2) - 252*epsq - 3*numpy.power(C1,2))*numpy.power(D,6)/720 ) 
        lon = lon0 + ( D - (1 + 2*T1 + C1)*numpy.power(D,3)/6 + (5 - 2*C1 + 28*T1 - 3*numpy.power(C1,2) + 8*epsq + 24*numpy.power(T1,2))*numpy.power(D,5)/120 )/numpy.cos(lat1)

        lat = lat*180/math.pi
        lon = lon*180/math.pi

        return(lon, lat)



