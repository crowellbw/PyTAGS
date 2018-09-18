#!/usr/bin/python

class Properties:
        def __init__(self,propfilename):
                self.dict={}
                infile=open(propfilename,'r')
                for line in infile:
                        if (line) and (line[0]!='#') and ('=' in line):
                                (key,val)=line.split('=')
                                self.dict[key]=val.strip()
                return

        def getstrainmode(self):
                if 'strainmode' in self.dict:
                        return self.dict['strainmode']
                return ''
        def gettemporalmode(self):
                if 'temporalmode' in self.dict:
                        return self.dict['temporalmode']
                return ''

        def getsitefile(self): 
                if 'sitefile' in self.dict:
                        return self.dict['sitefile']
                return 0

        def getdatadir(self): 
                if 'datadir' in self.dict:
                        return self.dict['datadir']
                return ''

        def getstartyear(self): 
                if 'startyear' in self.dict:
                        return self.dict['startyear']

        def getendyear(self): 
                if 'endyear' in self.dict:
                        return self.dict['endyear']

        def getmaxlon(self): 
                if 'maxlon' in self.dict:
                        return self.dict['maxlon']
        def getminlon(self): 
                if 'minlon' in self.dict:
                        return self.dict['minlon']
        def getmaxlat(self): 
                if 'maxlon' in self.dict:
                        return self.dict['maxlat']
        def getminlat(self): 
                if 'minlat' in self.dict:
                        return self.dict['minlat']
        def getgriddistweight(self): 
                if 'griddistweight' in self.dict:
                        return self.dict['griddistweight']

        def getnlon(self): 
                if 'nlon' in self.dict:
                        return self.dict['nlon']
        def getnlat(self): 
                if 'nlat' in self.dict:
                        return self.dict['nlat']
        def getdatadir(self): 
                if 'datadir' in self.dict:
                        return self.dict['datadir']
        def getdatasource(self): 
                if 'datasource' in self.dict:
                        return self.dict['datasource']
        def getdatasuffix(self): 
                if 'datasuffix' in self.dict:
                        return self.dict['datasuffix']
        def getmaxtrilength(self): 
                if 'maxtrilength' in self.dict:
                        return self.dict['maxtrilength']


        def getstationlength(self):
                if 'sitefile' in self.dict:
                        sf = self.dict['sitefile']
                        k=0
                        with open(sf) as f:
                                for i, l in enumerate(f):
                                     if not l.startswith("#"):
                                        k=k+1
                        stationlength = k
                        return stationlength





