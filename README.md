# PyTAGS
Python package for the Temporal Analysis of GNSS Strains

Brendan Crowell, University of Washington

Last edited September 18, 2018


This is still very much a work in progress, this will run with the properties in the current pytags.props file. Its currently setup for Python3.


For data, both SOPAC and UNAVCO/PBO time series formats should work. I tested SOPAC format very recently, but not the PBO format. I plan to add UNR and PANGA formats soon. You will need to download time series files and place them in a folder specified in the props file. Th


There are two strain modes, both of which work right now, triangulated and gridded strains. 

There are two temporal modes that work right now, daily time series and static velocity file (in GMT psvelo format). I do intend to add high-rate time series to this (reading from SAC or MSEED format), but right now I haven't coded it up. 


To run static velocity file, change temporal mode to 1 in 'pytags.props'. 
