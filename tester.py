#!/usr/bin/python27
from datetime import datetime, timedelta




mjd_1 = datetime(1858,11,17,0,0,0)+timedelta(days=50001)

tt=mjd_1.timetuple()
print tt.tm_year











