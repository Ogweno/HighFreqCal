#!/usr/bin/env python

# Import wild cards
import glob

import numpy as np

# Import the read function from obspy
from obspy.core import read, Stream, UTCDateTime

stime = UTCDateTime('2016-288T15:05:48')
etime = stime + 50.

# We will make a debug flag
debug = True

files = glob.glob('/msd/IU_HRV/2016/288/00_BH*.seed')
if debug:
    print(files)

# Make a stream object to hold data
st = Stream()
for curfile in files:
    if debug:
        print('Here is our current file: ' + curfile)
    st += read(curfile, starttime=stime, endtime=etime)
st.merge()
if debug:
    print(st)
    st.plot()
    
# Here we can do a 10 volt test
st2 = st.copy()
nstime = UTCDateTime('2016-288T15:05:53.0')
netime = nstime + 5.
st2.trim(starttime=nstime, endtime=netime)
# st2 Now has data for our first 10 V test
# Now we convert st2 into Volts

# Make two empty lists to save our results
mminus =[]
sminus =[]

# This is the first voltage 
for tr in st2:
    tr.data.astype(np.float64)
    # Here we are converting from counts to Volts
    tr.data = tr.data*(40. /(2.**26))
    print('Here is the mean voltage: ' + str(np.mean(tr.data)) + ' for: ' + str(tr.id))
    print('Here is the standard deviation:' + str(np.std(tr.data)) + ' for: ' + str(tr.id))
    mminus.append(np.mean(tr.data))
    sminus.append(np.std(tr.data))
    
# Now we can check the second voltage
st3 = st.copy()
nstime = UTCDateTime('2016-288T15:06:10.0')
netime = nstime + 15.
st3.trim(starttime=nstime, endtime=netime)

mplus = []
splus =[]

for tr in st3:
    tr.data.astype(np.float64)
    # Here we are converting from counts to Volts
    tr.data = tr.data*(40. /(2.**26))
    print('Here is the mean voltage: ' + str(np.mean(tr.data)) + ' for: ' + str(tr.id))
    print('Here is the standard deviation:' + str(np.std(tr.data)) + ' for: ' + str(tr.id))
    mplus.append(np.mean(tr.data))
    splus.append(np.std(tr.data))

# Now we have our voltages and stds.
# We will figure out the actual bit-weight since the above was a sanity check

# These are the p2pM and p2pStd
p2pM = []
p2pStd =[]
ids =[]

# Notice these will be in counts/V
for trs in zip(st2,st3):
    trs[0].data = trs[0].data.astype(np.float64)*(2.**26)/40.
    trs[1].data = trs[1].data.astype(np.float64)*(2.**26)/40.
    p2pM.append(abs(np.mean(trs[0].data)/20.) + abs(np.mean(trs[1].data)/20.))
    p2pStd.append(((np.std(trs[0].data)/20.)**2 + np.std((trs[1].data)/20.)**2)**.5)
    ids.append(trs[0].id)
    
# Here is our final calibration
for val in zip(p2pM, p2pStd,ids):
    print('Final sensitivity of ' + val[2] + ' is ' + str(val[0]) + '+/-' + str(2.5*val[1]))
    # Here I am doing a percentage deviation (I bet this is a bit poorly done)
    print('Deviation from nominal=' + str((((val[0]/(2.**26/40.)))-1.)*100.) + ' %')
    










