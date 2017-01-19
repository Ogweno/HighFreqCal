#!/usr/bin/env python

# Import wild cards
import glob

import numpy as np

# Import the read function from obspy
from obspy.core import read, Stream, UTCDateTime

stime = UTCDateTime('2017-011T16:58:52')
etime = stime + 700.

# We will make a debug flag
debug = True

files = glob.glob('/tr1/telemetry_days/XX_FBA2/2017/2017_011/00_EH*.seed')
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
nstime = UTCDateTime('2017-011T16:59:55.0')
netime = nstime + 20.
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
nstime = UTCDateTime('2017-011T17:01:00.0')
netime = nstime + 20.
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
	
# Since we converted to Volts in our sanity check above, we need to get back to counts	
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
    

# Ok, Now we move onto our test to determine sensitivity of the accelerometer 

# Start with the vertical 

print('')

stZ1 = st.copy()
stZ1 = stZ1.select(channel="EH0")

nstime = UTCDateTime('2017-011T17:04:20.0')
netime = nstime + 20.
stZ1.trim(starttime=nstime, endtime=netime)
# stZ1 Now has data for our first flip test
# Now we convert stZ1 into acceleration

# Make two empty lists to save our results
mplus =[]
splus =[]

# This is the first voltage 
for tr in stZ1:
    tr.data.astype(np.float64)
    # Here we are converting from counts to Volts and then from volts to acceleration 
    tr.data = tr.data*(40. /(2.**26))*((4.0*9.80665) /20.)
    print('Here is the mean acceleration: ' + str(np.mean(tr.data)) + ' for: ' + str(tr.id))
    print('Here is the standard deviation:' + str(np.std(tr.data)) + ' for: ' + str(tr.id))
    mplus.append(np.mean(tr.data))
    splus.append(np.std(tr.data))
    
print('')
    
# Now we can check the second voltage 
stZ2 = st.copy()
stZ2 = stZ2.select(channel="EH0")
print stZ2
nstime = UTCDateTime('2017-011T17:09:10.0')
netime = nstime + 20.
stZ2.trim(starttime=nstime, endtime=netime)

mminus = []
sminus =[]

for tr in stZ2:
    tr.data.astype(np.float64)
    # Here we are converting from counts to Volts and then from volts to acceleration 
    tr.data = tr.data*(40. /(2.**26))*((4.*9.80665) /20.)
    print('Here is the mean acceleration: ' + str(np.mean(tr.data)) + ' for: ' + str(tr.id))
    print('Here is the standard deviation:' + str(np.std(tr.data)) + ' for: ' + str(tr.id))
    mminus.append(np.mean(tr.data))
    sminus.append(np.std(tr.data))



# Now do the X-direction 

print('')

stY1 = st.copy()
stY1 = stY1.select(channel="EH1")

nstime = UTCDateTime('2017-011T17:05:40.0')
netime = nstime + 20.
stY1.trim(starttime=nstime, endtime=netime)
# stX1 Now has data for our first flip test
# Now we convert stX1 into acceleration

# Make two empty lists to save our results
mplus=[]
splus =[]

# This is the first voltage 
for tr in stY1:
    tr.data.astype(np.float64)
    # Here we are converting from counts to Volts and then from volts to acceleration 
    tr.data = tr.data*(40. /(2.**26))*((4.0*9.80665) /20.)
    print('Here is the mean acceleration: ' + str(np.mean(tr.data)) + ' for: ' + str(tr.id))
    print('Here is the standard deviation:' + str(np.std(tr.data)) + ' for: ' + str(tr.id))
    mplus.append(np.mean(tr.data))
    splus.append(np.std(tr.data))
    
print('')
    
# Now we can check the second voltage 
stY2 = st.copy()
stY2 = stY2.select(channel="EH1")

nstime = UTCDateTime('2017-011T17:07:45.0')
netime = nstime + 20.
stY2.trim(starttime=nstime, endtime=netime)

mminus = []
sminus =[]

for tr in stY2:
    tr.data.astype(np.float64)
    # Here we are converting from counts to Volts and then from volts to acceleration 
    tr.data = tr.data*(40. /(2.**26))*((4.0*9.80665) /20.)
    print('Here is the mean acceleration: ' + str(np.mean(tr.data)) + ' for: ' + str(tr.id))
    print('Here is the standard deviation:' + str(np.std(tr.data)) + ' for: ' + str(tr.id))
    mminus.append(np.mean(tr.data))
    sminus.append(np.std(tr.data))



# Get the final sensitivities for the Z-Axis
# These are the p2pM and p2pStd
p2pM = []
p2pStd =[]
ids =[]

# Notice these will be in V/(m/s^2)
for trs in zip(stZ1,stZ2):
	# Here, We just convert back to volts
    trs[0].data = trs[0].data.astype(np.float64)*(20./(4.0*9.80665))
    print trs[0].data
    trs[1].data = trs[1].data.astype(np.float64)*(20./(4.0*9.80665))
    p2pM.append(abs(np.mean(trs[0].data)/(9.79188087*2.)) + abs(np.mean(trs[1].data)/(9.79188087*2.)))
    p2pStd.append(((np.std(trs[0].data)/(9.79188087*2.))**2 + np.std((trs[1].data)/(9.79188087*2.))**2)**.5)
    ids.append(trs[0].id)
    
# Here is our final calibration for the Z axis
for val in zip(p2pM, p2pStd,ids):
    print('Final sensitivity of ' + val[2] + ' is ' + str(val[0]) + '+/-' + str(2.5*val[1]))
    # Here I am doing a percentage deviation (I bet this is a bit poorly done)
    print('Deviation from nominal=' + str((((val[0]/(20./(4.*9.80665))))-1.)*100.) + ' %')
    
    
    
# Repeat for the Y-Axis    
print('')

p2pM = []
p2pStd =[]
ids =[]

# Notice these will be in V/(m/s^2)
for trs in zip(stY1,stY2):
	# Here, We just convert back to volts
    trs[0].data = trs[0].data.astype(np.float64)*(20./(4.0*9.80665))
    trs[1].data = trs[1].data.astype(np.float64)*(20./(4.0*9.80665))
    p2pM.append(abs(np.mean(trs[0].data)/(9.79188087*2.)) + abs(np.mean(trs[1].data)/(9.79188087*2.)))
    p2pStd.append(((np.std(trs[0].data)/(9.79188087*2.))**2 + np.std((trs[1].data)/(9.79188087*2.))**2)**.5)
    ids.append(trs[0].id)
    
# Here is our final calibration for the Z axis
for val in zip(p2pM, p2pStd,ids):
    print('Final sensitivity of ' + val[2] + ' is ' + str(val[0]) + '+/-' + str(2.5*val[1]))
    # Here I am doing a percentage deviation (I bet this is a bit poorly done)
    print('Deviation from nominal=' + str((((val[0]/(20./(4.*9.80665))))-1.)*100.) + ' %')
       
    

