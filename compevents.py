#!/usr/bin/env python
from obspy.core import Stream, read, UTCDateTime
import matplotlib.pyplot as plt
import numpy as np

stas = ['ENG4','TST1']
locs = ['00','00']
chans = ['EH0','BH0']
sensors = ['EpiSensor', 'STS-2']

stF = Stream()

day = 22
year = 2017

stime = UTCDateTime('2017-022T04:38:00.0')


# Why not add an error for no episensor response
pazEpi = {'poles': [-981. + 1009j, -981. - 1009j , -3290.+1263j, -3290.-1263j],'zeros':[0.], 'gain': 2.45956*10**13, 'sensitivity': (2.**26)/40.*9.81}
pazSTS2 = {'gain': 5.96806*10**7, 'zeros': [0, 0], 'poles': [-0.035647 - 0.036879j,  
        -0.035647 + 0.036879j, -251.33, -131.04 - 467.29j, -131.04 + 467.29j],
        'sensitivity': 3.355500*10**10}


for slcs in zip(stas, locs, chans, sensors):
    st = read('/msd/XX_' + slcs[0] + '/' + str(year) + '/' + str(day).zfill(3) + '/' + slcs[1] + '_' + slcs[2] + '*')
    if slcs[3] == 'EpiSensor':
        paz = pazEpi
    elif slcs[3] == 'STS-2':
        paz = pazSTS2
    st.trim(stime-60., stime + 60.*60.*2.)
    st.detrend()
    st.taper(max_percentage=0.05, type='cosine')
    for tr in st:
        tr.simulate(paz_remove=paz)
        tr.filter("bandpass",freqmin = 1./50.,freqmax= 1./20., corners=4.)
    #st.trim(, stime+60.60*2.)
    stF += st


fig = plt.figure(1)
for tr in stF:
    t=np.arange(0, len(tr.data))/tr.stats.sampling_rate
    plt.plot(t,tr.data)
plt.show()


print('Here is our diff:' + str(stF[0].std()/stF[1].std()))
