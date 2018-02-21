#!/usr/bin/env python

from obspy.clients.fdsn import Client
from obspy import UTCDateTime
start_1="2018-01-17 00:00:00"
starttime = UTCDateTime(start_1)
endtime = starttime+600
client = Client(base_url='https://fdsnws.raspberryshakedata.com/')
waveform= client.get_waveforms('AM', 'R0D9C', '00', 'SHZ', starttime, endtime)

waveform.write(waveform[0].id + '.mssed', format='MSEED')
