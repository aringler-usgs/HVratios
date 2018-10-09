#!/usr/bin/env python
from obspy.clients.fdsn import Client
from obspy import read_inventory
from obspy.core import UTCDateTime 
from obspy.imaging.maps import plot_basemap
import matplotlib.pyplot as plt

# If you are going to play with this with US or something
# change resolution to "c"

from matplotlib.cm import get_cmap
codes = set([0,1])
cmap = get_cmap(name="viridis", lut=len(codes))


client = Client("IRIS")
nets = ['IC',"IU"]
#nets = ['IC']
stations = []
for idx, net in enumerate(nets):
    # Grab the list of station
    if 'inv' not in vars():
        inv = read_inventory('/APPS/metadata/SEED/' + net + '.dataless', format='SEED')
    else:
        inv += read_inventory('/APPS/metadata/SEED/' + net + '.dataless', format='SEED')
    #stations = getstalist(sp, stime, net)
    #if idx > 0:
    #fig = inv.plot(method='basemap')
    fig = inv.plot(method='basemap', show=False, label=False, size=32., color='blue', continent_fill_color='0.9', resolution='h', color_per_network=True, colormap=cmap)
for idx, net in enumerate(nets):    
    ax = fig.axes[0]
    ax.scatter([],[], 28, color=cmap(idx), label=net, marker="v")
stime = UTCDateTime('2001-001T00:00:00.0')
etime = UTCDateTime('2018-001T00:00:00.0')


cat = client.get_events(starttime=stime, endtime=etime, 
                                minmagnitude=6.5, maxdepth=50.)


#cat.plot(method='basemap', fig=fig, title=None)

lats = []
lons = []
labels = []
mags = []
colors = []
times = []

for event in cat:
    origin = event.preferred_origin() or event.origins[0]
    lats.append(origin.latitude)
    lons.append(origin.longitude)
    times.append(origin.time)
    magnitude = event.preferred_magnitude() or event.magnitudes[0]
    mag = magnitude.mag
    mags.append(mag)
    colors.append('k')



#fig = plot_basemap(lons, lats, 22., 'g', show=False, fig=fig, continent_fill_color='0.5', resolution='i')
fig = plot_basemap(lons, lats, 22., 'g', show=False, fig=fig, continent_fill_color='0.9', resolution='h')
ax = fig.axes[0]
ax.scatter([],[], 22., 'g', label='Earthquake', marker="o")
plt.legend(loc=9, ncol=3, bbox_to_anchor=(0.5,-0.01), scatterpoints=1)
fig.savefig('Map.pdf',format='PDF', dpi=600)
fig.savefig('Map.png',format='png', dpi=600)
