"""

Generate the following output:
- a synthetic earthquake catalogue (catdata.txt)
- the distribution of effective asperities (active_patches.txt)
- visualise fault slip for the 500-1000 years time window (slip_plot.pdf)

No input is required, as model settings are read from visc.py

Contributions from Adam Beall and Martijn van den Ende


"""


import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
import matplotlib as mpl
import numpy as np
from lib import catalogue
import visc


# Load model
params,p = visc.load()
p.read_output()

t_yr =  3600 * 24 * 365.0
W = params['W']


# Generate earthquake catalogue

qcat = catalogue.generate_catalogue(p)
qcat.calc_magnitudes()
arrCatData = np.vstack([qcat.t/t_yr,qcat.Mw,qcat.lengths,qcat.duration]).T
np.savetxt('catdata.txt',arrCatData)

# slice and reshape data

start_time = 500 * t_yr
end_time = params['maxtime']
end_time = 1000 * t_yr

warm_up =  start_time 

fx = p.ox['t'][:]
x_unique = p.ox["x"].unique()
sort_inds = np.argsort(x_unique)
x_unique = x_unique[sort_inds]
t_vals = np.sort(p.ox["t"].unique())

ind_warmup = np.where(t_vals >= warm_up)[0]
ind_end = np.where(t_vals >= end_time)[0]
if len(ind_end) == 0:
    ind_end = t_vals.size
else:
    ind_end = ind_end[0]

if len(ind_warmup) == 0:
    if t_vals.size > 100:
        ind_warmup = 100
    else:
        ind_warmup = 0
else:
    ind_warmup = ind_warmup[0]

Nx = len(x_unique)
Nt = len(t_vals) - 1

# slices = np.s_[Nx * ind_warmup:Nx * Nt]
data_shape = (ind_end - ind_warmup, Nx)

x = p.ox["x"][Nx * ind_warmup:Nx * ind_end].values.reshape(data_shape)[:, sort_inds]
time = p.ox["t"][Nx * ind_warmup:Nx * ind_end].values.reshape(data_shape)[:, sort_inds]
v = p.ox["v"][Nx * ind_warmup:Nx * ind_end].values.reshape(data_shape)[:, sort_inds]
# slip = p.ox["slip"][Nx * ind_warmup:Nx * ind_end].values.reshape(data_shape)[:, sort_inds]
slip = p.ox["slip"][Nx * ind_warmup:Nx * ind_end].values.reshape(data_shape)[:, sort_inds]
vel = p.ox["v"][Nx * ind_warmup:Nx * ind_end].values.reshape(data_shape)[:, sort_inds]

x = x[0,:]
time = time[:,0] / t_yr

slip0 = np.array(slip[0,:])

for i in range(time.size):
    slip[i,:] -= slip0 


# Calculate steadystate viscofrictional stress, to plot effective asperities
min_length = 315.68
eta_min = 18.0
eta_max = 20.0

etadata = np.loadtxt('eta_dist_%.2f.txt' %W)
etadata[:,0] -= np.average(etadata[:,0])


# Calculate effective asperities

arrL = []
etathresh = 60.0e6 * W / 1e-9
stress_drop = 0e6 * W / 1e-9
print("%.2e" %etathresh)
L = 0.0

arrActivepatchesX = [] 
arrActivepatchesY = [] 
prevpatch = False
includepatch = False

arrTempX = []
arrTempY = []

for i in range(1,etadata[:,1].size):
    eta = etadata[i,1]
    prevpatch = includepatch

    includepatch = eta > (etathresh - stress_drop)

    if includepatch:
        if not prevpatch:
            arrTempX.append(etadata[i,0])
            arrTempY.append(0.0)

        L += etadata[i,0] - etadata[i-1,0]
        arrTempX.append(etadata[i,0])
        arrTempY.append(1.0)

    if (not includepatch and L>0.0) or i == etadata[:,1].size -1:
        if L*1e3 > min_length:
            arrL.append(L*1e3)
            
            arrActivepatchesX = arrActivepatchesX+arrTempX
            arrActivepatchesY = arrActivepatchesY+arrTempY
            
            arrActivepatchesX.append(etadata[i,0])
            arrActivepatchesX.append(etadata[i,0])
            arrActivepatchesY.append(1.0)
            arrActivepatchesY.append(0.0)

        arrTempX = []
        arrTempY = []

        L = 0.0

np.savetxt('patch_sizes.txt',arrL)
np.savetxt('active_patches.txt',np.vstack([arrActivepatchesX,arrActivepatchesY]).T)






# Make slip plot

w, h = plt.figaspect(0.5)
plt.figure(figsize=(1.5*w,1.5*h))


velA = np.ma.masked_where(vel > 1e-6, vel)
velB = np.ma.masked_where(vel < 1e-6, vel)



#pcm = plt.pcolormesh((x[::1]+15e3)/1e3, slip[::1], np.log10(vel[::1,::1]),vmin=-11,vmax=1,rasterized=True)#,linewidths=1e-3,edgecolor='black')


pcm = plt.pcolormesh((x+15e3)/1e3, slip, np.log10(velB),vmin=-11,vmax=1,zorder = 1,rasterized=True)#,linewidths=1e-3,edgecolor='black')


catdata = np.loadtxt('catdata.txt')
arrStart = []
arrEnd = []
for i in range(len(catdata[:,0])):
    event_time = catdata[i,0] 
    if event_time > start_time/t_yr and event_time < end_time/t_yr:
        event_end_time = event_time + catdata[i,3]/60.0/24.0/365.0

        arrStart.append( event_time)
        arrEnd.append( event_end_time )


for i in range(len(arrStart)):
    arrEventTime = np.arange(arrStart[i],arrEnd[i],2. / 3600.0 / 24.0 / 365.0)
    arrI = np.zeros_like(arrEventTime)
    for j in range(len(arrEventTime)):
        arrI[j] = np.argmin( np.abs(arrEventTime[j] - time))



    for idx in arrI:
        plt.plot((x+15e3)/1e3,slip[int(idx),:],c="black",lw=0.25,zorder = 2)

pcm = plt.pcolormesh((x+15e3)/1e3, slip, np.log10(velA),vmin=-11,vmax=1,zorder=3,rasterized=True)#,linewidths=1e-3,edgecolor='black')



dt = 20.0
arrTime = np.arange(start_time,end_time,dt*t_yr)
arrI = np.zeros_like(arrTime)

for i,t in enumerate(arrTime):
    arrI[i] = np.argmin( np.abs(time*t_yr-t))
    print(arrI[i])


for i in arrI:
    plt.plot((x+15e3)/1e3,slip[int(i),:],c="white",lw=0.25,zorder = 4)



    

Blues = plt.get_cmap('Blues')

etaCap = np.minimum(etadata[:,1],etathresh)



pX, pY = etadata[0,:]
rect_scale = 2.0

for k in range(1,len(etadata[:,0])):
    if abs(etadata[k,1] - pY) > 1e-3:
        dx = etadata[k-1,0]-pX
        xbox = pX
        eta_val = np.log10(pY)
        eta_norm = (eta_val - eta_min) / (eta_max-eta_min)
        plt.gca().add_patch(Rectangle(((xbox*1e3+15e3)/1e3,-1.0 * rect_scale),dx,0.5*rect_scale,color=Blues(eta_norm)))


        pX = etadata[k,0]
        pY = etadata[k,1]


pX = etadata[0,0]
pY = etaCap[0]



for k in range(1,len(etadata[:,0])):
    if abs(etaCap[k] - pY) > 1e-3:
        dx = 0.5*(etadata[k-2,0]+ etadata[k-1,0])-pX
        xbox = pX
        eta_val = np.log10(pY)
        eta_norm = (eta_val - eta_min) / (eta_max-eta_min)

        pX = etadata[k,0]
        pY = etaCap[k]

        if eta_val == np.log10(etathresh) and dx*1e3 > min_length:
            plt.gca().add_patch(Rectangle(((xbox*1e3+15e3)/1e3,-0.5*rect_scale),dx,0.5*rect_scale,color='Black'))



plt.xlabel('Distance along fault (km)')
plt.ylabel('Accumulated slip (m)')

plt.ylim(-1*rect_scale,max(slip[-1,:]))
plt.xlim(0,30)
norm = mpl.colors.Normalize(vmin=eta_min, vmax=eta_max)
mappable = mpl.cm.ScalarMappable( norm=norm,cmap=Blues)
# annoying hack for old matplotlib version
mappable._A = []

plt.colorbar(mappable,ax=plt.gca(),label='Log(Viscosity)',orientation='horizontal',shrink=0.4,pad=0.1)
plt.colorbar(pcm,label='Slip velocity (m/s)',orientation='horizontal',shrink=0.4,pad=0.2)

plt.savefig('slip_plot.pdf' ,bbox_inches='tight',dpi=300)


