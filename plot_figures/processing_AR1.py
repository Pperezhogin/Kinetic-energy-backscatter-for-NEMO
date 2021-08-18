import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt

# input files
file1d = 'GYRE_1d_00010101_00301230_series_KEB.nc'
file360d = 'GYRE_360d_00010101_00301230_KEB.nc'
filemask = 'mesh_mask.nc'

# ---------------- functions ----------------
def plot_series(variable, ylabel, title, limits, num):
    ds = nc.Dataset(file1d)
    var = ds[variable]
    time = ds['time_centered']
    dt = time[1] - time[0]
    time = (time - min(time) + dt) / 86400 / 360
    
    plt.subplot(sub_y*100 + sub_x*10 + num)
    plt.plot(time,var)
    plt.xlabel('years')
    plt.ylabel(ylabel)
    plt.title(title)
    plt.ylim(limits)
    plt.grid()
    return

def plot_2dfield(variable, title, num):
    ds = nc.Dataset(file360d)
    field = read_average_2D(ds, variable, 10)
    lat = ds['nav_lat'][:]
    lon = ds['nav_lon'][:]
    
    plt.subplot(sub_y*100 + sub_x*10 + num)
    plt.contourf(lon, lat, field, 20, cmap='ocean')
    plt.colorbar()
    plt.xlabel('Longitude')
    plt.ylabel('Latitude') 
    plt.title(title)   
    return

def plot_average_xy(variable, xlabel, title, limits, num):
    ds = nc.Dataset(file360d)
    var = read_average_3D(ds, variable, 10)
    ds = nc.Dataset(filemask)
    depth = ds['gdept_1d'][:].squeeze()

    var_z = np.mean(var[:,1:-1,1:-1], axis=(1,2))
    
    plt.subplot(sub_y*100 + sub_x*10 + num)
    plt.semilogy(var_z,depth)
    plt.ylim(4000,5)
    plt.xlim(limits)
    plt.xlabel(xlabel)
    plt.ylabel('Depth, m') 
    plt.title(title)   
    plt.grid()
    return

def read_average_surf(dataset, variable, start_idx):
    return np.mean(dataset[variable][start_idx:,1,:,:], axis=0)
def read_average_2D(dataset, variable, start_idx):
    return np.mean(dataset[variable][start_idx:,:,:], axis=0)
def read_average_3D(dataset, variable, start_idx):
    return np.mean(dataset[variable][start_idx:,:,:,:], axis=0)

# ------------------ script -------------------
# number of subplots along x and y
sub_x = 3
sub_y = 3
plt.figure(figsize = (22,18), dpi = 300)

plot_2dfield('ediss_z', 'z-averaged $E_{diss}$, $m^2/s^3$', 1)
plot_2dfield('eback_z', 'z-averaged $E_{back}$, $m^2/s^3$', 2)
plot_average_xy('eback', '$m^2/s^3$', 'xy-averaged $E_{back}$', (0, 5e-8), 3)
plot_average_xy('ediss', '$m^2/s^3$', 'xy-averaged $E_{diss}$', (0, 5e-8), 4)
plot_series('avr_cdiss', '', 'average $c_{diss}$, i.e. $\overline{E_{back}} / \overline{E_{diss}}$', (0, 1), 5)
plot_series('SKEB_amp', 'amplitude factor', 'a posteriori correction', (0.5, 1.5), 6)


plt.draw()
plt.savefig('result_KEB.pdf')