import numpy as np
import netCDF4 as nc
import matplotlib.pyplot as plt

# input files
file1d = 'GYRE_1d_00010101_00301230_series_energy.nc'
file360d = 'GYRE_360d_00010101_00301230_ocean.nc'
file5d = 'GYRE_5d_00010101_00301230_surf.nc'
filemask = 'mesh_mask.nc'

# ---------------- functions ----------------
def plot_energy(num):
    ds = nc.Dataset(file1d)
    energy = ds['kin_energ']
    time = ds['time_centered']
    dt = time[1] - time[0]
    time = (time - min(time) + dt) / 86400 / 360
    
    plt.subplot(sub_y*100 + sub_x*10 + num)
    plt.plot(time,energy)
    plt.xlabel('years')
    plt.ylabel('$m^2/s^2$')
    plt.title('xyz-averaged resolved kinetic energy')
    plt.ylim(0, 0.01)
    return

def read_average_surf(dataset, variable, start_idx):
    return np.mean(dataset[variable][start_idx:,1,:,:], axis=0)
def read_average_3D(dataset, variable, start_idx):
    return np.mean(dataset[variable][start_idx:,:,:,:], axis=0)

def plot_surf_EKE(num):
    ds = nc.Dataset(file360d)
    u  = read_average_surf(ds, 'vozocrtx', 10)
    v  = read_average_surf(ds, 'vomecrty', 10)
    u2 = read_average_surf(ds, 'squaredu', 10)
    v2 = read_average_surf(ds, 'squaredv', 10)
    lat = ds['nav_lat_grid_T'][:]
    lon = ds['nav_lon_grid_T'][:]
    EKE = 0.5 * (u2+v2-np.square(u)-np.square(v)) # interpolation error neglected
    
    plt.subplot(sub_y*100 + sub_x*10 + num)
    plt.contourf(lon, lat, EKE, 15, vmin=0, vmax=0.3, cmap='ocean')
    plt.colorbar()
    plt.xlabel('Longitude')
    plt.ylabel('Latitude') 
    plt.title('surface EKE, $m^2/s^2$')   
    return

def plot_z_EKE(num):
    ds = nc.Dataset(file360d)
    u  = read_average_3D(ds, 'vozocrtx', 10)
    v  = read_average_3D(ds, 'vomecrty', 10)
    u2 = read_average_3D(ds, 'squaredu', 10)
    v2 = read_average_3D(ds, 'squaredv', 10)
    EKE = 0.5 * (u2+v2-np.square(u)-np.square(v)) # interpolation error neglected
    EKEz = np.mean(EKE[:,1:-1,1:-1], axis=(1,2)).squeeze()
    
    ds = nc.Dataset(filemask)
    depth = ds['gdept_1d'][:].squeeze()

    plt.subplot(sub_y*100 + sub_x*10 + num)

    plt.semilogy(EKEz, depth)
    plt.xlabel('EKE, $m^2/s^2$')
    plt.ylabel('depth, m')
    plt.grid()
    plt.xlim(0 , 0.03)
    plt.ylim(4000, 5)
    plt.title('Lateral mean EKE')   
    return

def plot_SST(num):
    ds = nc.Dataset(file360d)
    T = read_average_surf(ds, 'votemper', 10)
    lat = ds['nav_lat_grid_T'][:]
    lon = ds['nav_lon_grid_T'][:]
    
    plt.subplot(sub_y*100 + sub_x*10 + num)
    plt.contourf(lon, lat, T, np.linspace(8,30,23), vmin=11, vmax=24, cmap='jet')
    plt.colorbar()
    plt.xlabel('Longitude')
    plt.ylabel('Latitude')    
    plt.title('SST, $C^o$')
    return

def plot_snapshot_vorticity(num):
    ds = nc.Dataset(file5d)
    w = ds['surf_rot'][-1-72+18,:,:] # 30 March of the last year
    lat = ds['nav_lat'][:]
    lon = ds['nav_lon'][:]
    
    plt.subplot(sub_y*100 + sub_x*10 + num)
    plt.contourf(lon, lat, w, 30, vmin=-0.5, vmax=0.5, cmap='gray')
    plt.colorbar()
    plt.xlabel('Longitude')
    plt.ylabel('Latitude') 
    plt.title('relative voriticity in Coriolis units')   
    return

def plot_flux(num1, num2):
    rau0 = 1026.
    rcp = 3991.86795711963
    ds = nc.Dataset(file360d)
    u  = read_average_3D(ds, 'vozocrtx', 10)
    v  = read_average_3D(ds, 'vomecrty', 10)
    T  = read_average_3D(ds, 'votemper', 10)
    uT = read_average_3D(ds, 'veltempx', 10)
    vT = read_average_3D(ds, 'veltempy', 10)
    lat = ds['nav_lat_grid_T'][:]
    ds = nc.Dataset(filemask)
    umask = ds['umask'][:].squeeze()
    vmask = ds['vmask'][:].squeeze()
    depth = ds['gdept_1d'][:].squeeze()
    e3t = ds['e3t_1d'][:].squeeze()

    Qx = uT - np.multiply(u, T_to_u(T))
    Qy = vT - np.multiply(v, T_to_v(T))

    Qx = np.multiply(Qx, umask) * rau0 * rcp
    Qy = np.multiply(Qy, vmask) * rau0 * rcp

    (Q, llat) = meridional_flux(Qx, Qy, lat, grid_step = 26500, Nlat = 100)

    Qlat = integrate_flux(Q, e3t)

    Psi_MOC = MOC(u, v, lat, 26500, e3t, Nlat = 100)

    plt.subplot(sub_y*100 + sub_x*10 + num1)
    plt.plot(llat, Qlat)
    plt.xlabel('Latitude')
    plt.ylabel('Heat transport, $W$')
    plt.title('Meridional eddy heat transport')
    plt.grid()
    plt.ylim(-1.2e+14, 1.2e+14)

    plt.subplot(sub_y*100 + sub_x*10 + num2)
    Llat, Depth = np.meshgrid(depth, llat)
    plt.contourf(Depth, Llat, Q, 20, vmin = -2.5e+11, vmax = 2.5e+11, cmap = 'jet')
    plt.xlabel('Latitude')
    plt.ylabel('Depth, m')
    plt.title('Meridional eddy heat transport, $W/m$')
    plt.colorbar()
    plt.ylim(4000,5)
    plt.yscale('log')
    plt.contour(Depth, Llat, Psi_MOC, np.linspace(-5,5,21))
    return

def T_to_u(T):
    (nz, ny, nx) = T.shape
    u = np.zeros((nz, ny, nx))
    for i in range(0,nx-1):
        u[:,:,i] = 0.5 * (T[:,:,i] + T[:,:,i+1])
    return u

def T_to_v(T):
    (nz, ny, nx) = T.shape
    v = np.zeros((nz, ny, nx))
    for j in range(0,ny-1):
        v[:,j,:] = 0.5 * (T[:,j,:] + T[:,j+1,:])
    return v

def meridional_flux(Qx, Qy, lat, grid_step, Nlat):
    (nz, ny, nx) = Qx.shape
    
    divQ = np.zeros((nz, ny, nx))
    for i in range(1,nx):
        for j in range(1,ny):
            divQ[:,j,i] = (Qx[:,j,i] - Qx[:,j,i-1]) / grid_step \
                        + (Qy[:,j,i] - Qy[:,j-1,i]) / grid_step

    llat = np.linspace(np.min(lat),np.max(lat), Nlat)

    Q = np.zeros((Nlat,nz))
    for ilat in range(Nlat):
        for k in range(nz):
            Q[ilat,k] = np.sum(np.multiply(divQ[k,:,:],lat<llat[ilat])) * grid_step ** 2
    return Q, llat

def integrate_flux(Q, e3t):
    ds = nc.Dataset(filemask)
    (Nlat, nz) = Q.shape
    Qlat = np.zeros((Nlat,1)).squeeze()
    for ilat in range(Nlat):
        Qlat[ilat] = np.sum(np.multiply(Q[ilat,:], e3t))
    return Qlat

def MOC(u, v, lat, grid_step, e3t, Nlat):
    (Vint, llat) = meridional_flux(u, v, lat, grid_step, Nlat)
    (nz, ny, nx) = u.shape
    Psi = np.zeros((Nlat, nz))
    
    for k in range(nz-2, -1, -1):
        Psi[:,k] = Psi[:,k+1] + e3t[k] * Vint[:,k]
    
    Psi = Psi / 1.e+6 # to Sverdrups
    return Psi


# ------------------ script -------------------
# number of subplots along x and y
sub_x = 3
sub_y = 3
plt.figure(figsize = (22,18), dpi = 300)

plot_energy(1)
plot_surf_EKE(2)
plot_z_EKE(3)
plot_SST(4)
plot_snapshot_vorticity(5)
plot_flux(6,7)

plt.draw()
plt.savefig('result.pdf')