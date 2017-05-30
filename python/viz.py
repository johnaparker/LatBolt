import matplotlib.pyplot as plt
from matplotlib import animation
import numpy as np
import h5py

with h5py.File('temp.h5', 'r') as f:
    Ux = f["Ux"]
    Uy = f["Uy"]


    t0 = 0
    tf = Ux.shape[-1]

    
    i = 100
    NX,NY = Ux.shape[0:2]
    u_max  = 0.02
    x = np.arange(NX)
    y = np.arange(NY)
    X,Y = np.meshgrid(x,y, indexing='ij')
    tau = 1
    nu     = (2*tau-1)/6     # kinematic shear viscosity

    umag = np.sqrt( Ux[...,i].T**2 + Uy[...,i].T**2 ) / u_max
    plt.imshow(umag, extent=[0,1,0,1], vmin=0, vmax=1)
    bar = plt.colorbar()
    # stream = plt.streamplot(x,y, Ux[...,i], Uy[...,i], color=[1,1,1])
    skip = 4
    arrow = plt.quiver(x[::skip],y[::skip],Ux[::skip,::skip,i], Uy[::skip,::skip,i], color=[1,1,1])
    plt.xlabel('$x/l_x$')
    plt.ylabel('$y/l_y$')
    bar.set_label('$|\mathbf{u}|/u_\mathrm{max}$')
    td = 1/(nu*(2*np.pi/NX)**2 + (2*np.pi/NY)**2)


    def update(i):
        Ux_p = np.ravel(Ux[:,:,i].T)
        Uy_p = np.ravel(Uy[:,:,i].T)
        im.set_array(Ux_p**2 + Uy_p**2)
        arrow.set_UVC(Ux[::skip,::skip,i], Uy[::skip,::skip,i])
        return im,arrow

    im = plt.pcolormesh(Ux[...,0].T**2 + Uy[...,0].T**2, animated=True,  cmap='viridis')
    ani = animation.FuncAnimation(plt.gcf(), update, np.arange(t0,tf), interval=30, blit=True)

    plt.gca().set_aspect('equal')
    plt.show()

    data = Ux[40,40,:]
    plt.plot(np.abs(data))
    plt.show()
