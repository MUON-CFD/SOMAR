import numpy as np
try:
    from SOMAR import SOMAR
except: 
    pass
def generate_vortex_ring_velocity(points, Gamma, R, a, tilt, smoothing = "rankine", N_theta=200):
    """
    Compute the radial and vertical velocity components (u_r, u_z)
    induced by a Rankine-core vortex ring at given points.

    Parameters:
    - points: array of shape (N, 3), each row is (x, y, z)
    - Gamma: circulation
    - R: ring radius
    - a: core radius (Rankine core radius)
    - N_theta: number of segments to discretize the ring

    Returns:
    - u_x: x velocity at each input point
    - u_y: y velocity at each input point
    - u_z: vertical velocity at each input point
    """
    theta = np.linspace(0, 2 * np.pi, N_theta, endpoint=False)
    # define and rotate the circle around y axis
    x_ring = R * np.cos(theta) * np.cos(tilt)
    y_ring = R * np.sin(theta) 
    z_ring = - R * np.cos(theta) * np.sin(tilt) 
    


    dtheta = 2 * np.pi / N_theta
    u_x = np.zeros(points.shape[0])
    u_y = np.zeros(points.shape[0])
    u_z = np.zeros(points.shape[0])

    x = points[:,0]
    y = points[:,1]
    z = points[:,2]
       

    ux = np.zeros_like(x)
    uy = np.zeros_like(y)
    uz = np.zeros_like(z)
        
    for j in range(N_theta):
        x0 = x_ring[j]
        y0 = y_ring[j]
        z0 = z_ring[j]

        dx = x - x0
        dy = y - y0
        dz = z - z0
        r_sq = dx**2 + dy**2 + dz**2
        delta_sq = a**2*np.ones_like(r_sq)
        if (smoothing == 'rankine'):
            delta_sq[r_sq >= a**2] = 0.0
        
            
        r_smoothed = (r_sq + delta_sq)**1.5


        # Tangent vector at segment
        tx = -np.sin(theta[j]) * R * dtheta * np.cos(tilt)
        ty =  np.cos(theta[j]) * R * dtheta
        tz =  np.sin(theta[j]) * R * dtheta * np.sin(tilt)

        # Cross product s x r
        cx = ty * dz - tz * dy
        cy = tz * dx - tx * dz
        cz = tx * dy - ty * dx

        coeff = Gamma / (4 * np.pi * r_smoothed)

        ux += coeff * cx
        uy += coeff * cy
        uz += coeff * cz
        
        # print("theta", theta[j], "dx", dx, "dy", dy, "dz", dz, "r_sq", r_sq, "delta_sq", delta_sq, "r_smoothed", r_smoothed)
        # print("tx", tx, "ty", ty, "tz", tz)
        # print("cx", cx, "cy", cy, "cz", cz)
        # print("coeff", coeff)
        # print("ux", ux, "uy", uy, "uz", uz)
        # print("********************")
    
    
    
    
    
    return ux, uy, uz

def generate_rankine_buoyancy(coords,R,a,tilt, B=-1):
    """
    return the buoyancy field with B in the core and 0 outside
    """
    import matplotlib.pyplot as plt
    b = np.zeros_like(coords[...,0])
    rotateCoords = np.zeros_like(coords[...])
    rotateCoords[...,0] =  np.cos(tilt) * coords[...,0] - np.sin(tilt) * coords[...,2]
    rotateCoords[...,1] = coords[...,1]
    rotateCoords[...,2] = np.sin(tilt) * coords[...,0] + np.cos(tilt) * coords[...,2]
    r = np.sqrt(rotateCoords[...,0]**2 + rotateCoords[...,1]**2)
    arg = ((r-R)**2 + rotateCoords[...,2]**2)/a**2
    b[arg<1]=B
    
    return b

try:
    @SOMAR
    def initialize(u,v,w,b,Ucoords,Vcoords,Wcoords, Bcoords,R,a,Gamma,alpha):
        dx = Ucoords[1,0,0,0] - Ucoords[0,0,0,0]
        dy = Vcoords[0,1,0,1] - Vcoords[0,1,0,1]
        dz = Wcoords[0,0,1,2] - Wcoords[0,0,1,2]
        # x-direction 
        
        b[...] = np.reshape(generate_rankine_buoyancy(Bcoords, R, a, alpha),b[...].shape) 
        points = np.array([np.ravel(Ucoords[...,0])-dx/2,np.ravel(Ucoords[...,1]), np.ravel(Ucoords[...,2])]).T
        U, _, _ = generate_vortex_ring_velocity(points, Gamma, R, a, alpha)
        u[...] = np.reshape(U,u[...].shape)

        points =  np.array([np.ravel(Vcoords[...,0]),np.ravel(Vcoords[...,1]) - dy/2, np.ravel(Vcoords[...,2])]).T
        _, V, _ = generate_vortex_ring_velocity(points, Gamma, R, a, alpha)
        v[...] = np.reshape(V,v[...].shape)
        
        points =  np.array([np.ravel(Wcoords[...,0]),np.ravel(Wcoords[...,1]), np.ravel(Wcoords[...,2])-dz/2]).T
        _, _, W = generate_vortex_ring_velocity(points, Gamma, R, a, alpha)
        w[...] = np.reshape(W, w[...].shape)
        
        #points = np.array([np.stack(coords[...,0]),np.stack(coords[...,1]), np.stack(coords[...,2])])
except:
    pass
   
if __name__ == "__main__":
    # Example usage
    import numpy as np
    # coords=np.array([1.,0.,1.]).reshape((1,3))
    coords = np.mgrid[-1:1:9j, -1:1:9j, -1:1:9j]
    print("Coordinates shape:", coords.shape)
    print(coords[1,:,4,:])
    coords = np.reshape(coords, (3, -1)).T  # Reshape to (N, 3) where N is the number of points
    
    print("Reshaped coordinates shape:", coords.shape)
    R = 0.5
    a = 0.1
    tilt = np.pi / 4
    B = -1
    
    u,v,w = generate_vortex_ring_velocity(coords, Gamma=1.0, R=R, a=a, tilt=tilt, N_theta=2)
    v = v.reshape((9,9,9))
    print("Buoyancy field shape:", v.shape)
    import matplotlib.pyplot as plt
    plt.imshow(v[:,1, :], extent=(-1, 1, -1, 1), origin='lower')
    plt.show()

