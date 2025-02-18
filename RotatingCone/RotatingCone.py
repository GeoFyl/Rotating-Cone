
# Calculate moments of inertia for cone from mass, height and base radius 
def GetMomentsOfInertia(M, h, r):
    I1 = I2 = 3 * M * (r*r + (h*h)/4) / 20
    I3 = 6 * M * r*r / 20
    
    return I1, I2, I3

# Performs Runge-Kutta computation
def RungeKutta(Wx0, Wy0, Wz0, I1, I2, I3, stepSize, steps):
    
    # Init omega
    Wxn, Wyn, Wzn = Wx0, Wy0, Wz0

    # Find gamma terms from moments of inertia, and include the -h (step size) constant 
    Y1 = -stepSize * (I3 - I2) / I1
    Y2 = -stepSize * (I1 - I3) / I2
    Y3 = -stepSize * (I2 - I1) / I3

    for n in range(0, steps):

        Kx1 = Y1 * Wyn * Wzn
        Ky1 = Y2 * Wxn * Wzn
        Kz1 = Y3 * Wxn * Wyn
        
        Kx2 = Y1 * (Wyn + (Ky1/2)) * (Wzn + (Kz1/2))
        Ky2 = Y2 * (Wxn + (Kx1/2)) * (Wzn + (Kz1/2))
        Kz2 = Y3 * (Wxn + (Kx1/2)) * (Wyn + (Ky1/2))
    
        Kx3 = Y1 * (Wyn + (Ky2/2)) * (Wzn + (Kz2/2))
        Ky3 = Y2 * (Wxn + (Kx2/2)) * (Wzn + (Kz2/2))
        Kz3 = Y3 * (Wxn + (Kx2/2)) * (Wyn + (Ky2/2))
    
        Kx4 = Y1 * (Wyn + Ky3) * (Wzn + Kz2)
        Ky4 = Y2 * (Wxn + Kx3) * (Wzn + Kz2)
        Kz4 = Y3 * (Wxn + Kx3) * (Wyn + Ky2)
        
        Kx = (Kx1 + (2*Kx2) + (2*Kx3) + Kx4)/6
        Ky = (Ky1 + (2*Ky2) + (2*Ky3) + Ky4)/6
        Kz = (Kz1 + (2*Kz2) + (2*Kz3) + Kz4)/6
        
        Wxn += Kx
        Wyn += Ky
        Wzn += Kz
        
        print("Step %s: %s, %s, %s" % (n, Wxn, Wyn, Wzn))

    #return Wxn, Wyn, Wzn


# Define initial conditions, step size and number of steps
Wx0, Wy0, Wz0 = 1, 2, 3
mass = 5
radius = 2
height = 4

stepSize = 0.1
steps = 200

# Calculate moments of inertia
I1, I2, I3 = GetMomentsOfInertia(mass, height, radius)

# Do the computation
RungeKutta(Wx0, Wy0, Wz0, I1, I2, I3, stepSize, steps)


#input("\nPress Enter to exit.")