# Currently according to lecture notes, so for ellipsoidal object

# Calculate moments of inertia from semi-axis 
def GetMomentsOfInertia(M, a, b, c):
    I1 = M * (b*b + c*c)/5
    I2 = M * (a*a + c*c)/5
    I3 = M * (a*a + b*b)/5
    
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

    return Wxn, Wyn, Wzn


# Define initial conditions, step size and number of steps
Wx0 = 1
Wy0 = 1
Wz0 = 1

Mass = 3.14
a = 3
b = 2
c = 1

stepSize = 0.08
steps = 10

# Calculate moments of inertia
I1, I2, I3 = GetMomentsOfInertia(Mass, a, b, c)

# Do the computation
RungeKutta(Wx0, Wy0, Wz0, I1, I2, I3, stepSize, steps)


#input("\nPress Enter to exit.")