import math
import matplotlib.pyplot as plt
import matplotlib as mpl
#import numpy as np

# Calculate moments of inertia for cone from mass, height and base radius 
def GetMomentsOfInertia(M, h, r):
    I1 = I2 = 3 * M * (r*r + (h*h)/4) / 20
    I3 = 6 * M * r*r / 20
    
    return I1, I2, I3

# Performs Runge-Kutta computation
def RungeKutta():
    
    # Init omega
    Wxn, Wyn, Wzn = Wx0, Wy0, Wz0

    # Find gamma terms from moments of inertia, and include the -h (step size) constant 
    Y1 = -stepSize * (I3 - I2) / I1
    Y2 = -stepSize * (I1 - I3) / I2
    Y3 = -stepSize * (I2 - I1) / I3
    
    # Init arrays for storing W components for plotting
    WxArray, WyArray, WzArray = [Wxn], [Wyn], [Wzn]
    WmagArray = [math.sqrt(Wxn*Wxn + Wyn*Wyn + Wzn*Wzn)]

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
        
        Wxn += (Kx1 + (2*Kx2) + (2*Kx3) + Kx4)/6
        Wyn += (Ky1 + (2*Ky2) + (2*Ky3) + Ky4)/6
        Wzn += (Kz1 + (2*Kz2) + (2*Kz3) + Kz4)/6
        
        # Calculate magnitude of W (the angular speed)
        Wmag = math.sqrt(Wxn*Wxn + Wyn*Wyn + Wzn*Wzn)

        WxArray.append(Wxn)
        WyArray.append(Wyn)
        WzArray.append(Wzn)
        WmagArray.append(Wmag)

        print("Step %s: %s, %s, %s, magnitude: %s" % (n, Wxn, Wyn, Wzn, Wmag))
        
    return WxArray, WyArray, WzArray, WmagArray

def PlotTask1():
    # Populate array with values for plot x-axis (time)
    timeArray = []
    for i in range(0, steps + 1):
        timeArray.append(i * stepSize)

    # Plot graph
    plt.plot(timeArray, WxArray, label="ωx")
    plt.plot(timeArray, WyArray, label="ωy")
    plt.plot(timeArray, WzArray, label="ωz")
    plt.plot(timeArray, WmagArray, label="|ω|")
    plt.legend(loc="lower right")
    plt.title("Task 1 plots")
    plt.xlabel("time (s)")
    plt.ylabel("ωi (rad/s)")
    plt.show()    



# -----========== Main  ==========-----

# Define initial conditions, step size and number of steps
Wx0, Wy0, Wz0 = 3, 1, 2
mass = 5
radius = 2
height = 6

stepSize = 0.1
steps = 200

# Calculate moments of inertia
I1, I2, I3 = GetMomentsOfInertia(mass, height, radius)

# Do the computation
WxArray, WyArray, WzArray, WmagArray = RungeKutta()

# Produce plots of solutions for task 1
PlotTask1()

#input("\nPress Enter to exit.")