import math
import matplotlib.pyplot as plt

# -----========== Task 1 and 2  ==========-----

# Calculate moments of inertia for cone from mass, height and base radius 
def GetMomentsOfInertia(M, h, r):
    I1 = I2 = 3 * M * (r*r + (h*h)/4) / 20
    I3 = 6 * M * r*r / 20
    
    return I1, I2, I3

def RungeKutta():
    print("\n-------- Runge-Kutta --------\n")

    # Find gamma terms from moments of inertia, and include the -h (step size) constant 
    Y1 = -stepSize * (I3 - I2) / I1
    Y2 = -stepSize * (I1 - I3) / I2
    Y3 = -stepSize * (I2 - I1) / I3
    
    # Init arrays for storing W components for plotting
    Wx, Wy, Wz = [Wx0], [Wy0], [Wz0]
    Wmagnitude = [math.sqrt(Wx0*Wx0 + Wy0*Wy0 + Wz0*Wz0)]

    # Carry out computation
    for n in range(0, steps):

        Kx1 = Y1 * Wy[n] * Wz[n]
        Ky1 = Y2 * Wx[n] * Wz[n]
        Kz1 = Y3 * Wx[n] * Wy[n]
        
        Kx2 = Y1 * (Wy[n] + (Ky1/2)) * (Wz[n] + (Kz1/2))
        Ky2 = Y2 * (Wx[n] + (Kx1/2)) * (Wz[n] + (Kz1/2))
        Kz2 = Y3 * (Wx[n] + (Kx1/2)) * (Wy[n] + (Ky1/2))
    
        Kx3 = Y1 * (Wy[n] + (Ky2/2)) * (Wz[n] + (Kz2/2))
        Ky3 = Y2 * (Wx[n] + (Kx2/2)) * (Wz[n] + (Kz2/2))
        Kz3 = Y3 * (Wx[n] + (Kx2/2)) * (Wy[n] + (Ky2/2))
    
        Kx4 = Y1 * (Wy[n] + Ky3) * (Wz[n] + Kz2)
        Ky4 = Y2 * (Wx[n] + Kx3) * (Wz[n] + Kz2)
        Kz4 = Y3 * (Wx[n] + Kx3) * (Wy[n] + Ky2)
        
        # Calculate and append final values
        Wx.append(Wx[n] + (Kx1 + (2*Kx2) + (2*Kx3) + Kx4)/6)
        Wy.append(Wy[n] + (Ky1 + (2*Ky2) + (2*Ky3) + Ky4)/6)
        Wz.append(Wz[n] + (Kz1 + (2*Kz2) + (2*Kz3) + Kz4)/6)
        
        # Calculate magnitude of W (the angular speed)
        Wmag = math.sqrt(Wx[n+1]*Wx[n+1] + Wy[n+1]*Wy[n+1] + Wz[n+1]*Wz[n+1])
        Wmagnitude.append(Wmag)

        print("Step %s: %s, %s, %s, magnitude: %s" % (n, Wx[n+1], Wy[n+1], Wz[n+1], Wmag))
        
    return Wx, Wy, Wz, Wmagnitude

def PlotTask1():
    # Populate array with values for plot's x-axis (time)
    timeArray = []
    for i in range(0, steps + 1):
        timeArray.append(i * stepSize)

    # Plot graph
    plt.plot(timeArray, WxArray, label="ωx")
    plt.plot(timeArray, WyArray, label="ωy")
    plt.plot(timeArray, WzArray, label="ωz")
    plt.plot(timeArray, WmagArray, label="|ω|")
    #plt.legend(loc="lower left", bbox_to_anchor=(1, 0))
    plt.legend(loc="lower right")
    plt.title("Task 1 Angular Velocity")
    plt.xlabel("time (s)")
    plt.ylabel("ωi (rad/s)")   


# -----========== Task 3 and 4  ==========-----

def SemiImplicitEuler():
    print("\n-------- Semi-implicit Euler --------\n")    

    # Init acceleration, displacement and velocity
    a = -9.8 
    x = [0]
    v = [v0]
    
    # Carry out computation
    for n in range(0, steps):
        v.append(v[n] + stepSize * a)
        x.append(x[n] + stepSize * v[n+1])
        print("Step %s: v: %s, x: %s" % (n, v[n+1], x[n+1]))
    
    return x, v
   
def PlotTask3():
    # Populate array with values for plot's x-axis (time)
    timeArray = []
    for i in range(0, steps + 1):
        timeArray.append(i * stepSize)

    # Plot graph
    plt.figure()
    plt.plot(timeArray, v)
    plt.title("Task 3 Vertical Velocity")
    plt.xlabel("time (s)")
    plt.ylabel("velocity (m/s)")
   
    plt.figure()
    plt.plot(timeArray, x)
    plt.title("Task 3 Vertical Displacement")
    plt.xlabel("time (s)")
    plt.ylabel("displacement (m)")


# -----========== Task 5  ==========-----


# -----========== Main  ==========-----

# ---- Runge-Kutta ----

# Define initial conditions, step size and number of steps
Wx0, Wy0, Wz0 = 1, 2, 3
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

# ---- Semi-implicit Euler ----

# Define initial vertical velocity
v0 = 200

# Compute trajectory
x, v = SemiImplicitEuler()

# Produce plots of solutions for task 3
PlotTask3()


# Show the graphs
plt.show() 


#input("\nPress Enter to exit.")