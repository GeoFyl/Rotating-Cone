
# The differential equation
def f(x, y):
    return y - x*x + 1

# Performs Runge-Kutta computation
def RungeKutta(x0, y0, h, steps):
    
    xn = x0
    yn = y0

    for n in range(0, steps):

        k1 = f(xn, yn)

        k2 = f(xn + (h/2), yn + (h*k1/2))
    
        k3 = f(xn + (h/2), yn + (h*k2/2))
    
        k4 = f(xn + h, yn + (h*k3))
    
        k = (k1 + (2*k2) + (2*k3) + k4)/6
        
        yn += (h*k)
        xn += h
        
        print("Step %s: %s" % (n, yn))

    return yn


# Define initial conditions, step size and number of steps
x0 = 0
y0 = 0.5
h = 0.5
steps = 4

# Do the computation
RungeKutta(x0, y0, h, steps)



#input()