import matplotlib.image as image
import matplotlib.pyplot as plt
import numpy as np
import math

map_width = 3000
map_length = 1500

def transform(phi, lambdda):
    x = lambdda*math.cos(phi)
    y = phi
    return (y,x)

def inv_transform(x,y):
    phi = y
    lambdda = x/math.cos(phi)
    return (phi, lambdda)

def x_to_pixel(x):
    return round((x/(math.pi*2) + 0.5)*map_width)

def pixel_to_x(pixel):
    return (pixel/map_width - 0.5) *  math.pi * 2

def y_to_pixel(y):
    return round((y/math.pi  + 0.5)*map_length)

def pixel_to_y(pixel):
    return (pixel/map_length - 0.5) * math.pi 

def phi_to_pixel(phi):
    return round((phi + 0.5*math.pi)/math.pi*1023)

def lambdda_to_pixel(lambdda):
    return round((lambdda + math.pi)/(2*math.pi)*2047)

def plot_parallels():
    parallels = np.linspace(-0.4999*math.pi, 0.4999*math.pi, num=12) #These have a certain lattitude phi, while longitude changes
    lambddas = np.linspace(-math.pi, math.pi, map_width)
    xy = np.empty((12, map_width, 2))
    j=0
    for parallel in parallels:
        i = 0
        for lambdda in lambddas:
            (y,x) = transform(parallel, lambdda)
            adjusted = (y_to_pixel(y), x_to_pixel(x))
            xy[j][i] = adjusted
            i = i + 1
            
        j = j+1
        
    for j in range(0, len(parallels)):
        plt.plot(xy[j,:,1], xy[j,:,0], "w", alpha=0.3)
        
def plot_meridians():
    meridians = np.linspace(-math.pi, math.pi, num=12)  # These have a certain longitude lambdda, while latitude changes
    phis = np.linspace(-0.4999* math.pi, 0.4999* math.pi, map_length)
    xy_meridians = np.empty((12, map_length, 2))
    j = 0
    for meridian in meridians:
        i = 0
        for phi in phis:
            (y, x) = transform(phi, meridian)
            adjusted = (y_to_pixel(y), x_to_pixel(x))
            xy_meridians[j][i] = adjusted
            i = i + 1

        j = j + 1

    for j in range(0, len(meridians)):
        plt.plot(xy_meridians[j, :, 1], xy_meridians[j, :, 0], "w", alpha=0.3)

def avg_angle_dist():
    total = 0
    total_normal = 0
    steps = 5000
    lambddas = np.linspace(-math.pi, math.pi, num=steps)
    phis = np.linspace(-math.pi/2, math.pi/2, num=steps*2)
    for phi in phis:
        for lambdda in lambddas:
            # k = 1
            # h = math.sqrt(1 + lambdda**2*math.sin(phi)**2)
            angle = 2*math.atan(abs(0.5*lambdda*math.sin(phi)))
            total  = total + angle*math.cos(phi)
            total_normal = total_normal + math.cos(phi)
    return total/(total_normal)

      
 
world = image.imread("land_shallow_topo_2048.tif") # longitude lambdda by latitude phi, in this map


new_world = np.empty(  ( map_length , map_width , 3 ), dtype=np.uint8 ) # this will be (y,x)

angle_dist = np.empty(  ( map_length , map_width , 1 ))
area_dist = np.empty(  ( map_length , map_width , 1 ))

for x_pixel in range(0, map_width):
    for y_pixel in range(0, map_length):
        # print((x,y))
        (x,y) = (pixel_to_x(x_pixel), pixel_to_y(y_pixel))
        (phi, lambdda) = inv_transform(x,y)
        (phi_pixel, lambdda_pixel) = (phi_to_pixel(phi), lambdda_to_pixel(lambdda))
        # phi_list.append(phi)
        # print("angles: ", lambdda, phi)
        if (phi_pixel >= 0 and lambdda_pixel >=0 and phi_pixel <= 1023 and lambdda_pixel <= 2047):
            pixel_value = world[phi_pixel][lambdda_pixel]
            new_world[y_pixel][x_pixel] = pixel_value
            k = 1
            h = math.sqrt(1 + lambdda**2*math.sin(phi)**2)
            angle = 2*math.atan(abs(0.5*lambdda*math.sin(phi)))
            angle_dist[y_pixel][x_pixel] = angle        

def score(c):
    total = 0
    total_normal = 0
    steps = 1000
    lambddas = np.linspace(-math.pi, math.pi, num=steps*2)
    phis = np.linspace(-math.pi/2, math.pi/2, num=steps)
    for phi in phis:
        for lambdda in lambddas:
            k = 1
            h = (math.sqrt(1 + lambdda**2*math.sin(phi)**2))
            a1 = math.sqrt(k**2 + h**2 + 2)
            b1 = math.sqrt(k**2 + h**2 - 2)
            a = (a1 + b1)/2
            b = (a1-b1)/2
            if (a <= 0 or b <= 0):
                continue
            value = abs(math.log(c*a)) + abs(math.log(c*b))
            total  = total + value*math.cos(phi)
            total_normal = total_normal + math.cos(phi)
    return total/(total_normal)

gr = (math.sqrt(5) + 1) / 2
def gss(f, a, b, tol=1e-3):
    c = b - (b - a) / gr
    d = a + (b - a) / gr
    while abs(b - a) > tol:
        fc = f(c)
        fd = f(d)
        if fc < fd:  
            b = d
        else:
            a = c
        # We recompute both c and d here to avoid loss of precision which may lead to incorrect results or infinite loop
        c = b - (b - a) / gr
        d = a + (b - a) / gr
        
        
    return (((b + a) / 2), f((b + a) /2))
 
# score = score()
        
# avg = avg_angle_dist()

plt.figure(dpi = 400)
plt.axis('off')
plot_parallels()
plot_meridians()

plt.imshow(new_world)

hm = plt.imshow(angle_dist, cmap='hot', interpolation='nearest', alpha=0.7)  
ax = plt.gca()
cbar = ax.figure.colorbar(hm, ax=ax)