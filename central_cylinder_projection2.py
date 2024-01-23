import matplotlib.image as image
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import math


plt.figure(dpi=400)

map_width = 2000
map_length = 2000

def x_to_pixel(x):
    return round((x/( 2 *math.pi) + 0.5)*map_width)

def pixel_to_x(pixel):
    return (pixel/map_width - 0.5) * 2 *math.pi

def y_to_pixel(y):
    return round((y/( 2 *math.pi) + 0.5)*map_length)

def pixel_to_y(pixel):
    return (pixel/map_length - 0.5) *  2 *math.pi 

def phi_to_pixel(phi):
    return round((phi + 0.5*math.pi)/math.pi*1023)

def lambdda_to_pixel(lambdda):
    return round((lambdda + math.pi)/(2*math.pi)*2047)

def transform(phi, lambdda): #input is in real numbers aka pi and stuff
    x = lambdda
    y = math.tan(phi)
    return (y,x)
    
    
def inv_transform(x,y): #Input is in reals
    lambdda = x
    phi = math.atan(y)
    return (phi, lambdda)
    

world = image.imread("land_shallow_topo_2048.tif") # longitude theta by latitude phi, in this map

new_world = np.empty(  ( map_length , map_width , 3 ), dtype=np.uint8 ) # this will be (y,x)

area_dist = np.empty(  ( map_length , map_width , 1 ))
angle_dist = np.empty(  ( map_length , map_width , 1 ))


parallels = np.linspace(-0.5*math.pi, 0.5*math.pi, num=9) #These have a certain lattitude phi, while longitude changes
lambddas = np.linspace(-math.pi, math.pi, map_width)
xy = np.empty((9, map_width, 2))
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


meridians = np.linspace(-math.pi, math.pi, num=18)  # These have a certain longitude lambdda, while latitude changes
phis = np.linspace(-0.5 * math.pi, 0.5 * math.pi, map_length)
xy_meridians = np.empty((18, map_length, 2))
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


for x_pixel in range(0, map_width):
    for y_pixel in range(0, map_length):
        (phi, lambdda) = inv_transform(pixel_to_x(x_pixel),pixel_to_y(y_pixel))
        (phi_pixel, lambdda_pixel) = (phi_to_pixel(phi), lambdda_to_pixel(lambdda))
        if (phi_pixel >= 0 and lambdda_pixel >=0 and phi_pixel <= 1023 and lambdda_pixel <= 2047):
            pixel_value = world[phi_pixel][lambdda_pixel]
            new_world[y_pixel][x_pixel] = pixel_value
            sec = 1/math.cos(phi)
            s = sec**3
            omega = 2*math.asin(abs(sec - sec**2)/(sec + sec**2))
            area_dist[y_pixel][x_pixel] = s
            angle_dist[y_pixel][x_pixel] = omega
       
       
def score(c):
    total = 0
    total_normal = 0
    steps = 1000
    lambddas = np.linspace(-math.pi, math.pi, num=steps*2)
    phis = np.linspace(-math.pi/2, math.pi/2, num=steps)
    for phi in phis:
        for lambdda in lambddas:
            sec = 1/math.cos(phi)
            k = c*sec
            h = c*sec**2
            value = abs(math.log(k)) + abs(math.log(h))
            total  = total + value*math.cos(phi)
            total_normal = total_normal + math.cos(phi)
    return total/(total_normal)
 
# score = score()

def area_score(c):
    total = 0
    total_normal = 0
    steps = 1000
    lambddas = np.linspace(-math.pi, math.pi, num=steps*2)
    phis = np.linspace(-math.pi/2, math.pi/2, num=steps)
    for phi in phis:
        for lambdda in lambddas:
            sec = 1/math.cos(phi)
            k = c*sec
            h = c*sec**2
            value = abs(math.log(k*h))
            total  = total + value*math.cos(phi)
            total_normal = total_normal + math.cos(phi)
    return total/(total_normal)


def angle_score():
    total = 0
    total_normal = 0
    steps = 1000
    lambddas = np.linspace(-math.pi, math.pi, num=steps*2)
    phis = np.linspace(-math.pi/2, math.pi/2, num=steps)
    for phi in phis:
        for lambdda in lambddas:
            sec = 1/math.cos(phi)
            k = sec
            h = sec**2
            omega = 2*math.asin(abs(k - h)/(h + k))
            total  = total + omega*math.cos(phi)
            total_normal = total_normal + math.cos(phi)
    return total/(total_normal)

# area_score = area_score()
# angle_score = angle_score()

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

plt.axis('off')
plt.imshow(new_world)

# hm = plt.imshow(area_dist, cmap='RdPu', interpolation='nearest', alpha=0.7,  norm=colors.LogNorm())  
# ax = plt.gca()
# cbar = ax.figure.colorbar(hm, ax=ax)

hm = plt.imshow(angle_dist, cmap='hot', interpolation='nearest', alpha=0.7)  
ax = plt.gca()
cbar = ax.figure.colorbar(hm, ax=ax)














