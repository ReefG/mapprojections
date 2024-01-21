import matplotlib.image as image
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import math

map_width = 3000
map_length = 3000

def transform(phi, lambdda):
    B = math.cos(phi)*math.sin(lambdda)
    x = 0.5*math.log((1+B)/(1-B))
    y = math.atan(math.tan(phi)/math.cos(lambdda))
    return (y,x)

def inv_transform(x,y):
    phi = math.asin(math.sin(y)/math.cosh(x))
    lambdda = math.atan2(math.sinh(x), math.cos(y))     
    return (phi, lambdda)

def x_to_pixel(x):
    return round((x/6 + 0.5)*map_width)

def pixel_to_x(pixel):
    return (pixel/map_width - 0.5) * 6

def y_to_pixel(y):
    return round((y/6 + 0.5)*map_length)

def pixel_to_y(pixel):
    return (pixel/map_length - 0.5) * 6

def phi_to_pixel(phi):
    return round((phi + 0.5*math.pi)/math.pi*1023)

def lambdda_to_pixel(lambdda):
    return round((lambdda + math.pi)/(2*math.pi)*2047)

def plot_parallels():
    n =12
    parallels = np.linspace(-0.49*math.pi, 0.49*math.pi, num=n) #These have a certain lattitude phi, while longitude changes
    lambddas = np.linspace(-0.5*math.pi, 0.5*math.pi, map_width)

    xy = np.empty((n, map_width, 2))
    xy_2 = np.empty((n, map_width, 2))
    xy_3 = np.empty((n, map_width, 2))
    xy_4 = np.empty((n, map_width, 2))
    
    j = 0
    for parallel in parallels:
        i = 0
        for lambdda in lambddas:
            (y,x) = transform(parallel, lambdda)
            adjusted = (y_to_pixel(y), x_to_pixel(x))
            xy[j][i] = adjusted
            xy_2[j][i] = (y_to_pixel(-y), x_to_pixel(x))
            xy_3[j][i] = (y_to_pixel(-y - math.pi), x_to_pixel(x))
            xy_4[j][i] = (y_to_pixel(y + math.pi), x_to_pixel(x))
            i = i + 1
        
        j = j + 1
        
    plt.plot([0, map_width], [map_length/2, map_length/2],  "w", alpha=0.3)
    for j in range(0, len(parallels)):
        plt.plot(xy[j,:,1], xy[j,:,0], "w", alpha=0.3)
        plt.plot(xy_2[j,:,1], xy_2[j,:,0], "w", alpha=0.3)
        plt.plot(xy_3[j,:,1], xy_3[j,:,0], "w", alpha=0.3)
        plt.plot(xy_4[j,:,1], xy_4[j,:,0], "w", alpha=0.3)
        
        
        
def plot_meridians():
    n = 12
    meridians = np.linspace(-2*math.pi, 2*math.pi, num=n)  # These have a certain longitude lambdda, while latitude changes
    phis = np.linspace(-1/2* math.pi, 1/2*math.pi, map_length)
    xy_meridians = np.empty((n, map_length, 2))
    xy_meridians_top = np.empty((n, map_length, 2))
    xy_meridians_bot = np.empty((n, map_length, 2))
    j = 0
    for meridian in meridians:
        i = 0
        for phi in phis:
            (y, x) = transform(phi, meridian)
            adjusted = (y_to_pixel(y), x_to_pixel(x))
            adjusted_top = ((y_to_pixel(y + math.pi), x_to_pixel(x)))
            adjusted_bot = ((y_to_pixel(y - math.pi), x_to_pixel(x)))
            xy_meridians[j][i] = adjusted
            xy_meridians_top[j][i] = adjusted_top
            xy_meridians_bot[j][i] = adjusted_bot
            i = i + 1

        j = j + 1
    plt.plot([0, map_width], [y_to_pixel(-math.pi/2), y_to_pixel(-math.pi/2)],  "w", alpha=0.3)
    plt.plot([0, map_width], [y_to_pixel(math.pi/2), y_to_pixel(math.pi/2)],  "w", alpha=0.3)
    for j in range(0, len(meridians)):
        plt.plot(xy_meridians[j, :, 1], xy_meridians[j, :, 0], "w", alpha=0.3)
        plt.plot(xy_meridians_top[j, :, 1], xy_meridians_top[j, :, 0], "w", alpha=0.3)
        plt.plot(xy_meridians_bot[j, :, 1], xy_meridians_bot[j, :, 0], "w", alpha=0.3)
        
     
 
world = image.imread("land_shallow_topo_2048.tif") # longitude lambdda by latitude phi, in this map


new_world = np.empty(  ( map_length , map_width , 3 ), dtype=np.uint8 ) # this will be (y,x)

angle_dist = np.empty(  ( map_length , map_width , 1 ))
area_dist = np.empty(  ( map_length , map_width , 1 ))

phi_list = []
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
            B = math.cos(phi)*math.sin(lambdda)
            k = 1/math.sqrt(1-B**2)
            area_dist[y_pixel][x_pixel] = k**2   
            
def area_score(c):
    total = 0
    total_normal = 0
    steps = 1000
    lambddas = np.linspace(-math.pi, math.pi, num=steps)
    phis = np.linspace(-math.pi/2, math.pi/2, num=steps)
    for phi in phis:
        for lambdda in lambddas:
            B = math.cos(phi)*math.sin(lambdda)
            k = c*1/math.sqrt(1-B**2)
            value = abs(math.log(k**2))
            total  = total + value*math.cos(phi)
            total_normal = total_normal + math.cos(phi)
    return total/(total_normal)

def score(c):
    total = 0
    total_normal = 0
    steps = 1000
    lambddas = np.linspace(-math.pi, math.pi, num=steps)
    phis = np.linspace(-math.pi/2, math.pi/2, num=steps)
    for phi in phis:
        for lambdda in lambddas:
            B = math.cos(phi)*math.sin(lambdda)
            k = c*1/math.sqrt(1-B**2)
            value = abs(math.log(k)) + abs(math.log(k))
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
# area_score = area_score()

plt.figure(dpi = 400)
plot_parallels()
plot_meridians()
plt.axis('off')
plt.imshow(new_world)


hm = plt.imshow(area_dist, cmap='RdPu', interpolation='nearest', alpha=0.7, norm=colors.LogNorm())  
# hm = plt.imshow(area_dist, cmap='RdPu', interpolation='nearest', alpha=0.7)  
ax = plt.gca()
cbar = ax.figure.colorbar(hm, ax=ax)