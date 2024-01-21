import matplotlib.image as image
import matplotlib.pyplot as plt
import numpy as np
import math

map_width = 3000
map_length = 2000

phi_list =[]
arg_list = []

phi_1 = -math.pi/4
phi_2 = math.pi/3

n = 0.5*(math.sin(phi_1) + math.sin(phi_2))
C = math.cos(phi_1)**2 + 2*n*math.sin(phi_1)
p_0 = math.sqrt(C - 2*n*math.sin(0))/n

upper = (C + 2*n)/(n**2)
lower = (C - 2*n)/(n**2)

def transform(phi, lambdda):
    rho = math.sqrt(C -2*n*math.sin(phi))/n
    x = rho*math.sin(n*lambdda)
    y = p_0 - rho*math.cos(n*lambdda)
    return (y,x)

def inv_transform(x_pixel,y_pixel):
    x = (x_pixel/map_width - 0.5) * 9
    y = (y_pixel/map_length - 0.5) * 6 
    
    circle = x**2 + (p_0 - y)**2
    
    if (circle > upper or circle < lower):
        return (0,0, 0.5*math.pi)
    
    p = math.sqrt(x**2 + (p_0 - y)**2)
    theta = math.atan2(x,(p_0 - y))
    
    arg = (C - p**2*n**2)/(2*n)
    arg_list.append(arg)
    phi = math.asin(arg) + math.pi/2 #value from 0 to pi
    theta = theta/n + math.pi/(2*n) #from 0 to pi/n
    phi_list.append(phi)
    phi_pixel = round(phi/math.pi*1023)
    lambdda_pixel = round(theta*n/math.pi * 2047)
    
    return (phi_pixel, lambdda_pixel, phi - 0.5*math.pi)

def plot_parallels():
    parallels = np.linspace(-0.5*math.pi, 0.5*math.pi, num=10) #These have a certain lattitude phi, while longitude changes
    lambddas = np.linspace(-math.pi, math.pi, map_width)
    xy = np.empty((10, map_width, 2))
    j=0
    for parallel in parallels:
        i = 0
        for lambdda in lambddas:
            (y,x) = transform(parallel, lambdda)
            adjusted = ((y/6 + 0.5)*map_length, (x/9 + 0.5)*map_width,)
            xy[j][i] = adjusted
            i = i + 1
            
        j = j+1
        
    for j in range(0, len(parallels)):
        plt.plot(xy[j,:,1], xy[j,:,0], "w", alpha=0.3)
        
def plot_meridians():
    meridians = np.linspace(-math.pi, math.pi, num=14)  # These have a certain longitude lambdda, while latitude changes
    phis = np.linspace(-0.5 * math.pi, 0.5 * math.pi, map_length)
    xy_meridians = np.empty((14, map_length, 2))
    j = 0
    for meridian in meridians:
        i = 0
        for phi in phis:
            (y, x) = transform(phi, meridian)
            adjusted = ((y/6 + 0.5)*map_length, (x/9 + 0.5)*map_width,)
            xy_meridians[j][i] = adjusted
            i = i + 1

        j = j + 1

    for j in range(0, len(meridians)):
        plt.plot(xy_meridians[j, :, 1], xy_meridians[j, :, 0], "w", alpha=0.3)
    

 
world = image.imread("land_shallow_topo_2048.tif") # longitude theta by latitude phi, in this map


new_world = np.empty(  ( map_length , map_width , 3 ), dtype=np.uint8 ) # this will be (y,x)
angle_dist = np.empty(  ( map_length , map_width , 1 ))

for x_pixel in range(0, map_width):
    for y_pixel in range(0, map_length):
        x = (x_pixel/map_width - 0.5) * 9
        y = (y_pixel/map_length - 0.5) * 6 
        
        circle = x**2 + (p_0 - y)**2
        if (circle > upper or circle < lower):
            continue
        
        (phi, lambdda, phi_value) = inv_transform(x_pixel,y_pixel)
        if (phi >= 0 and  lambdda >=0 and phi <= 1023 and lambdda <= 2047):
            pixel_value = world[phi][lambdda]
            new_world[y_pixel][x_pixel] = pixel_value
            k = math.sqrt(C-2*n*math.sin(phi_value))/math.cos(phi_value)
            factor = 2*math.asin(abs(k - 1/k)/(k + 1/k))
            angle_dist[y_pixel][x_pixel] = factor
        

def avg_angle_dist():
    total = 0
    total_normal = 0
    steps = 1000
    lambddas = np.linspace(-math.pi, math.pi, num=steps)
    phis = np.linspace(-math.pi/2, math.pi/2, num=steps*2)
    for phi in phis:
        for lambdda in lambddas:
            k = math.sqrt(C-2*n*math.sin(phi))/math.cos(phi)
            omega = 2*math.asin(abs(k - 1/k)/(k + 1/k))
            total  = total + omega*math.cos(phi)
            total_normal = total_normal + math.cos(phi)
    return total/(total_normal)
 
avg_angle_dist()

def score():
    total = 0
    total_normal = 0
    steps = 1000
    lambddas = np.linspace(-math.pi, math.pi, num=steps*2)
    phis = np.linspace(-math.pi/2, math.pi/2, num=steps)
    for phi in phis:
        for lambdda in lambddas:
            k = math.sqrt(C-2*n*math.sin(phi))/math.cos(phi)
            # a1 = math.sqrt(k**2 + (1/k)**2 + 2)
            # b1 = math.sqrt(k**2 + (1/k)**2 - 2)
            # a = (a1 + b1)/2
            # b = (a1-b1)/2
            # if (a <= 0 or b <= 0):
            #     continue
            value = abs(math.log(k)) + abs(math.log(1/k))
            total  = total + value*math.cos(phi)
            total_normal = total_normal + math.cos(phi)
    return total/(total_normal)
 
score = score()
   
plot_parallels()      
plot_meridians()  
        
plt.imshow(new_world)
hm = plt.imshow(angle_dist, cmap='hot', interpolation='nearest', alpha=0.5)  

ax = plt.gca()
cbar = ax.figure.colorbar(hm, ax=ax)
