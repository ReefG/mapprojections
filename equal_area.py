import matplotlib.image as image
import matplotlib.pyplot as plt
import numpy as np
import math

map_width = 6284
map_length = 2000

def transform(phi, lambdda):
    x = lambdda
    y = math.sin(phi)
    return (y,x)

def inv_transform(x,y):
    lambdda = round(2047*x/map_width) 
    phi = math.asin((y - (map_length/2))/(map_length/2)) + 0.5*math.pi #gives angle from 0 to pi
    phi_pixel = round(1023*phi/(math.pi)) #gives pixel value from 0 to 1023
    return (phi_pixel, lambdda, phi - 0.5*math.pi)

def plot_parallels():
    parallels = np.linspace(-0.5*math.pi, 0.5*math.pi, num=12) #These have a certain lattitude phi, while longitude changes
    lambddas = np.linspace(-math.pi, math.pi, map_width)
    xy = np.empty((12, map_width, 2))
    j=0
    for parallel in parallels:
        i = 0
        for lambdda in lambddas:
            (y,x) = transform(parallel, lambdda)
            adjusted = ((y+1)*map_length/2, (x + math.pi)/(2*math.pi)*map_width)
            xy[j][i] = adjusted
            i = i + 1
            
        j = j+1
        
    for j in range(0, len(parallels)):
        plt.plot(xy[j,:,1], xy[j,:,0], "w", alpha=0.3)
        
def plot_meridians():
    meridians = np.linspace(-math.pi, math.pi, num=12)  # These have a certain longitude lambdda, while latitude changes
    phis = np.linspace(-0.5 * math.pi, 0.5 * math.pi, map_length)
    xy_meridians = np.empty((12, map_length, 2))
    j = 0
    for meridian in meridians:
        i = 0
        for phi in phis:
            (y, x) = transform(phi, meridian)
            adjusted = ((y+1)*map_length/2, (x + math.pi) / (2 * math.pi) * map_width)
            xy_meridians[j][i] = adjusted
            i = i + 1

        j = j + 1

    for j in range(0, len(meridians)):
        plt.plot(xy_meridians[j, :, 1], xy_meridians[j, :, 0], "w", alpha=0.3)

 
world = image.imread("land_shallow_topo_2048.tif") # longitude lambdda by latitude phi, in this map


new_world = np.empty(  ( map_length , map_width , 3 ), dtype=np.uint8 ) # this will be (y,x)

angle_dist = np.empty(  ( map_length , map_width , 1 ))

phi_list = []
for x in range(0, map_width):
    for y in range(0, map_length):
        # print((x,y))
        (phi, lambdda, phi_value) = inv_transform(x,y)
        # phi_list.append(phi)
        # print("angles: ", lambdda, phi)
        pixel_value = world[phi][lambdda]
        new_world[y][x] = pixel_value
        k = 1/math.cos(phi_value)
        omega = 2*math.asin(abs(k - 1/k)/(k + 1/k))
        angle_dist[y][x] = omega
        
        

        
# phi_list.sort()
# print(phi_list)
        
# for x in range(0, len(world[1])):
#     for y in range(0, len(world)):
#         if len(new_world[x][y]) < 3:
#             print((x,y), new_world[x][y])
        
# plt.imshow(world)

def avg_angle_dist():
    total = 0
    total_normal = 0
    steps = 1000
    lambddas = np.linspace(-math.pi, math.pi, num=steps)
    phis = np.linspace(-math.pi/2, math.pi/2, num=steps*2)
    for phi in phis:
        for lambdda in lambddas:
            k = 1/math.cos(phi)
            omega = 2*math.asin(abs(k - 1/k)/(k + 1/k))
            total  = total + omega*math.cos(phi)
            total_normal = total_normal + math.cos(phi)
    return total/(total_normal)

def score():
    total = 0
    total_normal = 0
    steps = 1000
    lambddas = np.linspace(-math.pi, math.pi, num=steps*2)
    phis = np.linspace(-math.pi/2, math.pi/2, num=steps)
    for phi in phis:
        for lambdda in lambddas:
            k = 1/math.cos(phi)
            value = abs(math.log(k)) + abs(math.log(1/k))
            total  = total + value*math.cos(phi)
            total_normal = total_normal + math.cos(phi)
    return total/(total_normal)
 
score = score()
avg = avg_angle_dist()


plot_parallels()
plot_meridians()

plt.imshow(new_world)

# hm = plt.imshow(angle_dist, cmap='hot', interpolation='nearest', alpha=0.7)  
# ax = plt.gca()
# cbar = ax.figure.colorbar(hm, ax=ax)