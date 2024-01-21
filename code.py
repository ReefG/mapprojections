import matplotlib.image as image
import matplotlib.pyplot as plt
import numpy as np
import math

map_width = 6284
map_length = 2000

def inv_transform(x,y):
    theta = round(2047*x/map_width) 
    phi = math.asin((y - (map_length/2))/(map_length/2)) + 0.5*math.pi #gives angle from 0 to pi
    phi_pixel = round(1023*phi/(math.pi)) #gives pixel value from 0 to 1023
    return (phi_pixel, theta)

 
world = image.imread("land_shallow_topo_2048.tif") # longitude theta by latitude phi, in this map


new_world = np.empty(  ( map_length , map_width , 3 ), dtype=np.uint8 ) # this will be (y,x)
# new_world = world.copy()

phi_list = []
for x in range(0, map_width):
    for y in range(0, map_length):
        # print((x,y))
        (phi, theta) = inv_transform(x,y)
        phi_list.append(phi)
        # print("angles: ", theta, phi)
        pixel_value = world[phi][theta]
        new_world[y][x] = pixel_value
        
# phi_list.sort()
# print(phi_list)
        
# for x in range(0, len(world[1])):
#     for y in range(0, len(world)):
#         if len(new_world[x][y]) < 3:
#             print((x,y), new_world[x][y])
        
        

plt.imshow(new_world)