import matplotlib.image as image
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import math

map_width = 2000
map_length = 2000

h = 0.00001

phi_1 = math.pi/4
phi_2 = (math.pi)/(3)
n = math.log(math.cos(phi_1)/math.cos(phi_2))/math.log(math.tan(math.pi/4 + phi_2/2)/math.tan(math.pi/4 + phi_1/2))
F = math.cos(phi_1)*(math.tan(math.pi/4 + phi_1/2)**n)/n
phi_0 = 0
rho_0 = F/math.tan(math.pi/4 + phi_0/2)**n
sign = 1
if (n >= 0): sign = 1
else: sign = -1


def transform(phi, lambdda):
    C = math.tan(math.pi/4 + phi/2)**n
    rho = F/C
    x = rho*math.sin(n*lambdda)
    y = rho_0-rho*math.cos(n*lambdda)
    return (y,x)

def inv_transform(x,y):
    rho = sign*math.sqrt(x**2 + (rho_0 - y)**2)
    if ((x,y) == (0,0)): return (0,0)
    theta = math.atan2(x, (rho_0-y))
    phi = 2*math.atan((F/rho)**(1/n)) - 0.5*math.pi
    lambdda = theta/n
    return (phi, lambdda)

def x_to_pixel(x):
    return round((x/16 + 0.5)*map_width)

def pixel_to_x(pixel):
    return (pixel/map_width - 0.5) * 16

def y_to_pixel(y):
    return round((y/16 + 0.5)*map_length)

def pixel_to_y(pixel):
    return ((pixel)/map_length - 0.5) * 16

def phi_to_pixel(phi):
    return round((phi + 0.5*math.pi)/math.pi*1023)

def lambdda_to_pixel(lambdda):
    return round((lambdda + math.pi)/(2*math.pi)*2047)

def get_scale_factors_at(phi, lambdda):
    (y_phi_h,x_phi_h) = transform(phi + h, lambdda)
    (y_phi_minush,x_phi_minush) = transform(phi - h, lambdda)
    
    (y_l_h,x_l_h) = transform(phi, lambdda + h)
    (y_l_minush,x_l_minush) = transform(phi, lambdda - h)
    
    x_phi = ( x_phi_h - x_phi_minush ) / (2*h)
    x_lambdda = (x_l_h - x_l_minush ) / (2*h)
    y_phi = ( y_phi_h - y_phi_minush ) / (2*h)
    y_lambdda = (y_l_h - y_l_minush) / (2*h)
    
    T = [[x_lambdda/(math.cos(phi)), x_phi], [y_lambdda/math.cos(phi), y_phi]]
    
    s = np.linalg.svd(T, compute_uv=False)
    return s, T

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
        
        standard_parallels = np.empty((2, map_width, 2))
        i_lambda = 0
        for lambdda in lambddas:
            (y1,x1) = transform(phi_1, lambdda)
            adjusted_1 = (y_to_pixel(y1), x_to_pixel(x1))
            (y2,x2) = transform(phi_2, lambdda)
            adjusted_2 = (y_to_pixel(y2), x_to_pixel(x2))
            standard_parallels[0][i_lambda] = adjusted_1
            standard_parallels[1][i_lambda] = adjusted_2
            i_lambda = i_lambda + 1
            
        # plt.plot(standard_parallels[0,:,1], standard_parallels[0,:,0], "r", alpha=0.3)
        # plt.plot(standard_parallels[1,:,1], standard_parallels[1,:,0], "r", alpha=0.3)
        
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

 
world = image.imread("land_shallow_topo_2048.tif") # longitude lambdda by latitude phi, in this map


new_world = np.empty(  ( map_length , map_width , 3 ), dtype=np.uint8 ) # this will be (y,x)

area_dist = np.empty(  ( map_length , map_width , 1 ))

# phi_list = []
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
            C = math.tan(phi/2 + math.pi/4)**(2*n+1)
            if (C != 0): 
                factor = 0.5*F**2*n**2*(1/math.cos(phi/2 + math.pi/4)**2)/C            
                area_dist[y_pixel][x_pixel] = factor

# def score():
#     total = 0
#     total_normal = 0
#     steps = 1000
#     lambddas = np.linspace(-math.pi, math.pi, num=steps*2)
#     phis = np.linspace(-math.pi/2, math.pi/2, num=steps)
#     for phi in phis:
#         for lambdda in lambddas:
#             C1 = (math.cos(phi)*math.tan(math.pi/4 + phi/2)**n)
#             if (C1 == 0): 
#                 continue
#             k = math.cos(phi_1)*math.tan(math.pi/4 + phi_1/2)**n/C1
#             h = k
#             C2 = math.tan(phi/2 + math.pi/4)**(2*n+1)
#             # if (C2 == 0): 
#             #     continue
#             s = 0.5*F**2*n**2*(1/math.cos(phi/2 + math.pi/4)**2)/C2     
#             a1 = math.sqrt(k**2 + h**2 + 2*s)
#             b1 = math.sqrt(k**2 + h**2 - 2*s)
#             a = (a1 + b1)/2
#             b = (a1-b1)/2
#             if (a <= 0 or b <= 0):
#                 continue
#             value = abs(math.log(a)) + abs(math.log(b))
#             total  = total + value*math.cos(phi)
#             total_normal = total_normal + math.cos(phi)
#     return total/(total_normal)

def score():
    total = 0
    total_normal = 0
    steps = 1000
    lambddas = np.linspace(-math.pi, math.pi, num=steps*2)
    phis = np.linspace(-math.pi*0.49999, math.pi*0.49999, num=steps)
    for phi in phis:
        for lambdda in lambddas:
            s, T = get_scale_factors_at(phi, lambdda)
            if (s[0] <= 0 or s[1] <= 0):
                continue
            value = abs(math.log(s[0])) + abs(math.log(s[1]))
            total  = total + value*math.cos(phi)
            total_normal = total_normal + math.cos(phi)
    return total/(total_normal)

def area_score2():
    total = 0
    total_normal = 0
    steps = 1000
    lambddas = np.linspace(-math.pi, math.pi, num=steps*2)
    phis = np.linspace(-math.pi*0.49999, math.pi*0.49999, num=steps)
    for phi in phis:
        for lambdda in lambddas:
            s = ((F*n*(1/math.cos(phi)))**2)/math.tan(math.pi/4 + phi/2)**(2*n)
            value = abs(math.log(s))
            total  = total + value*math.cos(phi)
            total_normal = total_normal + math.cos(phi)
    return total/(total_normal)
 
# score = score()

def area_score():
    total = 0
    total_normal = 0
    steps = 1000
    lambddas = np.linspace(-math.pi, math.pi, num=steps*2)
    phis = np.linspace(-math.pi*0.49999, math.pi*0.49999, num=steps)
    for phi in phis:
        for lambdda in lambddas:
            s, T = get_scale_factors_at(phi, lambdda)
            # if (s[0] <= 0 or s[1] <= 0):
            #     continue
            value = abs(math.log(s[0]*s[1]))
            total  = total + value*math.cos(phi)
            total_normal = total_normal + math.cos(phi)
    return total/(total_normal)

def s_min():
    steps =1000
    lambddas = np.linspace(-math.pi, math.pi, num=steps*2)
    phis = np.linspace(-math.pi*0.4999, math.pi/2, num=steps)
    smin = 100000000000000000
    for phi in phis:
        for lambdda in lambddas:
            s = ((F*n*(1/math.cos(phi)))**2)/math.tan(math.pi/4 + phi/2)**(2*n)
            if (s < smin): smin = s
    return smin

def area_score3():
    total = 0
    total_normal = 0
    steps = 1000
    lambddas = np.linspace(-math.pi, math.pi, num=steps*2)
    phis = np.linspace(-math.pi*0.4999, math.pi/2, num=steps)
    smin = s_min()           
    for phi in phis:
        for lambdda in lambddas:
            s =((F*n*(1/math.cos(phi)))**2)/math.tan(math.pi/4 + phi/2)**(2*n) / smin
            if (s <= 0): continue
            value = abs(math.log(s))
            total  = total + value*math.cos(phi)
            total_normal = total_normal + math.cos(phi)
    return total/(total_normal)

# area_score = area_score()

plt.figure(dpi = 400)
plot_parallels()
plot_meridians()
plt.axis('off')
plt.imshow(new_world)

# hm = plt.imshow(area_dist, cmap='RdPu', interpolation='nearest', alpha=0.7, norm=colors.LogNorm(vmin = 10**-1, vmax=10**2))  
# ax = plt.gca()
# cbar = ax.figure.colorbar(hm, ax=ax)