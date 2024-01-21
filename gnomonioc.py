import matplotlib.image as image
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import numpy as np
import math

map_width = 2000
map_length = 2000
h = 0.00001

phi_1 = 0

def transform(phi, lambdda):
    cosc = math.sin(phi_1)*math.sin(phi) + math.cos(phi_1)*math.cos(phi)*math.cos(lambdda)
    k = 1/cosc
    x =k*math.cos(phi)*math.sin(lambdda)
    y= k*(math.cos(phi_1)*math.sin(phi) - math.sin(phi_1)*math.cos(phi)*math.cos(lambdda))
    return (y,x)

def inv_transform(x,y):
    rho = math.sqrt(x**2 + y**2)
    if (rho == 0): return (phi_1, 0)
    c = math.atan(rho)
    phi = math.asin(math.cos(c)*math.sin(phi_1) + (y*math.sin(c)*math.cos(phi_1)/rho))
    if (phi == math.pi/2):
        return (phi, math.atan2(x, -y))
    if (phi == -math.pi/2): 
        return (phi, math.atan2(x,y))
    lambdda = math.atan2(x*math.sin(c), rho*math.cos(phi_1)*math.cos(c) - y*math.sin(phi_1)*math.sin(c))
    return (phi, lambdda)


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

def x_to_pixel(x):
    return round((x/8 + 0.5)*map_width)

def pixel_to_x(pixel):
    return (pixel/map_width - 0.5) * 8

def y_to_pixel(y):
    return round((y/8 + 0.5)*map_length)

def pixel_to_y(pixel):
    return ((pixel)/map_length - 0.5) *8

def phi_to_pixel(phi):
    return round((phi + 0.5*math.pi)/math.pi*1023)

def lambdda_to_pixel(lambdda):
    return round((lambdda + math.pi)/(2*math.pi)*2047)

def plot_parallels(n):
    parallels = np.linspace(-0.49999*math.pi, 0.49999*math.pi, num=n) #These have a certain lattitude phi, while longitude changes
    lambddas = np.linspace(-0.9999*math.pi, 0.9999*math.pi, map_width)
    xy = np.empty((n, map_width, 2))
    j=0
    for parallel in parallels:
        i = 0
        # (y,x) = transform(parallel, -0.9999*math.pi)
        # adjusted = (y_to_pixel(y), x_to_pixel(x))
        adjusted = (0,0)
        for lambdda in lambddas:
            cosc = math.sin(phi_1)*math.sin(parallel) + math.cos(phi_1)*math.cos(parallel)*math.cos(lambdda)
            if (cosc <= 0): 
                xy[j][i] = adjusted
                i = i + 1
                continue
            (y,x) = transform(parallel, lambdda)
            adjusted = (y_to_pixel(y), x_to_pixel(x))
            xy[j][i] = adjusted
            i = i + 1
        j = j+1
        
    for j in range(0, len(parallels)):
        plt.plot(xy[j,:,1], xy[j,:,0], "w", alpha=0.3)
        
def plot_meridians(n):
    meridians = np.linspace(-math.pi, math.pi, num=n)  # These have a certain longitude lambdda, while latitude changes
    phis = np.linspace(-0.5* math.pi, 0.5* math.pi, map_length)
    xy_meridians = np.empty((n, map_length, 2))
    j = 0
    for meridian in meridians:
        i = 0
        (y,x) = transform(-0.5*math.pi, meridian)
        adjusted = (y_to_pixel(y), x_to_pixel(x))
        for phi in phis:
            cosc = math.sin(phi_1)*math.sin(phi) + math.cos(phi_1)*math.cos(phi)*math.cos(meridian)
            if (cosc <= 0): 
                xy_meridians[j][i] = adjusted
                i = i + 1
                continue
            (y, x) = transform(phi, meridian)
            adjusted = (y_to_pixel(y), x_to_pixel(x))
            xy_meridians[j][i] = adjusted
            i = i + 1

        j = j + 1

    for j in range(0, len(meridians)):
        plt.plot(xy_meridians[j, :, 1], xy_meridians[j, :, 0], "w", alpha=0.3)
    
 
world = image.imread("land_shallow_topo_2048.tif") # longitude lambdda by latitude phi, in this map


new_world = np.empty(  ( map_length , map_width , 3 ), dtype=np.uint8 ) # this will be (y,x)

angle_dist = np.empty(  ( map_length , map_width , 1 ))
area_dist = np.empty(  ( map_length , map_width , 1 ))
scale_factors = np.empty((map_length , map_width, 2))
T_matrices = np.empty((map_length , map_width, 2, 2))

phi_list = []
for x_pixel in range(0, map_width):
    for y_pixel in range(0, map_length):
        # print((x,y))
        (x,y) = (pixel_to_x(x_pixel), pixel_to_y(y_pixel))
        # rho = math.sqrt(x**2 + y**2)
        # if (rho > 2): continue
        (phi, lambdda) = inv_transform(x,y)
        cosc = math.sin(phi_1)*math.sin(phi) + math.cos(phi_1)*math.cos(phi)*math.cos(lambdda)
        if (cosc < 0): continue
        (phi_pixel, lambdda_pixel) = (phi_to_pixel(phi), lambdda_to_pixel(lambdda))
        # phi_list.append(phi)
        # print("angles: ", lambdda, phi)
        if (phi_pixel >= 0 and lambdda_pixel >=0 and phi_pixel <= 1023 and lambdda_pixel <= 2047):
            pixel_value = world[phi_pixel][lambdda_pixel]
            new_world[y_pixel][x_pixel] = pixel_value
            # s, T = get_scale_factors_at(phi, lambdda)
            # scale_factors[y_pixel][x_pixel] = s
            # T_matrices[y_pixel][x_pixel] = T
            cosc = math.sin(phi_1)*math.sin(phi) + math.cos(phi_1)*math.cos(phi)*math.cos(lambdda)
            k = 1/cosc
            h = 1/cosc**2
            omega = 2*math.asin(abs(h-k)/(h+k))
            angle_dist[y_pixel][x_pixel] = omega
            area_dist[y_pixel][x_pixel] = h*k

def avg_angle_dist():
    total = 0
    total_normal = 0
    steps = 1000
    lambddas = np.linspace(-math.pi, math.pi, num=steps)
    phis = np.linspace(-math.pi/2, math.pi/2, num=steps*2)
    for phi in phis:
        for lambdda in lambddas:
            # if (1 + math.sin(phi_1)*math.sin(phi) + math.cos(phi_1)*math.cos(lambdda) == 0): continue
            # s, T = get_scale_factors_at(phi, lambdda)
            # omega = 2*math.asin(abs(s[0] - s[1])/(s[0] + s[1]))
            cosc = math.sin(phi_1)*math.sin(phi) + math.cos(phi_1)*math.cos(phi)*math.cos(lambdda)
            if (cosc < 0): continue
            k = 1/cosc
            h = 1/cosc**2
            omega = 2*math.asin(abs(h-k)/(h+k))
            total  = total + omega*math.cos(phi)
            total_normal = total_normal + math.cos(phi)
    return total/(total_normal)


def score(c):
    total = 0
    total_normal = 0
    steps = 1000
    lambddas = np.linspace(-math.pi, math.pi, num=steps*2)
    phis = np.linspace(-math.pi/2, math.pi/2, num=steps)
    for phi in phis:
        for lambdda in lambddas:
            cosc = math.sin(phi_1)*math.sin(phi) + math.cos(phi_1)*math.cos(phi)*math.cos(lambdda)
            if (cosc < 0): continue
            k = c/cosc
            h = c/cosc**2
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
    phis = np.linspace(-math.pi*0.5, math.pi*0.5, num=steps)
    for phi in phis:
        for lambdda in lambddas:
            cosc = math.sin(phi_1)*math.sin(phi) + math.cos(phi_1)*math.cos(phi)*math.cos(lambdda)
            if (cosc < 0): continue
            k = c/cosc
            h = c/cosc**2
            value = abs(math.log(k*h))
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

# area_score = area_score()

# angle_score = avg_angle_dist()

# plt.figure(dpi = 400)
# plt.axis('off')
# plot_parallels(16)
# plot_meridians(16)
# plt.imshow(new_world)

# hm = plt.imshow(angle_dist, cmap="hot", interpolation='nearest', alpha=0.7)  
# ax = plt.gca()
# cbar = ax.figure.colorbar(hm, ax=ax)

plt.figure(dpi = 400)
plt.axis('off')
plot_parallels(16)
plot_meridians(16)
plt.imshow(new_world)

hm = plt.imshow(area_dist, cmap="RdPu", interpolation='nearest', alpha=0.7,  norm=colors.LogNorm())  
ax = plt.gca()
cbar = ax.figure.colorbar(hm, ax=ax)