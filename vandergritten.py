import matplotlib.image as image
import matplotlib.pyplot as plt
import matplotlib.colors as colors

import numpy as np
import math

map_width = 2000
map_length = 2000
h = 0.001

def transform(phi, lambdda):
    if (phi == 0): 
        return (0, lambdda)
    
    sign_phi = 1
    if (phi >= 0): sign_phi = 1
    else: sign_phi = -1    
    
    if (phi > math.pi/2):
        return transform(math.pi/2, lambdda)
    
    if (phi < -math.pi/2):
        return transform(-math.pi/2, lambdda)
    
    theta = math.asin(abs(2*phi/math.pi))
    
    if (phi == math.pi/2 or phi == -math.pi/2 or lambdda == 0): 
        return (math.pi*sign_phi*math.tan(theta/2), 0)
    
    A = 0.5*abs(math.pi/lambdda - lambdda/math.pi)
    
    value = (math.sin(theta) + math.cos(theta) - 1)
    if (value == 0):
        return (0, lambdda)
    
    G = math.cos(theta)/value
    P = G*(2/math.sin(theta) - 1)
    Q = A**2 + G
    
    sign_lambdda = 1
    if (lambdda >= 0): sign_lambdda = 1
    else: sign_lambdda = -1
    
    x = sign_lambdda*math.pi*(A*(G - P**2) + math.sqrt(A**2*(G - P**2)**2 - (P**2 + A**2)*(G**2 - P**2)))/(P**2 + A**2)
    y = sign_phi*math.pi*abs(P*Q-A*math.sqrt((A**2 + 1)*(P**2 + A**2) - Q**2))/(P**2 + A**2)
    return (y,x)

def inv_transform(x,y):
    X = x/math.pi
    Y = y/math.pi
    c1 = -abs(Y)*(1 + X**2 + Y**2)
    c2 = c1 - 2*Y**2 + X**2
    c3 = -2*c1 + 1 + 2*Y**2 + (X**2 + Y**2)**2
    d = Y**2/c3 + (2*c2**3/(c3**3) - 9*c1*c2/(c3**2))/27
    a1 = (c1 - c2**2/(3*c3))/c3
    m1 = 2*math.sqrt(-a1/3)
    
    if (a1*m1 == 0): return (0,0)
    
    value = 3*d/(a1*m1)
    
    if (value > 1 or value < -1):
        return (0,0)
    
    theta_1 = 1/3*math.acos(value)

    
    sign_y = 1
    if (y >= 0): sign_y = 1
    else: sign_y = -1
    
    phi = sign_y*math.pi*(-m1*math.cos(theta_1 + math.pi/3) -c2/(3*c3))
    if (X == 0): return(phi, 0)
    
    lambdda = math.pi*(X**2 + Y**2 - 1 + math.sqrt(1 + 2*(X**2 - Y**2) + (X**2 + Y**2)**2))/(2*X)
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
    return round((x/7 + 0.5)*map_width)

def pixel_to_x(pixel):
    return (pixel/map_width - 0.5) * 7

def y_to_pixel(y):
    return round((y/7 + 0.5)*map_length)

def pixel_to_y(pixel):
    return ((pixel)/map_length - 0.5) * 7

def phi_to_pixel(phi):
    return round((phi + 0.5*math.pi)/math.pi*1023)

def lambdda_to_pixel(lambdda):
    return round((lambdda + math.pi)/(2*math.pi)*2047)

def plot_parallels():
    parallels = np.linspace(-0.5*math.pi, 0.5*math.pi, num=13) #These have a certain lattitude phi, while longitude changes
    lambddas = np.linspace(-math.pi, math.pi, map_width)
    xy = np.empty((13, map_width, 2))
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
    meridians = np.linspace(-math.pi, math.pi, num=13)  # These have a certain longitude lambdda, while latitude changes
    phis = np.linspace(-0.5* math.pi, 0.5* math.pi, map_length)
    xy_meridians = np.empty((13, map_length, 2))
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

angle_dist = np.empty(  ( map_length , map_width , 1 ))
area_dist = np.empty(  ( map_length , map_width , 1 ))
scale_factors = np.empty((map_length , map_width, 2))
T_matrices = np.empty((map_length , map_width, 2, 2))

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
            s, T = get_scale_factors_at(phi, lambdda)
            scale_factors[y_pixel][x_pixel] = s
            T_matrices[y_pixel][x_pixel] = T
            omega = 2*math.asin(abs(s[0] - s[1])/(s[0] + s[1]))
            angle_dist[y_pixel][x_pixel] = omega
            area_dist[y_pixel][x_pixel] = s[0]*s[1]

def avg_angle_dist():
    total = 0
    total_normal = 0
    steps = 1000
    lambddas = np.linspace(-math.pi, math.pi, num=steps)
    phis = np.linspace(-math.pi/2, math.pi/2, num=steps*2)
    for phi in phis:
        for lambdda in lambddas:
            s, T = get_scale_factors_at(phi, lambdda)
            omega = 2*math.asin(abs(s[0] - s[1])/(s[0] + s[1]))
            total  = total + omega*math.cos(phi)
            total_normal = total_normal + math.cos(phi)
    return total/(total_normal)

# total = 0
# total_area = 0
# for phi_pixel in range(0, 1023):
#     for lambdda_pixel in range(0, 2047):
#         phi = phi_pixel/1023*math.pi - 0.5*math.pi
#         lambdda = lambdda_pixel/2047*2*math.pi - math.pi
#         s, T = get_scale_factors_at(phi, lambdda)
#         omega = 2*math.asin(abs(s[0] - s[1])/(s[0] + s[1]))
#         total  = total + omega
#         s0 = s[0]
#         s1 = s[1]
#         factor = s0*s1
#         if (factor != 0):
#             total_area = total_area + math.log(factor)
# avg = total/(1023*2047)
# avg_area = total_area/(1023*2047)

def score(c):
    total = 0
    total_normal = 0
    steps = 1000
    lambddas = np.linspace(-math.pi, math.pi, num=steps*2)
    phis = np.linspace(-math.pi/2, math.pi/2, num=steps)
    for phi in phis:
        for lambdda in lambddas:
            s, T = get_scale_factors_at(phi, lambdda)
            s[0] = c*s[0]
            s[1] = c*s[1]
            if (s[0] <= 0 or s[1] <= 0):
                continue
            # ar = s[0]*s[1]
            # if (s[0] > 0 and s[1] > 0 and ar < 1):
            #     ar = ar
            value = abs(math.log(s[0])) + abs(math.log(s[1]))
            total  = total + value*math.cos(phi)
            total_normal = total_normal + math.cos(phi)
    return total/(total_normal)

def score_1():
    total = 0
    total_normal = 0
    steps = 1000
    lambddas = np.linspace(-math.pi, math.pi, num=steps*2)
    phis = np.linspace(-math.pi/2, math.pi/2, num=steps)
    for phi in phis:
        for lambdda in lambddas:
            s, T = get_scale_factors_at(phi, lambdda)
            if (s[0] <= 0 or s[1] <= 0):
                continue
            ar = s[0]*s[1]
            if (s[0] > 0 and s[1] > 0 and ar < 1):
                ar = ar
            value = abs(math.log(s[0])) + abs(math.log(s[1]))
            total  = total + value*math.cos(phi)
            total_normal = total_normal + math.cos(phi)
    return total/(total_normal)
 
# score = score()

def area_score(c):
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
            value = abs(math.log(c*c*s[0]*s[1]))
            total  = total + value*math.cos(phi)
            total_normal = total_normal + math.cos(phi)
    return total/(total_normal)

# area_score = area_score()

plt.figure(dpi = 400)
plt.axis('off')
plot_parallels()
plot_meridians()
plt.imshow(new_world)

hm = plt.imshow(angle_dist, cmap="hot", interpolation='nearest', alpha=0.7, vmax = 1)  
ax = plt.gca()
cbar = ax.figure.colorbar(hm, ax=ax)

# plt.figure(dpi = 400)
# plt.axis('off')
# plot_parallels()
# plot_meridians()
# plt.imshow(new_world)

# hm = plt.imshow(area_dist, cmap="RdPu", interpolation='nearest', alpha=0.7, norm=colors.LogNorm(vmin = 1, vmax = 10**3))  
# ax = plt.gca()
# cbar = ax.figure.colorbar(hm, ax=ax)