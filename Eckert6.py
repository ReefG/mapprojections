import matplotlib.image as image
import matplotlib.pyplot as plt
import numpy as np
import math

map_width = 3000
map_length = 1500

h = 0.0001

def get_x(phi, lambdda, theta):
    return lambdda*(1+math.cos(theta))/math.sqrt(2+math.pi)

def get_y(phi, lambdda, theta):
    return 2*theta/math.sqrt(2+math.pi)

def get_theta(phi):
    theta1 = phi
    delta_theta1 = 10000000000000
    while (abs(delta_theta1) > 0.01):
        delta_theta1 = -(theta1 + math.sin(theta1)-(1 + 0.5* math.pi)*math.sin(phi))/(1+math.cos(theta1))
        theta1 = theta1 + delta_theta1
    theta =  theta1
    return theta

def transform(phi, lambdda):
    theta = get_theta(phi)
    return (get_y(phi, lambdda, theta), get_x(phi, lambdda, theta))

def inv_transform(x,y):
    theta_inv = 0.5*math.sqrt(2+math.pi)*y
    phi = math.asin((theta_inv + math.sin(theta_inv))/(1+0.5*math.pi))
    lambdda = (x*math.sqrt(2 + math.pi))/(1 + math.cos(theta_inv))
    return (phi, lambdda)


def get_scale_factors_at(phi, lambdda):
    theta_minush = get_theta(phi - h)
    theta_h = get_theta(phi + h)
    theta = get_theta(phi)
    
    x_phi = ( get_x(phi + h, lambdda, theta_h) - get_x(phi - h, lambdda, theta_minush) ) / (2*h)
    x_lambdda = (get_x(phi, lambdda + h, theta) - get_x(phi, lambdda - h, theta) ) / (2*h)
    y_phi = ( get_y(phi + h, lambdda, theta_h) -get_y(phi - h, lambdda, theta_minush) ) / (2*h)
    y_lambdda = (get_y(phi, lambdda + h, theta) - get_y(phi, lambdda - h, theta) ) / (2*h)
    
    T = [[x_lambdda/(math.cos(phi)), x_phi], [y_lambdda/math.cos(phi), y_phi]]
    
    s = np.linalg.svd(T, compute_uv=False)
    return s, T

def x_to_pixel(x):
    return round((x/(4*math.pi/math.sqrt(2+math.pi)) + 0.5)*map_width)

def pixel_to_x(pixel):
    return (pixel/map_width - 0.5) * 4*math.pi/math.sqrt(2+math.pi)

def y_to_pixel(y):
    return round((y/(2* math.pi/math.sqrt(2+math.pi)) + 0.5)*map_length)

def pixel_to_y(pixel):
    return ((pixel)/map_length - 0.5) * 2 * math.pi/math.sqrt(2+math.pi)

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
            
def score(c):
    total = 0
    total_normal = 0
    steps = 1000
    lambddas = np.linspace(-math.pi, math.pi, num=steps*2)
    phis = np.linspace(-math.pi/2, math.pi/2, num=steps)
    for phi in phis:
        for lambdda in lambddas:
            s, T = get_scale_factors_at(phi, lambdda)
            value = abs(math.log(c*s[0])) + abs(math.log(c*s[1]))
            total  = total + value*math.cos(phi)
            total_normal = total_normal + math.cos(phi)
    return total/(total_normal)
 
# score = score()

# avg = avg_angle_dist()
plt.figure(dpi = 400)
plt.axis('off')
plot_parallels()
plot_meridians()
plt.imshow(new_world)

hm = plt.imshow(angle_dist, cmap="hot", interpolation='nearest', alpha=0.7)  
ax = plt.gca()
cbar = ax.figure.colorbar(hm, ax=ax)