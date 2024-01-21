import matplotlib.image as image
import matplotlib.pyplot as plt
import matplotlib.colors as colors

import numpy as np
import math

map_width = 2000
map_length = 2000
h = 0.001

phi_1 = math.acos(2/math.pi)

def transform(phi, lambdda):
    a = math.acos(math.cos(phi)*math.cos(lambdda/2))
    sinc_a = 1
    if (a != 0): sinc_a = math.sin(a)/a
    x = 1/2*(lambdda*math.cos(phi_1) + 2*math.cos(phi)*math.sin(lambdda/2)/sinc_a)
    y = 1/2*(phi + math.sin(phi)/sinc_a)
    return (y,x)


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
            value = abs(math.log(s[0])) + abs(math.log(s[1]))
            total  = total + value*math.cos(phi)
            total_normal = total_normal + math.cos(phi)
    return total/(total_normal)