# -*- coding: utf-8 -*-
"""
Created on Thu Sep  3 15:36:10 2020

@author: 용하
"""


import numpy as np
from scipy import interpolate
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import threading
from scipy.linalg import norm

##-----------패러미터--------------------------------------------------+
s=0.1 #폐곡선 표면부터 극판까지의 거리 
d=0.1 #극판 두께
a=0.5 #극판의 반지름 
div=100 #옆면을 나눈 수
div_z=50 #축을 나눈 수 

#----------그래프 조정-----------------------------------------------------+
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")

#-------------몸 곡면(구)-----------------------------------------+

def body_surface (radius, center_x, center_y, center_z): #반지름, 구의 중심 좌표 
    theta = np.linspace(0, 2 * np.pi, div)
    z = np.linspace(center_z-radius, center_z+radius, div_z)
    radi=[]
    z_array=[]
    for i in range(len(z)):
        r=(radius**2-abs(z[i]-center_z)**2)**0.5
        radi.append(r)
        z_array.append(np.array([z[i]]))
    
    radi=np.array(radi)
    z_array=np.array(z_array)
    x = np.outer(radi, np.cos(theta))+center_x
    y = np.outer(radi, np.sin(theta))+center_y
    ax.plot_surface(x, y, z_array, color='aqua',alpha=1, shade=2)


body_surface(10,1,1,1)
