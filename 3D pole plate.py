# -*- coding: utf-8 -*-
"""
Created on Mon Aug  3 11:14:30 2020

@author: 용하
"""

#------------모듈-------------------------------------------+
import numpy as np
from scipy import interpolate
from matplotlib import pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import threading
from scipy.linalg import norm
import stl
from stl import mesh
import matplotlib.tri as mtri
import time
import math 
##-----------패러미터--------------------------------------------------+
s=1 #폐곡선 표면부터 극판까지의 거리 
d_1=1 #세라믹 두께
d_2=1 #극판 두께
a=2 #극판의 반지름 
div=75 #옆면을 나눈 수
div_z=75 #축을 나눈 수 
div_c=150 #원판의 둘레를 나눈 수 
r_1=10 # skin
r_2=9 # bone 
r_3=8 # CSF
r_4=7 # grey matter
r_5=6 # white matter
now=time.localtime(time.time())
#----------그래프 조정-----------------------------------------------------+
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.set_xlabel("X")
ax.set_ylabel("Y")
ax.set_zlabel("Z")
plt.xlim(-10,10) #x 축 범위
plt.ylim(-10,10) 
ax.set_zlim(-10,10)

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
    ax.plot_surface(x, y, z_array, color='yellow',alpha=0.3)
    return x,y,z_array
##----------중요!!------------------------------------------------+
    
xi, yi, zi = body_surface(r_1,0,0,0) # xi, yi, zi 는 이 전체 알고리즘에서 surface의 array를 의미한다. 
'''ax.plot_trisurf(xi, yi, zi, linewidth=0.2, antialiased=True)  '''  
xi2,yi2,zi2 = body_surface(r_2,0,0,0)
xi3,yi3,zi3 = body_surface(r_3,0,0,0)
xi4,yi4,zi4 = body_surface(r_4,0,0,0)
xi5,yi5,zi5 = body_surface(r_5,0,0,0)

#------------폐곡선 라인 벡터 만들기 (반시계 방향)--------------------+
def line_vector(x,y,z): #x,y,z 는 모두 array 이다
    list_v=[]
    for i in range(len(z)): #len(z)=div_z
        vector=[]
        cur_v=[]
        for j in range(len(x[i])): #len(x[i])=div
            vector.append(np.array([xi[i][j],yi[i][j]])) #위치벡터
        for j in range(len(x[i])): #len(x[i])=div
            if j<div-1:
                cur_v.append(vector[j+1]-vector[j])
            
            else:
                cur_v.append(vector[0]-vector[j])
        list_v.append(cur_v)
    list_v=np.array(list_v)
    return list_v
list_v=line_vector(xi,yi,zi) #3차원 array 이다 div_z*div*2 개 짜리

#------------좌표추출-------------------------------------------+
def cor(p):
    if p.count('x')>=1:           
            s=p.split(',')
            b=[]
            for i in range(len(s)):
                b.append(s[i].strip())
            c=[]
            for i in range(len(b)):
                c.append(b[i][2:])
            cx=float(c[0])
            cy=float(c[1])
            cz=float(c[2])
            i,j=distance(cx,cy,cz)
            kx=xi[i][j]
            ky=yi[i][j]
            kz=zi[i][0]
            return kx,ky,kz
    else :
            s=p.split(',')
            b=[]
            for i in range(len(s)):
                b.append(s[i].rstrip('ged '))
            d=[]
            d.insert(0,b[0][8:])
            d.insert(1,b[1][11:])
            print(d)

#-------------거리가까운 메시의 좌표특정을 위한 거 ---------------------+
def distance2(x,y,z):
    d_1=[]
    for i in range(len(zi)):
        c=abs(z-zi[i])
        d_1.append(c)
    e_1=d_1.index(min(d_1))
    d_2=[]
    for i in range(len(xi[e_1])):
        c=(x-xi[e_1][i])**2+(y-yi[e_1][i])**2
        d_2.append(c)
    e_2=d_2.index(min(d_2))
    return e_1, e_2  #z가 몇번째인지 추출하고, x,y를 추출 하는 것이다  
'''리턴값을 e로해서, 좌표자체를 특정하는게 아니라 좌표가 리스트에서 존재하는 
위치를 특정하는거임 여기서 e_1은 z의 좌표 e_2 x,y 의 좌표 '''

def distance(x,y,z): #진짜 최단거리 찾는 방법 
    d_1=[]
    for i in range(div_z):
        for j in range(div):
            c=(x-xi[i][j])**2+(y-yi[i][j])**2+(z-zi[i][0])**2
            d_1.append(c)
    e=d_1.index(min(d_1))
    e1=math.floor(e/div)
    e2=e-div*e1
    return e1, e2  #z가 몇번째인지 추출하고, x,y를 추출 하는 것이다  
'''리턴값을 e로해서, 좌표자체를 특정하는게 아니라 좌표가 리스트에서 존재하는 
위치를 특정하는거임 여기서 e_1은 z의 좌표 e_2 x,y 의 좌표 '''

#-----------접선그리기--------------------------------------------+
def tan_yx(x,y,z): #dy/dx 
    i, j=distance(x,y,z)
    if 0 < j <div-1:
        m1=(yi[i][j+1]-yi[i][j])/(xi[i][j+1]-xi[i][j])
        m2=(yi[i][j]-yi[i][j-1])/(xi[i][j]-xi[i][j-1])
        m=(m1+m2)/2
        return m
    elif j == 0 :
        m1=(yi[i][1]-yi[i][0])/(xi[i][1]-xi[i][0])
        m2=(yi[i][0]-yi[i][-1])/(xi[i][0]-xi[i][-1])
        m=(m1 + m2)/2
        return m
    elif j == div-1 : #마지막 점의 경우 k+1 번째가 없기에 에러가 난다. 때문에 써준거 
        m1=(yi[i][0]-yi[i][j])/(xi[i][0]-xi[i][j])
        m2=(yi[i][j]-yi[i][j-1])/(xi[i][j]-xi[i][j-1])
        m=(m1 + m2)/2
        return m 

def tan_xz(x,y,z): # dx/dz 
    i, j=distance(x,y,z)
    if 0 < i <div_z-1:
        m1=(xi[i+1][j]-xi[i][j])/(zi[i+1]-zi[i])
        m2=(xi[i][j]-xi[i-1][j])/(zi[i]-zi[i-1])
        m=(m1+m2)/2
        return m
    elif i == 0 :
        m1=(xi[1][j]-xi[0][j])/(zi[1]-zi[0])
        m2=(xi[0][j]-xi[-1][j])/(zi[0]-zi[-1])
        m=(m1 + m2)/2
        return m
    elif i == div_z-1 : #마지막 점의 경우 k+1 번째가 없기에 에러가 난다. 때문에 써준거 
        m1=(xi[0][j]-xi[i][j])/(zi[0]-zi[i])
        m2=(xi[i][j]-xi[i-1][j])/(zi[i]-zi[i-1])
        m=(m1 + m2)/2
        return m 
    
def tan_yz(x,y,z): # dy/dz 
    i, j=distance(x,y,z)
    if 0 < i <div_z-1:
        m1=(yi[i+1][j]-yi[i][j])/(zi[i+1]-zi[i])
        m2=(yi[i][j]-yi[i-1][j])/(zi[i]-zi[i-1])
        m=(m1+m2)/2
        return m
    elif i == 0 :
        m1=(yi[1][j]-yi[i][j])/(zi[1]-zi[i])
        m2=(yi[i][j]-yi[-1][j])/(zi[i]-zi[-1])
        m=(m1 + m2)/2
        return m
    elif i == div_z-1 : #마지막 점의 경우 k+1 번째가 없기에 에러가 난다. 때문에 써준거 
        m1=(yi[0][j]-yi[i][j])/(zi[0]-zi[i])
        m2=(yi[i][j]-yi[i-1][j])/(zi[i]-zi[i-1])
        m=(m1 + m2)/2
        return m 

#----------직선 함수-----------------------------------------------+
def line_1(m,x,y,z,b,color):
    ''' m 은 기울기, b 는 중심점부터 끝 까지의 길이
    좌우 끝점의 좌표를 리턴한다'''
    x_1=x-(b)*(1+m**2)**-0.5
    x_2=x+(b)*(1+m**2)**-0.5
    y_1= line(m,x_1,x,y)
    y_2= line(m,x_2,x,y)
    plt.plot([x_1,x_2],[y_1,y_2],[z,z], color)
    return x_1, y_1, z, x_2, y_2,z 
#------------라인 함수--------------------------------------------+ 
def line(m,x, xt, yt):
    return m*(x-xt)+yt   

#-----------접선벡터----------------------------------------------+
def vector_xy(x,y,z):
    i, j=distance(x,y,z)
    if 0 < j <div-1:
        v1=np.array([xi[i][j+1]-xi[i][j],yi[i][j+1]-yi[i][j],0])
        v2=np.array([xi[i][j]-xi[i][j-1],yi[i][j]-yi[i][j-1],0])
        v=(v1+v2)/2
        if norm(v) == 0:
            return np.array([0,1,0])
        else :
            return v
    elif j == 0 :
        v1=np.array([xi[i][1]-xi[i][0],yi[i][1]-yi[i][0],0])
        v2=np.array([xi[i][0]-xi[i][-1],yi[i][0]-yi[i][-1],0])
        v=(v1 + v2)/2
        if norm(v) == 0:
            return np.array([0,1,0])
        else :
            return v
    elif j == div-1 : #마지막 점의 경우 k+1 번째가 없기에 에러가 난다. 때문에 써준거 
        v1=np.array([xi[i][0]-xi[i][j],yi[i][0]-yi[i][j],0])
        v2=np.array([xi[i][j]-xi[i][j-1],yi[i][j]-yi[i][j-1],0])
        v=(v1 + v2)/2
        if norm(v) == 0:
            return np.array([0,1,0])
        else :
            return v

def vector_z(x,y,z): #나중에 추가 수정해야 할 부분으로 보인다. z=0, div_z 에서 포인트가 중첩되는것 때매 에러가 나는건데 이걸 궂이 해결할 필요가 있나 싶기도 하고.. 
    i, j=distance(x,y,z)
    if 0 < i <div_z-1:
        v1=np.array([xi[i+1][j]-xi[i][j],yi[i+1][j]-yi[i][j],zi[i+1][0]-zi[i][0]])
        v2=np.array([xi[i][j]-xi[i-1][j],yi[i][j]-yi[i-1][j],zi[i][0]-zi[i-1][0]])
        v=(v1+v2)/2
        return v
    elif i == 0 :
        v1=np.array([xi[1][j]-xi[0][j],yi[1][j]-yi[0][j],zi[1][0]-zi[0][0]])
        v=v1
        return np.array([1,0,0])
    elif i == div_z-1 : #마지막 점의 경우 k+1 번째가 없기에 에러가 난다. 때문에 써준거 
        v2=np.array([xi[i][j]-xi[i-1][j],yi[i][j]-yi[i-1][j],zi[i][0]-zi[i-1][0]])
        v=v2
        return np.array([-1,0,0])
    
#-----------노말 단위벡터--------------------------------------+
def n_vector(x,y,z):
    v_xy=vector_xy(x,y,z)
    v_z=vector_z(x,y,z)
    n_v=np.cross(v_xy, v_z)
    mag = norm(n_v)
    n_v=n_v/mag
    return n_v


#-----------평평한 극판---------------------------------------+
def cylinder(x,y,z,s,d_1,d_2,r): 
    '''x,y,z 는 폐곡선 위의 좌표, s는 아랫면 까지의 거리, d_1은 아래쪽 세라믹 두께, d_2는 위쪽 극판두께, r은 반지름'''
    t = np.array([0,d_1,d_1+d_2]) # 층을 3개로 강제로 나눈다
    theta = np.linspace(0, 2 * np.pi, div_c) #2pi 를 50개로 분해 
    radi = np.linspace(0, r, 2) #반지름을 중점과 끝점으로만 분해 
    v=n_vector(x,y,z) #법선벡터
    n1=vector_xy(x,y,z)
    
    mag1 = norm(n1)
    n1=n1/mag1 #평면단위벡터 1
    
    n2=np.cross(v,n1) 
    mag2=norm(n2)
    n2=n2/mag2 #평면단위벡터 2
    p0=np.array([x,y,z])+s*v # 밑면 중심의 위치 벡터
   

    #use meshgrid to make 2d arrays
    radi,theta1 = np.meshgrid(radi, theta)
    t, theta2 = np.meshgrid(t, theta)

    #generate coordinates for surface
    # "Tube"
    X1, Y1, Z1 = [p0[i] + v[i] * t + r * np.sin(theta2) * n1[i] + r * np.cos(theta2) * n2[i]for i in [0, 1, 2]]
    # "Bottom"
    X2, Y2, Z2 = [p0[i] + radi[i] * np.sin(theta1) * n1[i] + radi[i] * np.cos(theta1) * n2[i]for i in [0, 1, 2]]
    # "Middle"
    X3, Y3, Z3 = [p0[i] + v[i]*d_1 + radi[i] * np.sin(theta1) * n1[i] + radi[i] * np.cos(theta1) * n2[i]for i in [0, 1, 2]]
    # "Top"
    X4, Y4, Z4 = [p0[i] + v[i]*(d_1+d_2) + radi[i] * np.sin(theta1) * n1[i] + radi[i] * np.cos(theta1) * n2[i]for i in [0, 1, 2]]
    

    ax.plot_surface(X1, Y1, Z1, color='red', alpha=0.5)
    ax.plot_surface(X2, Y2, Z2, color='red', alpha=0.5)
    ax.plot_surface(X3, Y3, Z3, color='red', alpha=0.5)
    ax.plot_surface(X4, Y4, Z4, color='red', alpha=0.5)
    return X1, Y1, Z1, X2, Y2, Z2, X3, Y3, Z3, X4, Y4, Z4 # 모든좌표를 다 추출한것

#----------------------------------
def pick_p(a): #원판 둘레의 포인트를 뽑기위한 함수
    b=[]
    for i in range(len(a)):
        b.append(a[i][1])
    b=np.array(b)
    return b

def pick_face(x,y,z,v,s): #좌표에 해당하는 facet을 뽑는 작업 v는 해당하는 노말 벡터 
    '''x,y,z 는 아랫면의 좌표, v는 법선 벡터, s는 폐곡면으로부터 아랫면까지의 거리 '''
    pz=z-s*v[2] #사실 이부분이 논리상 문제가 조금 있기는 하다. s를 확정 지을수 없기 때문에.. 추후 수정필요 
    i=0
    while pz >= zi[i][0]: # z 축 레이어 특정, 내가 원하는 건 i-1 이다. 
        i+=1
        if i==div_z-1:
            break
    v_0=np.array([x-xi[i-1][0], y-yi[i-1][0], z-zi[i-1][0]])
    v_2=np.cross(v_0,v)
    if not any([v[0], v[1]]) == 0:
        
        if v_2[2] < 0 : 
            j=0
            while v_2[2] < 0: #xy 쪽을 추출하는 거다. 내가 원하는 건 j-1이다.  
                v_1=np.array([x-xi[i-1][j+1], y-yi[i-1][j+1], z-zi[i-1][0]])
                v_2=np.cross(v_1,v)
                j+=1
                if j==div-1:
                    j-=1
                    break
    
        else :
            j=0
            while v_2[2] >= 0:
                v_1=np.array([x-xi[i-1][j-2], y-yi[i-1][j-2], z-zi[i-1][0]])
                v_2=np.cross(v_1,v)
                j-=1
                if j==-1*(div-1):
                    j+=1
                    break
    
    else: 
        dis_list=[]
        for k in range(div):
            distance=(x-xi[i-1][k])**2+(y-yi[i-1][k])**2
            dis_list.append(distance)
        m1=dis_list.index(min(dis_list))
        del dis_list[m1]
        m2=dis_list.index(min(dis_list))
        if 0<m1<div-1: 
            if m1==m2 :
                j=m1+1
            elif m1>m2 :
                j=m1
        elif m1==0 :
            if m1==m2 :
                j=m1+1
            else:
                j=m1
        elif m1==div-1 :
            if m2==0:
                j=0
            else:
                j=div-1
    
    
    v_0=np.array([x-xi[i][0], y-yi[i][0], z-zi[i][0]])
    v_2=np.cross(v_0,v)
    if not any([v[0], v[1]]) == 0: #법선 벡터가[0,0,-1] 인 경우에는 버그가 난다. 그것을 수정한 것 
        
        if v_2[2] < 0 : 
            j1=0
            while v_2[2] < 0: #xy 쪽을 추출하는 거다. 내가 원하는 건 j-1이다. 
                v_1=np.array([x-xi[i][j1+1], y-yi[i][j1+1], z-zi[i][0]])
                v_2=np.cross(v_1,v)
                j1+=1
                if j1==div-1:
                    j1-=1
                    break
                
        else :
            j1=0
            while v_2[2] >= 0:
                v_1=np.array([x-xi[i][j1-2], y-yi[i][j1-2], z-zi[i][0]])
                v_2=np.cross(v_1,v)
                j1-=1
                if j1==-1*(div-1):
                    j1+=1
                    break
        
    else: 
        dis_list=[]
        for k in range(div):
            distance=(x-xi[i][k])**2+(y-yi[i][k])**2
            dis_list.append(distance)
        m1=dis_list.index(min(dis_list))
        del dis_list[m1]
        m2=dis_list.index(min(dis_list))
        if 0<m1<div-1: 
            if m1==m2 :
                j1=m1+1
            elif m1>m2 :
                j1=m1
        elif m1==0 :
            if m1==m2 :
                j1=m1+1
            else:
                j1=m1
        elif m1==div-1 :
            if m2==0:
                j1=0
            else:
                j1=div-1

    return i, j, j1 # i는 z축 넘버링, j, j1은 각각 아랫쪽 윗쪽 x,y 넘버링 
    
def face(p1,p2,p3,p4) : #점 4개로 평면에 방정식에 해당하는 a,b,c,d 를 구하는 것 
    p1=list(p1)
    p2=list(p2)
    if not p1==p2:
        p1=np.array(p1)
        p2=np.array(p2)
        v_1=p2-p1
    else : 
        p1=np.array(p1)
        p2=np.array(p2)
        v_1=p4-p2    
    v_2=p3-p2
    v_n=np.cross(v_1, v_2)
    a=v_n[0]
    b=v_n[1]
    c=v_n[2]
    d=v_n[0]*p2[0]+v_n[1]*p2[1]+v_n[2]*p2[2]
    return a, b, c, d
    

def cross_point(v,a,b,c,d,x,y,z): #x,y,z 는 원판위에 좌표, v 는 법선 벡터
    h1=abs(a*x+b*y+c*z-d)/(a**2+b**2+c**2)**0.5 #평면과 원판위 점 사이의 거리
    cos=abs(np.dot(np.array([a,b,c]),v))/(a**2+b**2+c**2)**0.5
    h=h1/cos
    p0=np.array([x,y,z])
    p1=p0-h*v
    return p1
#------------3d flat plate---------------------------------------+
def flat_plate(x,y,z,s,d_1,d_2,r):
    '''x,y,z 는 폐곡선 위의 찍은 점의 좌표, s는 찍은 점부터 극판 까지의 거리
    d 는 극판 두께, r은 극판 반지름 '''
    e_z, e =distance(x,y,z)
    x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4= cylinder(x,y,z,s,d_1,d_2,r) #4가 세라믹의 좌표임 
    x_list=pick_p(x2)
    y_list=pick_p(y2)
    z_list=pick_p(z2)
    v=n_vector(x,y,z)
    p=[] 
    for i in range(len(x_list)):
        p.append(np.array([x_list[i],y_list[i],z_list[i]]))
    p=np.array(p) # 아랫면의 둘레에 있는 위치벡터들의 리스트 
    ap_list=[] #폐곡선 위의 점들의 리스트를 만들 것이다 
    info_list = [] # 이후의 함수를 위한 것 
    facet_num = []
    for i in range(len(p)):
        t0,t1,t2 =pick_face(p[i][0],p[i][1],p[i][2],v,s)
          #위에서 구한 메쉬의 4개의 포인트 
        p1=np.array([xi[t0-1][t1-1],yi[t0-1][t1-1],zi[t0-1][0]])
        p2=np.array([xi[t0-1][t1],yi[t0-1][t1],zi[t0-1][0]])
        p3=np.array([xi[t0][t2],yi[t0][t2],zi[t0][0]])
        p4=np.array([xi[t0][t2-1],yi[t0][t2-1],zi[t0][0]])
        a,b,c,d=face(p1,p2,p3,p4)
        ap=cross_point(v,a,b,c,d,p[i][0],p[i][1],p[i][2]) 
        ap_list.append(ap) 
        info=[t0,t1,t2,ap[2]] #축 순서, 큰쪽 작은쪽 xy 순서, z 좌표
        info_list.append(info)
        facet_num.append([t0,t2]) #FACET 3번 포인트의 노드 넘버를 뽑는다. 
    
    ap_list=np.array(ap_list)  #폐곡선에 크로스되는 점들의 xyz 어레이
    
    side_face(x_list,y_list,z_list,ap_list) #하이드로겔의 옆면을 만든것. 아웃풋 해야함 
 #   cross_point3(ap_list, facet_num, p) # 구 face 위의 grid 와 크로스되는 지점들을 구하는것 
    
    info2_list=[]
    for j in range(int(len(info_list)/2)): 
        info2=[info_list[j][0],min(info_list[j][1],info_list[j][2]),
               max(info_list[div_c-1-j][1],info_list[div_c-1-j][2]),info_list[j][3]]
        info2_list.append(info2)
    
    base_sur=[]
    
    for k in range(len(info2_list)):
        cro_point=cross_point2(info2_list[k][0],info2_list[k][1],info2_list[k][2],info2_list[k][3])
        cro_point.insert(0,ap_list[div_c-1-k])
        cro_point.append(ap_list[k])
        cro_point=np.array(cro_point) #원판의 시작점과 끝점, 중간에 교차되는 점들을 모두 구한 array다
        base_sur.append(cro_point)
             
    base_sur=np.array(base_sur)
    base_face(base_sur) #아랫면을 만든것 아웃풋 해야함 
    
    output2(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,ap_list)
    near_point(ap_list)
    num=delit (facet)
    output3(x,y,z,v)
    print(len(dic))
    print(len(facet))
    export1(dic,facet,num)

def delit (facet):
    many_list=[]
    for i in range(len(facet)):
        many_list.append(i+1)    
    
    prb=[]    
    for i in range(len(facet)):
        if facet[i+1][1]==3:
            info=[]
            info.append(facet[i+1][2])
            info.append(facet[i+1][3])
            info.append(facet[i+1][4])
            info2=set(info)
            info2=list(info2)
            if not  len(info)==len(info2):
                prb.append(i+1)
            
        elif facet[i+1][1]==4:
            info=[]
            info.append(facet[i+1][2])
            info.append(facet[i+1][3])
            info.append(facet[i+1][4])
            info.append(facet[i+1][5])
            info2=set(info)
            info2=list(info2)
            if not  len(info)==len(info2):
                prb.append(i+1)
    

    a_sub_b = [x for x in many_list if x not in prb]
    return a_sub_b

#-----------니어 포인트---------------------------------------------------#
def near_point(ap_list):
    many=len(dic)
    many2=len(facet)
    point_num=[]
    for i in range(len(ap_list)):
        e1,e2=distance(ap_list[i][0],ap_list[i][1],ap_list[i][2])
        point_num.append(e1*div+e2+1)
    
    for i in range(len(point_num)):
        if i<len(point_num)-1: 
            facet[many2+2*i+1]=[-1,3,many-(3*div_c+3)+i+1,many-(3*div_c+3)+i+2,point_num[i]] #ap_list_num=many-(3*div_ci+2)+i+1
            facet[many2+2*i+2]=[-1,3,many-(3*div_c+3)+i+2,point_num[i],point_num[i+1]]
       
            
#-----------아웃풋2,3: 위에 서 구한 실린더의 노드와 폐곡선과 만나는 노드, 페이스들을 모조리 추출한것----+
def output2(x1,y1,z1,x2,y2,z2,x3,y3,z3,x4,y4,z4,ap_list):
    point=[]
        
    for i in range(div_c) :
        point.append([x1[i][0],y1[i][0],z1[i][0]])
    for i in range(div_c) : # 합치면 안된다. 노드의 순서가 달라짐 
        point.append([x1[i][1],y1[i][1],z1[i][1]])
    for i in range(div_c) : # 합치면 안된다. 노드의 순서가 달라짐 
        point.append([x1[i][2],y1[i][2],z1[i][2]])
        
    point.append([x2[0][0],y2[0][0],z2[0][0]])
    point.append([x3[0][0],y3[0][0],z3[0][0]])
    point.append([x4[0][0],y4[0][0],z4[0][0]])
    
    many=len(dic)
    for i in range(len(point)):
        dic[many+i+1]=point[i]
    
    many2=len(facet)
    facet_list=[]
    for i in range(2): #층이 세개여서 면으로 이루어진 층은 세개의 층이다
        for j in range(div_c): 
            if j < div_c-1: #j가 0~48까지, 즉 노드 넘버로는 1~49, 극판의 동그란 면 만드는 것 
                facet_list.append([-1,3,many+(i)*div_c+j+1,many+(i)*div_c+j+2,many+3*div_c+i+1])
    
    for j in range(div_c): 
            if j < div_c-1: #i=2 일때만 따로 만들어 준 것이다. 바운데리를 나눠야 해서 
                facet_list.append([click_num+1,3,many+(2)*div_c+j+1,many+(2)*div_c+j+2,many+3*div_c+2+1])
                
            
    for i in range(2): #
        for j in range(div_c): 
                if j < div_c-1: #j가 0~48까지, 즉 노드 넘버로는 1~49, 극판의 옆면만드는 것 
                    facet_list.append([-1,4,many+i*div_c+j+1,many+i*div_c+j+2,many+(i+1)*div_c+j+2,many+(i+1)*div_c+j+1])
            
    for i in range(len(facet_list)):
        facet[many2+i+1]=facet_list[i]
        
'''def cross_point3(ap_list, facet_num,p): #facet_num 은 각 ap_list가 속해있는 facet의 넘버, p는 극판 아랫면의 리스트  
    many=len(dic)
    many2=len(facet)
    for i in range(len(facet_num)):    
        if i<len(facet_num)-1:
            if facet_num[i][0]==facet_num[i+1][0]:
                if not facet_num[i][1]==facet_num[i+1][1]: # xy축 방향의 크로스 되는 점이 새롭게 생길 때
                    z_num=facet_num[i][0]
                    xy_num=min(facet_num[i][1],facet_num[i+1][1])
                    direction=0 #기본 오른쪽 방향
                    if xy_num==facet_num[i][1]: # 같으면 왼쪽 방향, 즉 i i+1순서로 커지는 경우에, 내가 구하고자 하는건 그 반대인 왼쪽 방향임 
                        direction=-1
                    v1=np.array([xi[z_num][xy_num],yi[z_num][xy_num],zi[z_num][0]])
                    v2=np.array([xi[z_num-1][xy_num],yi[z_num-1][xy_num],zi[z_num][0]])
                    v, t= cross_plane_vector(ap_list[i],ap_list[i+1],p[i],p[i+1],v1,v2)
                    dic[many+i+1]=v
                    del facet[div*(facet_num[i][0]-1)+facet_num[i][1]] #이거 바꿔야 할 수도 
                    del facet[div*(facet_num[i+1][0]-1)+facet_num[i+1][1]]
                    facet[many2+8*(i+1)-7]=[1,3,div*(z_num-1)+facet_num[i][1],
                                      div*(z_num-1)+facet_num[i][1]+1,many-(3*div_c+2)+i+1] #many-(3*div_c+2)+i+1는 ap_list[i]의 넘버링을 의미 
                    facet[many2+8*(i+1)-6]=[1,3,div*(z_num)+facet_num[i][1],div*(z_num)+facet_num[i][1]+1,many-(3*div_c+2)+i+1]
                    facet[many2+8*(i+1)-5]=[1,3,div*(z_num-1)+xy_num,many+i+1,many-(3*div_c+2)+i+1]
                    facet[many2+8*(i+1)-4]=[1,3,div*(z_num)+xy_num,many+i+1,many-(3*div_c+2)+i+1]
                    facet[many2+8*(i+1)-3]=[1,3,div*(z_num)+xy_num+1+direction,
                                      div*(z_num-1)+xy_num+1+direction,many-(3*div_c+2)+i+1]
                    facet[many2+8*(i+1)-2]=[3,3,many-(3*div_c+2)+i+1,many+i+1,many-(2*div_c+2)+i+1]
                    facet[many2+8*(i+1)-1]=[3,3,many-(3*div_c+2)+i+2,many+i+1,many-(2*div_c+2)+i+2]
                    facet[many2+8*(i+1)]=[3,3,many-(2*div_c+2)+i+1,many+i+1,many-(2*div_c+2)+i+2]
                    
                elif facet_num[i][1]==facet_num[i+1][1]:
                    facet[many2+i+1]=[2,4,many-(3*div_c+2)+i+1,many-(3*div_c+2)+i+2,
                                      many-(2*div_c+2)+i+1,many-(2*div_c+2)+i+2] #ap_list i, i+1, 아랫판 i, i+1
            
            elif not facet_num[i][0]==facet_num[i+1][0]:
                if facet_num[i][1]==facet_num[i+1][1]:
                    xy_num=facet_num[i][1]
                    z_num=min(facet_num[i][0],facet_num[i+1][0])
                    direction=0 #기본 위쪽 방향 
                    if z_num==facet_num[i][0]: # 같으면 아래쪽 방향, 즉 i i+1순서로 커지는 경우에, 내가 구하고자 하는건 그 반대인 아래쪽 방향임 
                        direction=-1
                    v1=np.array([xi[z_num][facet_num[i][1]],yi[z_num][facet_num[i][1]],zi[z_num][0]])
                    v2=np.array([xi[z_num][facet_num[i][1]-1],yi[z_num][facet_num[i][1]-1],zi[z_num][0]])
                    v, t = cross_plane_vector(ap_list[i],ap_list[i+1],p[i],p[i+1],v1,v2)
                    
                    dic[many+i+1]=v
                    del facet[div*(facet_num[i][0]-1)+facet_num[i][1]]
                    del facet[div*(facet_num[i+1][0]-1)+facet_num[i+1][1]]
                    facet[many2+8*(i+1)-7]=[1,3,div*(facet_num[i][0]-1)+facet_num[i][1],
                                      div*(facet_num[i][0])+facet_num[i][1],many-(3*div_c+2)+i+1] #many-(3*div_c+2)+i+1는 ap_list[i]의 넘버링을 의미 
                    facet[many2+8*(i+1)-6]=[1,3,div*(facet_num[i][0]-1)+facet_num[i][1]+1,
                                            div*(facet_num[i][0])+facet_num[i][1]+1,many-(3*div_c+2)+i+1]
                    facet[many2+8*(i+1)-5]=[1,3,div*(z_num)+xy_num,many+i+1,many-(3*div_c+2)+i+1]
                    facet[many2+8*(i+1)-4]=[1,3,div*(z_num)+xy_num+1,many+i+1,many-(3*div_c+2)+i+1]
                    facet[many2+8*(i+1)-3]=[1,3,div*(z_num+direction)+xy_num,
                                      div*(facet_num[i][0]+direction)+xy_num+1,many-(3*div_c+2)+i+1]
                    facet[many2+8*(i+1)-2]=[3,3,many-(3*div_c+2)+i+1,many+i+1,many-(2*div_c+2)+i+1]
                    facet[many2+8*(i+1)-1]=[3,3,many-(3*div_c+2)+i+2,many+i+1,many-(2*div_c+2)+i+2]
                    facet[many2+8*(i+1)]=[3,3,many-(2*div_c+2)+i+1,many+i+1,many-(2*div_c+2)+i+2]
                
                elif not facet_num[i][1]==facet_num[i+1][1]:
                    xy_num=min(facet_num[i][1],facet_num[i+1][1])   
                    z_num=min(facet_num[i][0],facet_num[i+1][0])
                    v1=np.array([xi[z_num][xy_num],yi[z_num][xy_num],zi[z_num][0]])
                    v2=np.array([xi[z_num][xy_num+1],yi[z_num][xy_num+1],zi[z_num][0]])
                    v3=np.array([xi[z_num+1][xy_num],yi[z_num+1][xy_num],zi[z_num][0]])
                    v4=np.array([xi[z_num][xy_num-1],yi[z_num][xy_num-1],zi[z_num][0]])
                    v5=np.array([xi[z_num-1][xy_num],yi[z_num-1][xy_num],zi[z_num][0]])
                    xy_direction=1
                    z_direction=1
                    if xy_num==facet_num[i][1]:
                        xy_direction=0
                    if z_num==facet_num[i][0]:
                        z_direction=0
                    
                    del facet[div*(facet_num[i][0]-1)+facet_num[i][1]]
                    del facet[div*(facet_num[i+1][0]-1)+facet_num[i+1][1]]
                    
                    c1, t1 = cross_plane_vector(ap_list[i],ap_list[i+1],p[i],p[i+1],v1,v2)
                    c2, t2 = cross_plane_vector(ap_list[i],ap_list[i+1],p[i],p[i+1],v1,v3)
                    c3, t3 = cross_plane_vector(ap_list[i],ap_list[i+1],p[i],p[i+1],v1,v4)
                    c4, t4 = cross_plane_vector(ap_list[i],ap_list[i+1],p[i],p[i+1],v1,v5)
                    
                    #중앙점의 위를 자나냐 아래를 지나냐를 기준으로 나눈것. 위면 1, 아래면 0 이다. 
                    if t3==0 and t4==0 : 
                        case=1
                    if t4==0 and t1==0 : 
                        case=1
                    if t1==0 and t2==0 : 
                        case=0
                    if t2==0 and t3==0 : 
                        case=0
                    
                    direction=[xy_direction,z_direction,case] # 총 8가지의 케이스가 존재하게 된다.
                    
                    facet[many2+8*(i+1)-7]=[1,3,div*(facet_num[i][0]-1)+facet_num[i][1],
                                      div*(facet_num[i][0])+facet_num[i][1],many-(3*div_c+2)+i+1] #many-(3*div_c+2)+i+1는 ap_list[i]의 넘버링을 의미 
                    facet[many2+8*(i+1)-6]=[1,3,div*(facet_num[i][0]-1)+facet_num[i][1]+1,
                                            div*(facet_num[i][0])+facet_num[i][1]+1,many-(3*div_c+2)+i+1]
                    facet[many2+8*(i+1)-5]=[1,3,div*(z_num)+xy_num,many+i+1,many-(3*div_c+2)+i+1]
                    facet[many2+8*(i+1)-4]=[1,3,div*(z_num)+xy_num+1,many+i+1,many-(3*div_c+2)+i+1]
                    facet[many2+8*(i+1)-3]=[1,3,div*(z_num+direction)+xy_num,
                                      div*(facet_num[i][0]+direction)+xy_num+1,many-(3*div_c+2)+i+1]'''
                    
                    
                    
                    
def cross_plane_vector(p1,p2,p3,p4,v1,v2):
    a,b,c,d=face(p1,p2,p3,p4)
    vector=v1-v2
    abc=np.array([a,b,c])
    t=-1*(np.dot(abc,v2)+d)/np.dot(abc,vector)
    if not 0<=t<=1:
        t==0
    v=v2+t*vector
    return v, t

def output3(x,y,z,v):
    point=np.array([x,y,z])
    point1=point+(s/2)*v
    point2=point+(s+d_1/2)*v
    point3=point+(s+d_1+d_2/2)*v
    region[(3*click_num+1)+5]=[point1[0],point1[1],point1[2],6]
    region[(3*click_num+2)+5]=[point2[0],point2[1],point2[2],7]
    region[(3*click_num+3)+5]=[point3[0],point3[1],point3[2],8]
    
    
def export1(dic,face,num):
    f=open('./export{}{}{}{}.poly'.format(now.tm_mon,now.tm_mday,now.tm_hour,now.tm_min),'w')
    f.write('# Part 1 - node list\n')
    f.write('# node count, 3 dim, no attribute, no boundary marker\n')
    f.write('{} 3 0 0\n'.format(len(dic)))
    f.write('# Node index, node coordinates\n')
    for i in range(len(dic)):
        f.write('{} {} {} {}\n'.format(i+1,dic[i+1][0],dic[i+1][1],dic[i+1][2]))
    f.write('# Part 2 - facet list\n')
    f.write('# facet count, have boundary marker\n')
    f.write('{} 1\n'.format(len(num)))
    f.write('# facets\n')
    for i in range(len(num)):
        if facet[num[i]][1]==4:
            f.write('1 0 {}\n{} {} {} {} {}\n'.format(facet[num[i]][0],facet[num[i]][1],facet[num[i]][2],facet[num[i]][3],facet[num[i]][4],facet[num[i]][5]))
        elif facet[num[i]][1]==3:
            f.write('1 0 {}\n{} {} {} {}\n'.format(facet[num[i]][0],facet[num[i]][1],facet[num[i]][2],facet[num[i]][3],facet[num[i]][4]))
    f.write('# Part 3 - hole list\n')
    f.write('0 # no hole\n')
    f.write('# Part 4 - region list\n')
    f.write('{} # number of region\n'.format(len(region)))
    for i in range(len(region)):
        f.write('{} {} {} {} {}\n'.format(i+1, region[i+1][0],region[i+1][1],region[i+1][2],region[i+1][3]))

    
def export2(dic,facet):
    f=open('./export{}{}{}{}.node'.format(now.tm_mon,now.tm_mday,now.tm_hour,now.tm_min),'w')
    f.write('{}  3  0  0\n'.format(len(dic)))
    for i in range(len(dic)):
        f.write('   {}    {}  {}  {}\n'.format(i,dic[i+1][0],dic[i+1][1],dic[i+1][2]))
#-------------------------------------------------------------------------------------+
def side_face(x_list,y_list,z_list,ap_list): #밑면이랑 폐곡선 크로스되는 지점 이용해서 메쉬 짠것 
    x_grid=[]
    y_grid=[]
    z_grid=[]
    x_grid.append(x_list)
    y_grid.append(y_list)
    z_grid.append(z_list)
    x2_list=[]
    y2_list=[]
    z2_list=[]
    for i in range(len(ap_list)):
        x2_list.append(ap_list[i][0])
        y2_list.append(ap_list[i][1])
        z2_list.append(ap_list[i][2])
    x2_list=np.array(x2_list)
    y2_list=np.array(y2_list)
    z2_list=np.array(z2_list)
    x_grid.append(x2_list)
    y_grid.append(y2_list)
    z_grid.append(z2_list)
    x_grid=np.array(x_grid)
    y_grid=np.array(y_grid)
    z_grid=np.array(z_grid)
    
    ax.plot_surface(x_grid, y_grid, z_grid, color='black', alpha=0.5)

def cross_point2(zt, ht, lt, z) : # 폐곡선 위의 원에 수평한 선이 그리드들과 만나는 점들 
    p2_list=[]
    j_list=list(range(lt,ht))
    
    for i in j_list:
        v1=np.array([xi[zt][i],yi[zt][i],zi[zt][0]])
        v2=np.array([xi[zt-1][i],yi[zt-1][i],zi[zt-1][0]]) #point 2 임
        vector=v1-v2
        k=(z-zi[zt-1][0])/vector[2]
        p2=v2+k*vector
        p2_list.append(p2)
        
    return p2_list      

def base_face(base_sur):
    x_list=[]
    y_list=[]
    z_list=[]
    for i in range(len(base_sur)):
        x_l=[]
        y_l=[]
        z_l=[]
        for j in range(len(base_sur[i])):
            x_l.append(base_sur[i][j][0])
            y_l.append(base_sur[i][j][1])
            z_l.append(base_sur[i][j][2])
        
        x_list.append(x_l)
        y_list.append(y_l)
        z_list.append(z_l)
    
    x_grid, y_grid, z_grid =fit2(x_list, y_list, z_list)

    ax.plot_surface(x_grid, y_grid, z_grid, color='blue', alpha=1)
#--------array의 갯수를 맞추기 위한 함수(최댓값으로 맞춘다)--------------------+
def fit2(a,b,c): #a는 여기서 x array list를, b는 y array list, c는 z array list를 의미한다
    how_many=[]
    re_a=[]
    re_b=[]
    re_c=[]
    for i in a:
        how_many.append(len(i))
    many= max(how_many)  # 가장 많은 어레이의 갯수를 찾고, 그 개수에 맞추어 나머지를 인터폴 레이트 함 

    for i in range(len(a)):
        f=interpolate.interp1d(a[i],b[i])
        a_new=np.linspace(a[i][0],a[i][-1],many,endpoint=True)
        re_a.append(a_new)
        re_b.append(f(a_new))
    grid_a=np.array(re_a)
    grid_b=np.array(re_b)
    
    for i in range(len(c)):
        re_c.append([c[i][0]])
    grid_c=np.array(re_c)
    
    return grid_a, grid_b, grid_c


##-----------------tri_mesh----------------------------------------------------------------+
#----------------body stl--------------------------------------------------+
'''triang=mtri.Triangulation(xi,yi,zi)
data = np.zeros(len(triang.triangles), dtype=mesh.Mesh.dtype)
mobius_mesh = mesh.Mesh(data, remove_empty_areas=False)
mobius_mesh.x[:] = xi[triang.triangles]
mobius_mesh.y[:] = yi[triang.triangles]
mobius_mesh.z[:] = zi[triang.triangles]
mobius_mesh.save('mysurface.stl')'''

def output1(xi,yi,zi): #body에 대한 아웃 노드와 페이스들 
    point=[]
    for i in range(div_z):
        for j in range(div):
            point.append(np.array([xi[i][j],yi[i][j],zi[i][0]]))
    point=np.array(point)

    numbering=[]
    for i in range(div_z):
        for j in range(div):
            numbering.append(i*div+j+1)
    dic={}
    for i in range(len(numbering)):
        dic[numbering[i]]=point[i]
            
    facet_list=[]
    facet={}
    for i in range(div_z-1):# 순서대로 바운더리 마커, 노드수, 노드 넘버 
        if i==0:
            for j in range(div):
                if j < div-1: #j가 0~98까지, 즉 노드 넘버로는 1~99
                    facet_list.append([-1,3,div*i+j+2,div*(i+1)+j+2,div*(i+1)+j+1]) #맨앞에 1은 바운데리 마커, 3은 노드 수  
                
        
        elif 0<i<div_z-2 :
            for j in range(div):
                if j < div-1: #j가 0~98까지, 즉 노드 넘버로는 1~99
                    facet_list.append([-1,4,div*i+j+1,div*i+j+2,div*(i+1)+j+2,div*(i+1)+j+1])
               
                    
        elif i==div_z-2:
            for j in range(div):
                if j < div-1: #j가 0~98까지, 즉 노드 넘버로는 1~99
                    facet_list.append([-1,3,div*i+j+1,div*i+j+2,div*(i+1)+j+2])
               
    
    for i in range(len(facet_list)):
        facet[i+1]=facet_list[i]
    
    region={}
    region[1]=[(r_1+r_2)/2,0,0,1]
    region[2]=[(r_2+r_3)/2,0,0,2]
    region[3]=[(r_3+r_4)/2,0,0,3]
    region[4]=[(r_4+r_5)/2,0,0,4]
    region[5]=[0,0,0,5]
    return dic, facet, region

def add_output2(xi2,yi2,zi2):
    many=len(dic)
    point=[]
    for i in range(div_z):
        for j in range(div):
            point.append(np.array([xi2[i][j],yi2[i][j],zi2[i][0]]))
    point=np.array(point)
    for i in range(len(point)):
        dic[many+i+1]=point[i]
             
    facet_list=[]
    
    for i in range(div_z-1):# 순서대로 바운더리 마커, 노드수, 노드 넘버 
        if i==0:
            for j in range(div):
                if j < div-1: #j가 0~98까지, 즉 노드 넘버로는 1~99
                    facet_list.append([-1,3,many+div*i+j+2,many+div*(i+1)+j+2,many+div*(i+1)+j+1]) #맨앞에 4은 바운데리 마커, 3은 노드 수  
        
        elif 0<i<div_z-2 :
            for j in range(div):
                if j < div-1: #j가 0~98까지, 즉 노드 넘버로는 1~99
                    facet_list.append([-1,4,many+div*i+j+1,many+div*i+j+2,many+div*(i+1)+j+2,many+div*(i+1)+j+1])
            
        elif i==div_z-2:
            for j in range(div):
                if j < div-1: #j가 0~98까지, 즉 노드 넘버로는 1~99
                    facet_list.append([-1,3,many+div*i+j+1,many+div*i+j+2,many+div*(i+1)+j+2])
                       
    many2=len(facet)
    for i in range(len(facet_list)):
        facet[many2+i+1]=facet_list[i]
        
#-----------먼저 몸의 노드와 페이스들을 모조리 뽑는다-------------+    
dic, facet, region =output1(xi,yi,zi)
add_output2(xi2,yi2,zi2)
add_output2(xi3,yi3,zi3)
add_output2(xi4,yi4,zi4)
add_output2(xi5,yi5,zi5)
print(len(dic))
print(len(facet))
#-----------노말벡터 그리기-----------------------------------+
def normal(m,lx,ly,z):
    ''' 아래판 좌표와 노말 벡터의 기울기를 입력하면 폐곡선 위의 좌표를 찾아서 출력 '''
    ex_l=[] #l(폐곡선 좌표)에서 선별된 친구들 
    d=[]
    f=[] #폐곡선 위의 점중에서, 엡실론 보다 작은 거리에 있는 친구들 
    h=[]
    for i in range(len(l)):
        v_1=np.array([lx-l[i][0],ly-l[i][1],0])
        v_2=np.array([cur_v[i][0],cur_v[i][1],0])
        cross=np.cross(v_1,v_2)
        if cross[2] > 0: #법선 벡터랑 폐곡선 벡터랑 양수인 값만 뽑아낸다 
            ex_l.append([l[i][0],l[i][1]])
   
    for i in range(len(ex_l)):
        x=ex_l[i][0]
        y=ex_l[i][1]
        c=abs((y-ly-m*(x-lx))*(m**2 + 1)**-0.5)
        d.append(c)
        epsilon=0.01
        if c < epsilon:
            f.append([x,y])
    if len(f) > 0 : 
        for i in f :
            d_1=(i[0]-lx)**2+(i[1]-ly)**2 #엡실론 내부의 좌표로부터 lx,ly 까지의 거리 
            h.append(d_1)
        k=h.index(min(h))
        return f[k][0], f[k][1],z
            
    else : #엡실론 내부에 해당하는 좌표가 없는경우에는 그냥 최소거리 출력 
        k=d.index(min(d))
        return ex_l[k][0], ex_l[k][1],z 


#---------------2차원 극판-----------------------------------------+
def pan(x, y, z, s, d, r):
    '''x,y 는 폐곡선 위의 찍은 점의 좌표, s는 찍은 점부터 극판 까지의 거리
    d 는 극판 두께, r은 극판 반지름 '''
    us=d + s # 점부터 극판 위쪽까지의 거리
    m=tan_yx(x,y,z)
    x_1, y_1, z =sam_dis(s,x,y,z)
    x_2, y_2, z =sam_dis(us,x,y,z)
    lx_1, ly_1, z, lx_2, ly_2, z = line_1(m,x_1,y_1,z,r,'r') #아래판의 좌우 끝 좌표 
    ux_1, uy_1, z, ux_2, uy_2, z = line_1(m,x_2,y_2,z,r,'r')  #위 판의 좌우 끝 좌표
    tx_1,ty_1,z=normal(-1/m,lx_1,ly_1,z)
    tx_2,ty_2,z=normal(-1/m,lx_2,ly_2,z)
    return lx_1, ly_1, lx_2, ly_2, ux_1, uy_1, ux_2, uy_2, tx_1, ty_1, tx_2,ty_2
#---------------3차원 평면 극판 -----------------------------------+
def pan_flat(x,y,z,s,d,r): #pan 3d 랑 비슷한 논리로 짰다. 문제는 이게 맞는질 모르겠네
    o_z, o =distance(x,y,z)
    i_z, j_z=curv_z(x,y,z,a)
    i=0
    list_mx_sur=[] # 폐곡선 위의 선들의 리스트 
    list_my_sur=[]
    list_dx_sur=[] # 아래쪽 극판 선들의 리스트
    list_dy_sur=[]
    list_ux_sur=[] # 위쪽 극판 선들의 리스트
    list_uy_sur=[]
    list_z=[]
    while o_z+i<= i_z:
        h_1=zi[o_z+i]-zi[o_z]
        r_u=abs((a**2-h_1**2))**0.5 #높이랑 반지름 이용해서 극판의 길이부분 구함 
        dx_1, dy_1, dx_2, dy_2, ux_1, uy_1, ux_2, uy_2, tx_1, ty_1, tx_2,ty_2 = pan(x,y,zi[o_z+i],s,d,r_u)
        list_dx_sur.append([dx_1,dx_2])
        list_dy_sur.append([dy_1,dy_2])
        list_ux_sur.append([ux_1,ux_2])
        list_uy_sur.append([uy_1,uy_2])
        list_z.append([zi[o_z+i]])
        
        mi, mj = curv(x,y,zi[o_z+i],r_u)
        mx_list, my_list= pick_sur(mi, mj)  
        list_mx_sur.append(mx_list)
        list_my_sur.append(my_list)
        
        i+=1
        
    list_dx_sur.reverse()
    list_dy_sur.reverse()
    list_ux_sur.reverse()
    list_uy_sur.reverse()
    list_z.reverse()
    list_mx_sur.reverse()
    list_my_sur.reverse()
    
    j=-1
    while o_z+j>= j_z :
        h_2=zi[o_z]-zi[o_z+j]
        r_d=abs((a**2-h_2**2))**0.5
        list_dx, list_dy, list_ux, list_uy = pan_2(x,y,zi[o_z+j],s,d,r_d)
        list_dx_sur.append(list_dx)
        list_dy_sur.append(list_dy)
        list_ux_sur.append(list_ux)
        list_uy_sur.append(list_uy)
        list_z.append([zi[o_z+j]])
        
        mi, mj = curv(x,y,zi[o_z+j],r_d)
        mx_list, my_list= pick_sur(mi, mj)  
        list_mx_sur.append(mx_list)
        list_my_sur.append(my_list)
        
        j-=1
              
    z_grid=np.array(list_z)
    
#---------------곡선 길이로, 해당 좌표 구하는 함수------------------------------------+
def curv(x,y,z,r):
    length_r=0
    length_l=0
    i_z, i=distance(x,y,z)
    j_z, j=distance(x,y,z)
    
    while length_r < r: # 길이 r 만큼 오른쪽으로 떨어진 곳의 좌표  
        if i < div-1:
            length_r+=((tan_yx(xi[i_z][i],yi[i_z][i],z)**2 + 1)**0.5)*abs((xi[i_z][i+1]-xi[i_z][i]))
            i+=1
        else : 
            length_r+=((tan_yx(xi[i_z][i],yi[i_z][i],z)**2 + 1)**0.5)*abs((xi[i_z][0]-xi[i_z][i]))
            i=0
    
    while length_l < r: # 길이 r 만큼 왼쪽으로 떨어진 곳의 좌표 
        if j > 0:
            length_l+=((tan_yx(xi[j_z][j],yi[j_z][j],z)**2 + 1)**0.5)*abs((xi[j_z][j]-xi[j_z][j-1]))
            j-=1
        else : 
            length_l+=((tan_yx(xi[j_z][j],yi[j_z][j],z)**2 + 1)**0.5)*abs((xi[j_z][j]-xi[j_z][div-1]))
            j=div-1
    
    return i,j


#-------------폐곡선 위에서 노말벡터만큼 떨어진 점 좌표 구하기-----------+
def sam_dis(a,x,y,z):
    m=tan_yx(x,y,z)
    x_1=x+(a*m*(m**2 + 1)**-0.5) #x,y 로 부터 기울기 -1/m 이고 a만큼 거리가 떨어진 점
    x_2=x-(a*m*(m**2 + 1)**-0.5)
    y_1=line(-1/m,x_1,x,y)
    y_2=line(-1/m,x_2,x,y)
    v_1=np.array([x_1,y_1,0]) 
    v_2=np.array([x_2,y_2,0])
    v=np.array([x,y,0])
    n_1=v_1-v #법선 벡터 1 
    n_2=v_2-v #법선 벡터 2
    i, j=distance(x,y,z)
    c=np.array([list_v[i][j][0],list_v[i][j][1],0])
    cross_1=np.cross(n_1,c)
    cross_2=np.cross(n_2,c)
    
    if cross_1[2]>0: #법선벡터1이랑 폐곡선 벡터랑 외적한 값이 양수면 그린다
        return  x_1, y_1, z
        
    elif cross_2[2]>0:
        return x_2, y_2, z

#---------------폐곡선에서 특정 거리만큼 떨어진 곳의 라인----------------+
def cur_line(x,y,z,s,r): #s 는 폐곡선에서 떨어진 거리 
    e, k=distance(x,y,z)
    i,j=curv(x,y,z,r)
    l_x=[]
    l_y=[]
    if k <= i :
        while k <= i :
            x_1, y_1, z_1 = sam_dis(s,xi[e][k],yi[e][k],z)
            k+=1
            l_x.append(x_1)
            l_y.append(y_1)
           
    else : 
        while k <= div-1 :
            x_1, y_1, z_1 = sam_dis(s,xi[e][k],yi[e][k],z) 
            k+=1
            l_x.append(x_1)
            l_y.append(y_1)
           
        while k <= i+div :
            x_1, y_1, z_1 = sam_dis(s,xi[e][k-div],yi[e][k-div],z) 
            k+=1
            l_x.append(x_1)
            l_y.append(y_1)
          
    r_x=[]
    r_y=[]
    
    e, k=distance(x,y,z)
    if k >= j :    
        while k >= j :
            x_2, y_2, z_2 = sam_dis(s,xi[e][k],yi[e][k],z)
            k-=1
            r_x.append(x_2)
            r_y.append(y_2)
          
    else:
        while k >= 0 :
            x_2, y_2, z_2 = sam_dis(s,xi[e][k],yi[e][k],z) 
            k-=1
            r_x.append(x_2)
            r_y.append(y_2)
           
        while k >= j-div :
            x_2, y_2, z_2 = sam_dis(s,xi[e][k+div],yi[e][k+div],z) 
            k-=1
            r_x.append(x_2)
            r_y.append(y_2)
        
    r_x.reverse()
    r_y.reverse()
    del r_x[-1]
    del r_y[-1]
    list_x=[]
    list_y=[]
    list_x=r_x+l_x
    list_y=r_y+l_y
    array_x=np.array(list_x)
    array_y=np.array(list_y)
    return array_x, array_y
    
#---------------2차원 곡면 극판-----------------------------------------+
def pan_2(x, y, z, s, d, r):
    '''x,y 는 폐곡선 위의 찍은 점의 좌표, s는 찍은 점부터 극판 까지의 거리
    d 는 극판 두께, r은 극판 반지름 '''
    us=d + s # 점부터 극판 위쪽까지의 거리
    list_dx, list_dy = cur_line(x,y,z,s,r)
    list_ux, list_uy = cur_line(x,y,z,us,r)
    return list_dx, list_dy, list_ux, list_uy

def pan_2m(x, y, z, r):
    ''' 3d 를 위해서 함수를 두개로 나눔. 판 위에서의 라인을 리스트로 뽑아내는것 '''
    i,j=curv(x,y,z,r)
    list_mx, list_my = cur_line(x,y,z,0,r)
    return list_mx, list_my

#----------z축 curv 함수------------------------------------------------+
def curv_z(x,y,z,a): # curv 함수와 비슷한 논리 로 짬 
    length_u=0
    length_d=0
    i_z, i=distance(x,y,z)
    j_z, j=distance(x,y,z)
    while length_u < a: # 길이 a 만큼 위로 떨어진 곳의 좌표  
        if i_z < div_z-1:
            length_u+=abs((zi[i_z+1]-zi[i_z]))
            i_z+=1
        
        else : #만약 z가 범위를 벗어나면 와일 구문이 무한 루프를 돔. 그래서 껴넌것 
            break
    
    while length_d < a: # 길이 a 만큼 아래로 떨어진 곳의 좌표 
        if j_z > 0:
            length_d+=abs((zi[j_z]-zi[j_z-1]))
            j_z-=1
            
        else: 
            break
    return i_z-1 ,j_z+1 #마지막 루프에서 한번씩 더 빠지고 더해지므로 그것을 보정한 것   

#------------3차원 곡면 극판--------------------------------------------+
def pan_3d(x,y,z,s,d,a): #s 는 찍은 점부터 극판까지의 거리, d는 극판의 두께, a는 극판 반지름 
    o_z, o =distance(x,y,z)
    i_z, j_z=curv_z(x,y,z,a)
    i=0
    list_mx_sur=[] # 폐곡선 위의 선들의 리스트 
    list_my_sur=[]
    list_dx_sur=[] # 아래쪽 극판 선들의 리스트
    list_dy_sur=[]
    list_ux_sur=[] # 위쪽 극판 선들의 리스트
    list_uy_sur=[]
    list_z=[]
    while o_z+i<= i_z:
        h_1=zi[o_z+i][0]-zi[o_z][0] # zi도 2차원 array 이므로 이렇게 뽑아야 한다
        r_u=abs((a**2-h_1**2))**0.5 #높이랑 반지름 이용해서 극판의 길이부분 구함 
        list_dx, list_dy, list_ux, list_uy = pan_2(x,y,zi[o_z+i][0],s,d,r_u)
        list_dx_sur.append(list_dx)
        list_dy_sur.append(list_dy)
        list_ux_sur.append(list_ux)
        list_uy_sur.append(list_uy)
        list_z.append(zi[o_z+i])
        
        mi, mj = curv(x,y,zi[o_z+i][0],r_u) #폐곡선 위에서 좌표 뽑기
        mx_list, my_list= pick_sur(mi, mj,o_z+i)  
        list_mx_sur.append(mx_list)
        list_my_sur.append(my_list)
        
        i+=1
        
    list_dx_sur.reverse()
    list_dy_sur.reverse()
    list_ux_sur.reverse()
    list_uy_sur.reverse()
    list_z.reverse()
    list_mx_sur.reverse()
    list_my_sur.reverse()
    
    j=-1
    while o_z+j>= j_z :
        h_2=zi[o_z][0]-zi[o_z+j][0]
        r_d=abs((a**2-h_2**2))**0.5
        list_dx, list_dy, list_ux, list_uy = pan_2(x,y,zi[o_z+j][0],s,d,r_d)
        list_dx_sur.append(list_dx)
        list_dy_sur.append(list_dy)
        list_ux_sur.append(list_ux)
        list_uy_sur.append(list_uy)
        list_z.append(zi[o_z+j])
        
        mi, mj = curv(x,y,zi[o_z+j][0],r_d)
        mx_list, my_list= pick_sur(mi, mj,o_z+j)  
        list_mx_sur.append(mx_list)
        list_my_sur.append(my_list)
        
        j-=1
              
    z_grid=np.array(list_z)
    #fit 함수를 통해 array의 개수를 최댓갓으로 통일하여 메쉬를 짤수있게 만든다
    dx_grid, dy_grid = fit(list_dx_sur, list_dy_sur) 
    ux_grid, uy_grid = fit(list_ux_sur, list_uy_sur)
    mx_grid, my_grid = fit(list_mx_sur, list_my_sur)
    
    #면 그리기
    ax.plot_surface(dx_grid, dy_grid, z_grid, color='r', alpha=1)
    ax.plot_surface(ux_grid, uy_grid, z_grid, color='r', alpha=1)
    ax.plot_surface(mx_grid, my_grid, z_grid, color='black', alpha=1)
    #옆면 그리기
    side_sur(dx_grid, dy_grid, z_grid, ux_grid, uy_grid, 'r')
    side_sur(dx_grid, dy_grid, z_grid, mx_grid, my_grid, 'black')
       
#------- 폐곡면위의 점들 뽑을때----------------------------------+
def pick_sur(i,j,k):
    if i+1>=j:
        x_list=xi[k][j:i+1]
        y_list=yi[k][j:i+1]
    else : #i+1 < j
        x_list1=xi[k][j:div]
        x_list2=xi[k][:i+1]
        x_list=x_list1+x_list2
        y_list1=yi[k][j:div]
        y_list2=yi[k][:i+1]
        y_list=y_list1+y_list2
    
    return x_list, y_list
    
    
#--------array의 갯수를 맞추기 위한 함수(최댓값으로 맞춘다)--------------------+
def fit(a,b): #a는 여기서 x array list를, b는 y array list를 의미한다
    how_many=[]
    re_a=[]
    re_b=[]
    for i in a:
        how_many.append(len(i))
    many= max(how_many)  # 가장 많은 어레이의 갯수를 찾고, 그 개수에 맞추어 나머지를 인터폴 레이트 함 
    for i in range(len(a)): 
        tck, ui =interpolate.splprep([a[i],b[i]], s=0, per=False)
        #per = 0 이면 폐곡선이 아니라 곡선 형태로 인터폴 레이션 한다. 
        xii,yii=interpolate.splev(np.linspace(0,1,many),tck)
        re_a.append(xii)
        re_b.append(yii)
    grid_a=np.array(re_a)
    grid_b=np.array(re_b)
    return grid_a, grid_b

#---------극판 옆면 만드는 함수------------------------------+
def side_sur(a_1, b_1, c, a_2, b_2, col):
    '''a_1, b_1 은 한쪽 극판의 x 어레이와 y 어레이를 의미한다. 
    어레이에서 각각 추출하여 세로운 어레이를 만드는 것
    근데 생각해보니깐 z 어레이도 2개 만들어야 할것 같은데...'''
   
    a1_list=pick_side(a_1)
    a2_list=pick_side(a_2)
    a_list=[]
    a_list.append(a1_list)
    a_list.append(a2_list)
    x_grid=np.array(a_list)
    
    b1_list=pick_side(b_1)
    b2_list=pick_side(b_2)
    b_list=[]
    b_list.append(b1_list)
    b_list.append(b2_list)
    y_grid=np.array(b_list)
        
    c_list=pick_side(c)
    c_list_d=[]
    c_list_d.append(c_list)
    c_list_d.append(c_list)
    z_grid=np.array(c_list_d)
    
    ax.plot_surface(x_grid, y_grid, z_grid, color=col, alpha=1)
    
def pick_side(a) : #각 어레이의 첫점과 끝점을 추출하여, 극판 부분의 둘레의 좌표를 빼내는 거다
    a_list1=[]
    a_list2=[]
    for i in range(len(a)):
        a_list1.append(a[i][0])
        a_list2.append(a[i][-1])
    a_list2.reverse()
    a_list=a_list1+a_list2
    a_list.append(a[0][0])
    return a_list
    

#-------------------클릭----------------------------------------------+
time = None #클릭 후 시간을 세기위한 장치 
def onclick(event):
    time_interval = 0.25 #0.25초 이내에 더블클릭해야 인식함 
    global time
    if event.button==3: #우클릭시
        p=ax.format_coord(event.xdata,event.ydata) 
        #matplotlib 내장함수. 클릭 위치의 좌표 string으로 추출 
        kx,ky,kz=cor(p)
        print(p)
        
        if time is None:
            time = threading.Timer(time_interval, on_singleclick, [event,kx,ky,kz,d_1,d_2,a]) #arg를 튜플형태로 넣어서 싱글클릭에 넣는듯? 
            time.start()
            
        if event.dblclick:
            time.cancel()
            ax.scatter(kx, ky, kz, color='green')
            on_dblclick(event,kx,ky,kz,s,d_1,d_2,a)
            

#--------------------MAIN--------------------------------------------+    
##----------- 따블 클릭할 때 ---------------------------------------+
click_num=0
def on_dblclick(event,x,y,z,s,d_1,d_2,a):
    global time
    print("You double-clicked", event.button, event.xdata, event.ydata)
    time = None
    flat_plate(x,y,z,s,d_1,d_2,a)
    global click_num 
    click_num=click_num+1
    print(click_num)
    return click_num
    
##----------- 싱글 클릭할 때 ---------------------------------------+

def on_singleclick(event,x,y,z,d_1,d_2,a):
    global time
    print("You single-clicked", event.button, event.xdata, event.ydata)
    time = None
    pass

cid = fig.canvas.mpl_connect('button_press_event', onclick)

#-------------------------------------------------------------------------
