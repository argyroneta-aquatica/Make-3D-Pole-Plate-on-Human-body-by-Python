# -*- coding: utf-8 -*-
"""
Created on Wed Jul 29 11:17:20 2020

@author: 용하
"""

#------------모듈-------------------------------------------+
import numpy as np
from scipy import interpolate
from matplotlib import pyplot as plt

##-----------패러미터--------------------------------------------------+
s=0.1 #폐곡선 표면부터 극판까지의 거리 
d=0.1 #극판 두께
r=0.5 #극판의 반지름 
div = 400 # 폐곡선 나눈 수
#-------------폐곡선-----------------------------------------+
x = np.array([9.5, 13, 14, 12.5, 10])
y = np.array([12, 10, 14.5, 14, 12])

# x 와 x[0] array 를 합치는 것(처음이랑 끝을 이어서 폐곡선 만든 것) 
x = np.r_[x, x[0]]
y = np.r_[y, y[0]]

# x,y 사이사이 값에 추가 적인 값을 넣어, 스무스한 곡선을 만드는 것 
tck, u = interpolate.splprep([x, y], s=0, per=True)

# 스무스한 곡선에서 1000의 점 빼기 
xi, yi = interpolate.splev(np.linspace(0, 1, div), tck)
#------------폐곡선 라인 벡터 만들기 (반시계 방향)--------------------+
vector_1=[]
for i in range(len(xi)):
    vector_1.append([xi[i],yi[i]])
l=np.array(vector_1)
cur_v=[]
for i in range(len(l)):
    if i<div-1:
        cur_v.append(l[i+1]-l[i])
    
    else:
        cur_v.append(l[0]-l[i])
#---------폐곡선 그래프 그리기-----------------------------------+
fig, ax = plt.subplots(1, 1)
plt.xlim(9,15)
plt.ylim(9,15)
ax.plot(xi, yi, '-b')


#------------좌표추출-------------------------------------------+
def cor(p):
    s=p.split('y') # 나눌게 애매해서 중간에 y로 나눔   
    b=[]
    for i in range(len(s)):
         b.append(s[i].strip())
    c=[]
    c.append(b[0][2:])
    c.append(b[1][1:])
    cx=float(c[0])
    cy=float(c[1])
    k=distance(cx, cy)
    kx=xi[k]
    ky=yi[k]
    return kx, ky


#-------------거리가까운 메시의 좌표특정을 위한 거 ---------------------+
def distance(x,y):
    d=[]
    for i in range(len(xi)):
        c=(x-xi[i])**2+(y-yi[i])**2
        d.append(c)
    e=d.index(min(d))
    return e 
'''리턴값을 e로해서, 좌표자체를 특정하는게 아니라 좌표가 리스트에서 존재하는 
위치를 특정하는거임 '''

#-----------접선그리기--------------------------------------------+
def tan(x,y):
    k=distance(x,y)
    if 0 < k < div-1:
        m1=(yi[k+1]-yi[k])/(xi[k+1]-xi[k])
        m2=(yi[k]-yi[k-1])/(xi[k]-xi[k-1])
        m=(m1+m2)/2
        return m
    elif k == 0 :
        m1=(yi[1]-yi[0])/(xi[1]-xi[0])
        m2=(yi[-1]-yi[-2])/(xi[-1]-xi[-2])
        m=(2*m1 + m2)/3
        return m
    elif k == div-1 : #마지막 점의 경우 k+1 번째가 없기에 에러가 난다. 때문에 써준거 
        m1=(yi[1]-yi[0])/(xi[1]-xi[0])
        m2=(yi[-1]-yi[-2])/(xi[-1]-xi[-2])
        m=(m1 + 2*m2)/3
        return m 

def line_1(m,x,y,b,color):
    ''' m 은 기울기, b 는 중심점부터 끝 까지의 길이
    좌우 끝점의 좌표를 리턴한다'''
    x_1=x-(b)*(1+m**2)**-0.5
    x_2=x+(b)*(1+m**2)**-0.5
    y_1= line(m,x_1,x,y)
    y_2= line(m,x_2,x,y)
    plt.plot([x_1,x_2],[y_1,y_2], color)
    return x_1, y_1, x_2, y_2 
#------------라인 함수--------------------------------------------+ 
def line(m,x, xt, yt):
    return m*(x-xt)+yt   

#-----------노말벡터 그리기-----------------------------------+
def normal(m,lx,ly):
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
        return f[k][0], f[k][1]
            
    else : #엡실론 내부에 해당하는 좌표가 없는경우에는 그냥 최소거리 출력 
        k=d.index(min(d))
        return ex_l[k][0], ex_l[k][1] 


#---------------2차원 극판-----------------------------------------+
def pan(x, y, s, d, r):
    '''x,y 는 폐곡선 위의 찍은 점의 좌표, s는 찍은 점부터 극판 까지의 거리
    d 는 극판 두께, r은 극판 반지름 '''
    us=d + s # 점부터 극판 위쪽까지의 거리
    m=tan(x,y)
    x_1, y_1 =sam_dis(s,x,y)
    x_2, y_2 =sam_dis(us,x,y)
    lx_1, ly_1, lx_2, ly_2 = line_1(m,x_1,y_1,r,'r') #아래판의 좌우 끝 좌표 
    ux_1, uy_1, ux_2, uy_2 = line_1(m,x_2,y_2,r,'r')  #위 판의 좌우 끝 좌표
    plt.plot([lx_1,ux_1],[ly_1,uy_1],'r')
    plt.plot([lx_2,ux_2],[ly_2,uy_2],'r') #판끼리 연결함 
    tx_1,ty_1=normal(-1/m,lx_1,ly_1)
    tx_2,ty_2=normal(-1/m,lx_2,ly_2)
    plt.plot([lx_1,tx_1],[ly_1,ty_1],'black')
    plt.plot([lx_2,tx_2],[ly_2,ty_2],'black')
    
#---------------곡선 길이로, 해당 좌표 구하는 함수------------------------------------+
def curv(x,y,r):
    length_r=0
    length_l=0
    i=distance(x,y)
    j=distance(x,y)
    
    while length_r < r: # 길이 r 만큼 오른쪽으로 떨어진 곳의 좌표  
        if i < div-1:
            length_r+=((tan(xi[i],yi[i])**2 + 1)**0.5)*abs((xi[i+1]-xi[i]))
            i+=1
            
        else : 
            length_r+=((tan(xi[i],yi[i])**2 + 1)**0.5)*abs((xi[0]-xi[i]))
            i=0
        
    
    while length_l < r: # 길이 r 만큼 왼쪽으로 떨어진 곳의 좌표 
        if j > 0:
            length_l+=((tan(xi[j],yi[j])**2 + 1)**0.5)*abs((xi[j]-xi[j-1]))
            j-=1
            
        else : 
            length_l+=((tan(xi[j],yi[j])**2 + 1)**0.5)*abs((xi[j]-xi[div-1]))
            j=div-1
        
    return i,j


#-------------폐곡선 위에서 노말벡터만큼 떨어진 점 좌표 구하기-----------+
def sam_dis(a,x,y):
    m=tan(x,y)
    x_1=x+(a*m*(m**2 + 1)**-0.5) #x,y 로 부터 기울기 -1/m 이고 a만큼 거리가 떨어진 점
    x_2=x-(a*m*(m**2 + 1)**-0.5)
    y_1=line(-1/m,x_1,x,y)
    y_2=line(-1/m,x_2,x,y)
    v_1=np.array([x_1,y_1,0]) 
    v_2=np.array([x_2,y_2,0])
    v=np.array([x,y,0])
    n_1=v_1-v #법선 벡터 1 
    n_2=v_2-v #법선 벡터 2
    e=distance(x,y)
    c=np.array([cur_v[e][0],cur_v[e][1],0])
    cross_1=np.cross(n_1,c)
    cross_2=np.cross(n_2,c)
    
    if cross_1[2]>0: #법선벡터1이랑 폐곡선 벡터랑 외적한 값이 양수면 그린다
        return  x_1, y_1
        
    elif cross_2[2]>0:
        return x_2, y_2

#---------------폐곡선에서 특정 거리만큼 떨어진 곳의 라인----------------+
def cur_line(x,y,s,r):
    k=distance(x,y)
    i,j=curv(x,y,r)
    l_x=[]
    l_y=[]
    if k <= i :
        while k <= i :
            x_1, y_1 = sam_dis(s,xi[k],yi[k])
            k+=1
            l_x.append(x_1)
            l_y.append(y_1)
           
    else : 
        while k <= div-1 :
            x_1, y_1 = sam_dis(s,xi[k],yi[k]) 
            k+=1
            l_x.append(x_1)
            l_y.append(y_1)
           
        while k <= i+div :
            x_1, y_1 = sam_dis(s,xi[k-div],yi[k-div]) 
            k+=1
            l_x.append(x_1)
            l_y.append(y_1)
          
    r_x=[]
    r_y=[]
    
    k=distance(x,y)
    if k >= j :    
        while k >= j :
            x_2, y_2 = sam_dis(s,xi[k],yi[k])
            k-=1
            r_x.append(x_2)
            r_y.append(y_2)
          
    else:
        while k >= 0 :
            x_2, y_2 = sam_dis(s,xi[k],yi[k]) 
            k-=1
            r_x.append(x_2)
            r_y.append(y_2)
           
        while k >= j-div :
            x_2, y_2 = sam_dis(s,xi[k+div],yi[k+div]) 
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
    ax.plot(list_x, list_y, 'r')
    
#---------------2차원 곡면 극판-----------------------------------------+
def pan_2(x, y, s, d, r):
    '''x,y 는 폐곡선 위의 찍은 점의 좌표, s는 찍은 점부터 극판 까지의 거리
    d 는 극판 두께, r은 극판 반지름 '''
    us=d + s # 점부터 극판 위쪽까지의 거리
    cur_line(x,y,s,r)
    cur_line(x,y,us,r)
    i,j=curv(x,y,r)
    lx_1, ly_1= sam_dis(s,xi[i],yi[i])
    lx_2, ly_2= sam_dis(s,xi[j],yi[j])
    ux_1, uy_1= sam_dis(us,xi[i],yi[i])
    ux_2, uy_2= sam_dis(us,xi[j],yi[j])
    plt.plot([lx_1,ux_1],[ly_1,uy_1],'r')
    plt.plot([lx_2,ux_2],[ly_2,uy_2],'r') #판끼리 연결함 
    plt.plot([lx_1,xi[i]],[ly_1,yi[i]],'black')
    plt.plot([lx_2,xi[j]],[ly_2,yi[j]],'black') #아래판이랑 폐곡선 연결 
    
    
#--------------------MAIN--------------------------------------------+    
##----------- 따블 클릭할 때 ---------------------------------------+
def onclick(event):
    if event.dblclick:
        p=ax.format_coord(event.xdata,event.ydata) 
        #matplotlib 내장함수. 클릭 위치의 좌표 string으로 추출 
        print(p)         
        kx,ky=cor(p)
        ax.scatter(kx, ky, color='green')
        if event.button==3: #우클릭시 
            pan_2(kx,ky,s,d,r)
        elif event.button==1: #좌클릭시 
            pan(kx,ky,s,d,r)
        plt.show()
        
cid = fig.canvas.mpl_connect('button_press_event', onclick)
#-------------------------------------------------------------------------