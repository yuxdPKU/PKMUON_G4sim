import matplotlib.pyplot as plt

'''
# 定义两个点的坐标

#smear
x1_s, y1_s, z1_s = 1.285598, 1.401191, 13.340000
x2_s, y2_s, z2_s = 0.945970, 1.537330, 13.340000
#l1: (1.285598 1.401191 13.340000)
#l2: (0.945970 1.537330 13.340000)

#no smear
x1, y1, z1 = 1.331506, 1.359462, 13.340000
x2, y2, z2 = 1.327560, 1.374350, 13.340000
#l1: (1.331506 1.359462 13.340000)
#l2: (1.327560 1.374350 13.340000)

# 创建一个新的图形
#fig, ax = plt.subplots()
fig = plt.figure()
ax = plt.axes(projection='3d')

# 画出两个点的连线
#ax.plot([x1, x2], [y1, y2],color='red')
#ax.plot([x1_s, x2_s], [y1_s, y2_s],color='blue')
ax.scatter3D(x, y, z)

# 显示图形
plt.show()
'''


# 定义两个点的坐标
#no smear
x = [0,1.327560]
y = [0,1.374350]
z = [0,13.340000]

x_s = [0,1.331506]
y_s = [0,1.359462]
z_s = [0,13.340000]

'''
x = [1.331506,1.327560]
y = [1.359462,1.374350]
z = [13.340000,13.340000]

x_s = [1.285598,0.945970]
y_s = [1.401191,1.537330]
z_s = [13.340000,13.340000]
'''

x1 = [0,1.285598]
y1 = [0,1.401191]
z1 = [0,13.340000]

x1_s = [0,0.945970]
y1_s = [0,1.537330]
z1_s = [0,13.340000]


# 创建一个新的三维图形
fig = plt.figure()
ax = plt.axes(projection='3d')

# 画出两个点的连线
ax.plot(x, y, z, color='red')
ax.plot(x_s, y_s, z_s, color='blue')

ax.plot(x1, y1, z1, color='yellow')
ax.plot(x1_s, y1_s, z1_s, color='green')



# 显示图形
plt.show()