#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d.axes3d as p3
import matplotlib.animation as animation
from matplotlib.backend_bases import NavigationToolbar2, Event

home = NavigationToolbar2.home

def new_home(self, *args, **kwargs):
    s = 'home_event'
    event = Event(s, self)
    event.foo = 100
    self.canvas.callbacks.process(s, event)
    home(self, *args, **kwargs)

NavigationToolbar2.home = new_home

prev_next = NavigationToolbar2.forward

def new_next(self, *args, **kwargs):
    s = 'next_event'
    event = Event(s, self)
    event.foo = 100
    self.canvas.callbacks.process(s, event)
    prev_next(self, *args, **kwargs)

NavigationToolbar2.forward = new_next

prev_prev = NavigationToolbar2.back

def new_prev(self, *args, **kwargs):
    s = 'prev_event'
    event = Event(s, self)
    event.foo = 100
    self.canvas.callbacks.process(s, event)
    prev_prev(self, *args, **kwargs)

NavigationToolbar2.back = new_prev

file_name = '/home/jackdra/LQCD/Scripts/Python_Analysis/Configs/test.txt'

x1_list,x2_list,n_list = [[]],[[]],[[]]
with open(file_name,'r') as f:
    for line in f:
        strp_line = line.strip()
        if len(strp_line) == 0: continue
        if strp_line[0]=='#':
            x1_list.append([])
            x2_list.append([])
            n_list.append([])
        else:
            split_line = list(map(float,strp_line.split(' ')))
            x1_list[-1].append(split_line[:3])
            x2_list[-1].append(split_line[3:6])
            n_list[-1].append(int(split_line[6]))
x1_max = np.max(list(map(len,x1_list)))
for ic,ix1 in enumerate(x1_list):
    this_len = len(ix1)
    if this_len < x1_max:
        x1_list[ic] += [(0,0,0) for _ in range(this_len,x1_max)]
        x2_list[ic] += [(0,0,0) for _ in range(this_len,x1_max)]
        n_list[ic] += [1 for _ in range(this_len,x1_max)]
x1_list,x2_list,x3_list = np.array(x1_list),np.array(x2_list),np.array(n_list)
color_list = ['None','blue','red']


# Attaching 3D axis to the figure
fig = plt.figure()
ax = p3.Axes3D(fig)

# Fifty lines of random 3-D lines
# data = [Gen_RandLine(25, 3) for index in range(50)]
# data
# Creating fifty line objects.
# NOTE: Can't pass empty arrays into 3d version of plot()
# x1_list = np.swapaxis(x1_list,0,1)
# x2_list = np.swapaxis(x2_list,0,1)
this_data = [x1_list,x2_list,n_list]
lines = []
for ic,(ix1,ix2,icn) in enumerate(zip(x1_list[0],x2_list[0],n_list[0])):
    if ic == 0:
        lines.append(ax.plot([ix1[0],ix1[0]+ix2[0]],
                             [ix1[1],ix1[1]+ix2[1]],
                             [ix1[2],ix1[2]+ix2[2]],color = color_list[icn],label=str(ic+1))[0])
    else:
        lines.append(ax.plot([ix1[0],ix1[0]+ix2[0]],
                             [ix1[1],ix1[1]+ix2[1]],
                             [ix1[2],ix1[2]+ix2[2]],color = color_list[icn])[0])
# Setting the axes properties
ax.set_xlim3d([-2, 2])
ax.set_xlabel('X')

ax.set_ylim3d([-2, 2])
ax.set_ylabel('Y')

ax.set_zlim3d([-2, 2])
ax.set_zlabel('Z')

ax.set_title('3D Test')
ax.legend()
# Creating the Animation object

def update_lines(num, dataLines,lines):
    for icl,iline in enumerate(lines):
        ix1 = x1_list[num][icl]
        ix2 = x2_list[num][icl]
        icn = n_list[num][icl]
        if icl == 0:
            iline.set_label(str(num+1))
        iline.set_color(color_list[icn])
        iline.set_data( [[ix1[0],ix1[0]+ix2[0]],
                        [ix1[1],ix1[1]+ix2[1]]])
        iline.set_3d_properties([ix1[2],ix1[2]+ix2[2]])
        ax.legend()
    return lines


pause = False
t=0
def simData():
    t_max = len(x1_list)-1
    global t
    t = 0
    while t < t_max:
        if not pause:
            t = t +1
        yield t

    # print(evt.foo)

def onClick(event):
    global pause
    pause ^= True
# fig.canvas.mpl_connect('button_press_event', onClick)


def handle_home(evt):
    global pause
    pause ^= True
    if pause:
        print('\r paused',end='      ')
    else:
        print('\r un-paused',end='     ')

def handle_next(evt):
    global t
    print('\r Moving next',end='      ')
    t = (t + 1)%(len(x1_list)-1)

def handle_prev(evt):
    global t
    print('\r Moving previous',end='      ')
    t = (t - 1)%(len(x1_list)-1)

fig.canvas.mpl_connect('home_event', handle_home)
fig.canvas.mpl_connect('prev_event', handle_prev)
fig.canvas.mpl_connect('next_event', handle_next)
# home = NavigationToolbar2.home
line_ani = animation.FuncAnimation(fig, update_lines, simData,fargs=(this_data,lines),
                                   interval=1000, blit=False,repeat=True)

plt.show()
