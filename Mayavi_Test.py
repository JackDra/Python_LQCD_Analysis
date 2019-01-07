#!/usr/bin/env python

%gui qt

from mayavi import mlab
import mayavi as ma
import numpy as np
# mlab.test_contour3d()

n_turn = 100
compress = 10
t = np.linspace(0,n_turn*2*np.pi,n_turn*10)
x,y,z = np.sin(t),np.cos(t),t/compress
mlab.plot3d(x,y,z)
mlab.plot3d([x[0],0],[y[0],0],[z[0],z[0]])
mlab.plot3d([0,0],[0,0],[z[0],-1])
mlab.points3d([0],[0],[-1],mode='cube')

r = 1
theta,phi = np.mgrid[0:2*np.pi:20j,0:np.pi:20j]
x,y,z = r*np.cos(theta)*np.sin(phi),r*np.sin(theta)*np.sin(phi),r*np.cos(phi)
mlab.mesh(x,y,z,representation='wireframe')


points = np.array([[0,0,0],[0,1,0],[1,1,0],[1,0,0],[0.5,0.5,1]])
t = [[0,1,2],[0,3,2],[0,1,4],[1,2,4],[2,3,4],[3,0,4]]
x,y,z = points.T
mlab.triangular_mesh(x,y,z,t)

vals = np.random.random((2<<6,2<<6))
mlab.imshow(vals)


x,y,z = np.mgrid[-5:5:64j,-5:5:64j,-5:5:64j]
mlab.contour3d(0.5*x**2 + y**2 + 2*z**2)
mlab.volume_slice(x,y,z,0.5*x*x + y*y + 2*z*z)

mlab.test_quiver3d()


x,y,z = np.mgrid[-2:3:50j,-2:3:50j,-2:3:50j]
r = np.sqrt(x**2 + y**2 + z**2)
u = y*np.sin(r)/(r+0.001)
v = -x*np.sin(r)/(r+0.001)
w = np.ones_like(z)*0.1
mlab.flow(x,y,z,u,v,w,seedtype='plane')


def lorenz(x,y,z,s=10,r=28,b=8/3):
    u = s*(y-x)
    v = r*x-y-x*z
    w = x*y-b*z
    return u,v,w
x,y,z = np.mgrid[-50:50:20j,-50:50:20j,-10:60:20j]
u,v,w = lorenz(x,y,z)
mlab.quiver3d(x,y,z,u,v,w,scale_factor=0.01,mask_points=5)
mlab.flow(x,y,z,u,v,w)

from scipy.integrate import odeint
def lorenz_ode(state,t):
    return np.array(lorenz(*state))

t = np.linspace(0,50,2000)
x,y,z = odeint(lorenz_ode,(10,50,50),t).T
mlab.colorbar()
mlab.plot3d(x,y,z,t,tube_radius=None)
mlab.axes()
mlab.plot3d()
mlab.show_pipeline()


import time
@mlab.animate
def anim():
    x,y = np.mgrid[0:3:1,0:3:1]
    s = mlab.surf(x,y,x*0.01,representation='wireframe')
    fig = mlab.gcf()
    while True:
        for i in range(5):
            s1 = slice(0,3,1/(i+2))
            x,y = np.mgrid[s1,s1]
            sc = x*x*0.05*(i+1)
            s.mlab_source.reset(x=x,y=y,scalars = sc)
            fig.scene.reset_zoom()
            mlab.process_ui_events()
            time.sleep(0.02)
            yield
anim()

import time
@mlab.animate
def sin_animate():
    x,y = np.mgrid[-3:3:100j,-3:3:100j]
    z = np.sin(x**2 + y**2)
    s = mlab.surf(x,y,z)
    t = 0
    while True:
        t += 0.01
        s.mlab_source.set(scalars = np.sin(t+(x**2 + y**2)))
        yield
sin_animate()





def lorenz(x,y,z,s=10,r=28,b=8/3):
    u = s*(y-x)
    v = r*x-y-x*z
    w = x*y-b*z
    return u,v,w
x,y,z = np.mgrid[-50:50:20j,-50:50:20j,-10:60:20j]
u,v,w = lorenz(x,y,z)
mlab.quiver3d(x,y,z,u,v,w,scale_factor=0.01,mask_points=5)
mlab.flow(x,y,z,u,v,w)

a = np.random.random((4,4))
mlab.surf(a)
fig = mlab.get_engine()
surface = fig.scenes[0].children[0].children[0].children[0].children[0].children[0]
surface.contour.contours = [0.5195050456410771]
surface.actor.mapper.scalar_range = np.array([0.1108337 , 0.92817639])
surface.actor.mapper.progress = 1.0
# surface.actor.mapper.input_connection = <tvtk.tvtk_classes.algorithm_output.AlgorithmOutput object at 0x7fecec3f27d8>
# surface.actor.mapper.input_connection = <tvtk.tvtk_classes.algorithm_output.AlgorithmOutput object at 0x7fecec3f27d8>
surface.enable_contours = True

src = mlab.pipeline.array2d_source(a)
warp = mlab.pipeline.warp_scalar(src)
poly = mlab.pipeline.poly_data_normals(warp)
mlab.pipeline.surface(poly)

mlab.test_flow()
from tvtk.api import tvtk
from scipy import special
x,y = np.mgrid[-10:10:20j,-10:10:20j]
r = np.sqrt(x**2 + y**2)
z = 5*special.j0(r)
spoints = tvtk.StructuredPoints(origin=(-12.5,-12.5,0),
                                spacing=(0.5,0.5,1),
                                dimensions=(20,20,1))
spoints.point_data.scalars = z.T.ravel()
spoints.point_data.scalars.name = 'scalar'
src = mlab.pipeline.add_dataset(spoints)
warp = mlab.pipeline.warp_scalar(src)
poly = mlab.pipeline.poly_data_normals(warp)
surf = mlab.pipeline.surface(poly)

x,y,z = np.mgrid[-5:5:128j,-5:5:128j,-5:5:128j]
scalars = np.sin(x*y*z)/(x*y*z)
spoints = tvtk.StructuredPoints(origin=(-5,-5,-5),
                                spacing=(10/127,10/127,10/127),
                                dimensions=(128,128,128))
s = scalars.T
spoints.point_data.scalars = s.ravel()
spoints.point_data.scalars.name = 'scalars'
src = mlab.pipeline.add_dataset(spoints)
cut = mlab.pipeline.scalar_cut_plane(src)
contour = mlab.pipeline.iso_surface(src)

r,th,z = np.mgrid[1:10:25j,0:2*np.pi:51j,0:5:25j]
x,y = np.cos(th)*r,np.sin(th)*r
scalar = x**2 + y**2 + z**2

pts = np.empty(z.shape + (3,))
pts[...,0] = x
pts[...,1] = y
pts[...,2] = z

pts = pts.transpose(2,1,0,3).copy()
pts = pts.reshape(pts.size//3,3)

sgrid = tvtk.StructuredGrid(dimensions=x.shape)
sgrid.points = pts
sgrid.point_data.scalars = np.ravel(scalar.T.copy())
sgrid.point_data.scalars.name = 'scalar'
src = mlab.pipeline.add_dataset(sgrid)
plane = mlab.pipeline.grid_plane(src)
plane.grid_plane.axis = 'x'
plane_y = mlab.pipeline.grid_plane(src)
plane_y.grid_plane.axis = 'y'
plane_x = mlab.pipeline.grid_plane(src)
plane_x.grid_plane.axis = 'z'
c_plane = mlab.pipeline.contour_grid_plane(src)
c_plane.enable_contours = False
iso = mlab.pipeline.iso_surface(src)


points = np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1]])
triangles = np.array([[0,1,3],[0,3,2],[1,2,3],[0,2,1]])
mesh = tvtk.PolyData()
mesh.points = points
mesh.polys = triangles

temp = np.array([10,20,30,40])
mesh.point_data.scalars = temp
mesh.point_data.scalars.name = 'temp'


velo = np.array([[0,0,0],[1,0,0],[0,1,0],[0,0,1]])
mesh.point_data.vectors = velo
mesh.point_data.vectors.name = 'velo'

src = mlab.pipeline.add_dataset(mesh)
surf = mlab.pipeline.surface(src)
vec = mlab.pipeline.vectors(src)
