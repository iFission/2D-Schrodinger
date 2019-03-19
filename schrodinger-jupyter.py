#%% [markdown]
# # Schrödinger Assignment
#%%
import numpy as np
import scipy.constants as c


def deg_to_rad(deg):
    return round(deg / 180 * np.pi, 5)


def rad_to_deg(rad):
    return round(rad / np.pi * 180, 5)


#%%
def spherical_to_cartesian(r, theta, phi):
    x = r * np.sin(theta) * np.cos(phi)
    y = r * np.sin(theta) * np.sin(phi)
    z = r * np.cos(theta)
    x, y, z = round(x, 5), round(y, 5), round(z, 5)
    return round(x, 5), round(y, 5), round(z, 5)


def cartesian_to_spherical(x, y, z):
    r = np.sqrt(np.power(x, 2) + np.power(y, 2) + np.power(z, 2))
    theta = np.arccos(z / r)
    phi = np.arctan2(y, x)
    r, theta, phi = round(r, 5), round(theta, 5), round(phi, 5)
    return r, theta, phi


#%%
import numpy as np
from scipy.constants import pi
import scipy.constants as c

j = 1j


def angular_wave_func(m, l, theta, phi):
    if l == 0 and m == 0:
        Y = np.sqrt(1 / (4 * pi))
    elif l == 1:
        if m == 1:
            Y = -np.sqrt(3 / (8 * c.pi)) * np.sin(theta) * np.exp(phi * j)
        elif m == 0:
            Y = np.sqrt(3 / (4 * pi)) * np.cos(theta)
        elif m == -1:
            Y = np.sqrt(3 / (8 * c.pi)) * np.sin(theta) * np.exp(phi * -j)
    elif l == 2:
        if m == 0:
            Y = np.sqrt(5 / (16 * pi)) * (3 * np.power(np.cos(theta), 2) - 1)
        elif m == 1:
            Y = -np.sqrt(15 /
                         (8 * pi)) * np.cos(theta) * np.sin(theta) * np.exp(
                             phi * j)
        elif m == -1:
            Y = np.sqrt(15 /
                        (8 * pi)) * np.cos(theta) * np.sin(theta) * np.exp(
                            phi * -j)
    return np.round(Y, 5)


#%%
import numpy as np
import scipy.constants as c
a = c.physical_constants['Bohr radius'][0]


def radial_wave_func(n, l, r):
    if n == 1 and l == 0:
        R = 2 * np.exp(-r / a)
    elif n == 2:
        if l == 0:
            R = 1 / np.sqrt(2) * (1 - r / (2 * a)) * np.exp(-r / (2 * a))
        elif l == 1:
            R = 1 / np.sqrt(24) * (r / a) * np.exp(-r / (2 * a))
    elif n == 3:
        if l == 0:
            R = 2 / (81 * np.sqrt(3)) * (27 - 18 * (r / a) + 2 * np.power(
                (r / a), 2)) * np.exp(-r / (3 * a))
        elif l == 1:
            R = 8 / (27 * np.sqrt(6)) * (1 - r / (6 * a)) * (r / a) * np.exp(
                -r / (3 * a))
        elif l == 2:
            R = 4 / (81 * np.sqrt(30)) * np.power(
                (r / a), 2) * np.exp(-r / (3 * a))
    elif n == 4:
        if l == 0:
            R = 1 / 4 * (1 - 3 / 4 * (r / a) + 1 / 8 + np.power(
                (r / a), 2) - 1 / 192 * np.power(
                    (r / a), 3)) * np.exp(-r / (4 * a))
    return np.round(R, 5)


#%%
def frange(start, end, points):
    ls = []
    step = (end - start) / (points - 1)

    # used points as counter as opposed to while loop (while start <= end), to avoid floating point error
    for i in range(points):
        ls.append(float(start))
        start += step
    return ls


def mgrid2d(xstart, xend, xpoints, ystart, yend, ypoints):
    return [[[x] * ypoints for x in frange(xstart, xend, xpoints)],
            [frange(ystart, yend, ypoints)] * xpoints]


#%%
def frange(start, end, points):
    ls = []
    step = (end - start) / (points - 1)

    # used points as counter as opposed to while loop (while start <= end), to avoid floating point error
    for i in range(points):
        ls.append(float(start))
        start += step
    return ls


def mgrid3d(xstart, xend, xpoints, ystart, yend, ypoints, zstart, zend,
            zpoints):
    x_list = frange(xstart, xend, xpoints)
    y_list = frange(ystart, yend, ypoints)
    z_list = frange(zstart, zend, zpoints)
    x = [[[x] * zpoints] * ypoints for x in x_list]
    y = [[[y] * zpoints for y in y_list]] * xpoints
    z = [[z_list] * ypoints] * xpoints
    return x, y, z


#%%

import numpy as np


def real_angular_wave_func(m, l, theta, phi):
    if m < 0:
        return np.real(
            1j / np.sqrt(2) * (angular_wave_func(m, l, theta, phi) - (
                (-1)**m) * angular_wave_func(-m, l, theta, phi)))
    elif m == 0:
        return np.real(angular_wave_func(m, l, theta, phi))
    else:
        return np.real(
            1 / np.sqrt(2) * (angular_wave_func(-m, l, theta, phi) + (
                (-1)**m) * angular_wave_func(m, l, theta, phi)))


def hydrogen_wave_func(n, l, m, roa, Nx, Ny, Nz):
    xx, yy, zz = tuple(
        np.round(
            np.array(mgrid3d(-roa, roa, Nx, -roa, roa, Ny, -roa, roa, Nz)), 5))

    # '''vectorize cartesian to spherical function'''
    # cartesian_to_spherical_vectorized = np.vectorize(cartesian_to_spherical)
    # r, theta, phi = cartesian_to_spherical_vectorized(xx, yy, zz)
    r, theta, phi = np.vectorize(cartesian_to_spherical)(xx, yy, zz)

    # return r, theta, phi
    ang = np.vectorize(real_angular_wave_func)(m, l, theta, phi)
    rad = np.vectorize(radial_wave_func)(
        n, l, r * a)  # radius multiply by a (Bohr radius)
    real_wave_func = np.absolute(ang * rad)**2
    return tuple(
        np.round([xx, yy, zz, real_wave_func],
                 5))  # use np's round function, convert nparray to tuple


#%% [markdown]
# ## Week 9:
#%% [markdown]
# Use Mayavi to plot the real orbitals of your assigned hydrogen function. The real orbitals will be a linear combination of your complex wave functions.

#%%
# Code to save the data to a file so that
# you don't have to keep on computing it:

print('Test ')
x, y, z, mag = hydrogen_wave_func(4, 1, -1, 40, 100, 100, 100)
print('x, y, z:')
print(x, y, z)
print('mag:')
print(mag)
print(x, y, z, mag)
x.dump('x_test.dat')
y.dump('y_test.dat')
z.dump('z_test.dat')
mag.dump('den_test.dat')

#%%
from mayavi import mlab

mu, sigma = 0, 0.1
x = np.load('x_test.dat')
y = np.load('y_test.dat')
z = np.load('z_test.dat')

density = np.load('den_test.dat')
figure = mlab.figure('DensityPlot')
pts = mlab.contour3d(density, contours=40, opacity=0.4)
mlab.axes()
mlab.show()

#%%

import numpy as np

from traits.api import HasTraits, Instance, Array, on_trait_change
from traitsui.api import View, Item, HGroup, Group

from tvtk.api import tvtk
from tvtk.pyface.scene import Scene

from mayavi import mlab
from mayavi.core.api import PipelineBase, Source
from mayavi.core.ui.api import SceneEditor, MayaviScene, MlabSceneModel
# Create some data
data = np.load('den_test.dat')


# The object implementing the dialog
class VolumeSlicer(HasTraits):
    # The data to plot
    data = Array()

    # The 4 views displayed
    scene3d = Instance(MlabSceneModel, ())
    scene_x = Instance(MlabSceneModel, ())
    scene_y = Instance(MlabSceneModel, ())
    scene_z = Instance(MlabSceneModel, ())

    # The data source
    data_src3d = Instance(Source)

    # The image plane widgets of the 3D scene
    ipw_3d_x = Instance(PipelineBase)
    ipw_3d_y = Instance(PipelineBase)
    ipw_3d_z = Instance(PipelineBase)

    _axis_names = dict(x=0, y=1, z=2)

    #---------------------------------------------------------------------------
    def __init__(self, **traits):
        super(VolumeSlicer, self).__init__(**traits)
        # Force the creation of the image_plane_widgets:
        self.ipw_3d_x
        self.ipw_3d_y
        self.ipw_3d_z

    #---------------------------------------------------------------------------
    # Default values
    #---------------------------------------------------------------------------
    def _data_src3d_default(self):
        return mlab.pipeline.scalar_field(
            self.data, figure=self.scene3d.mayavi_scene)

    def make_ipw_3d(self, axis_name):
        ipw = mlab.pipeline.image_plane_widget(
            self.data_src3d,
            figure=self.scene3d.mayavi_scene,
            plane_orientation='%s_axes' % axis_name)
        return ipw

    def _ipw_3d_x_default(self):
        return self.make_ipw_3d('x')

    def _ipw_3d_y_default(self):
        return self.make_ipw_3d('y')

    def _ipw_3d_z_default(self):
        return self.make_ipw_3d('z')

    #---------------------------------------------------------------------------
    # Scene activation callbaks
    #---------------------------------------------------------------------------
    @on_trait_change('scene3d.activated')
    def display_scene3d(self):
        outline = mlab.pipeline.outline(
            self.data_src3d,
            figure=self.scene3d.mayavi_scene,
        )
        self.scene3d.mlab.view(40, 50)
        # Interaction properties can only be changed after the scene
        # has been created, and thus the interactor exists
        for ipw in (self.ipw_3d_x, self.ipw_3d_y, self.ipw_3d_z):
            # Turn the interaction off
            ipw.ipw.interaction = 0
        self.scene3d.scene.background = (0, 0, 0)
        # Keep the view always pointing up
        self.scene3d.scene.interactor.interactor_style = tvtk.InteractorStyleTerrain(
        )

    def make_side_view(self, axis_name):
        scene = getattr(self, 'scene_%s' % axis_name)

        # To avoid copying the data, we take a reference to the
        # raw VTK dataset, and pass it on to mlab. Mlab will create
        # a Mayavi source from the VTK without copying it.
        # We have to specify the figure so that the data gets
        # added on the figure we are interested in.
        outline = mlab.pipeline.outline(
            self.data_src3d.mlab_source.dataset,
            figure=scene.mayavi_scene,
        )
        ipw = mlab.pipeline.image_plane_widget(
            outline, plane_orientation='%s_axes' % axis_name)
        setattr(self, 'ipw_%s' % axis_name, ipw)

        # Synchronize positions between the corresponding image plane
        # widgets on different views.
        ipw.ipw.sync_trait('slice_position',
                           getattr(self, 'ipw_3d_%s' % axis_name).ipw)

        # Make left-clicking create a crosshair
        ipw.ipw.left_button_action = 0

        # Add a callback on the image plane widget interaction to
        # move the others
        def move_view(obj, evt):
            position = obj.GetCurrentCursorPosition()
            for other_axis, axis_number in self._axis_names.items():
                if other_axis == axis_name:
                    continue
                ipw3d = getattr(self, 'ipw_3d_%s' % other_axis)
                ipw3d.ipw.slice_position = position[axis_number]

        ipw.ipw.add_observer('InteractionEvent', move_view)
        ipw.ipw.add_observer('StartInteractionEvent', move_view)

        # Center the image plane widget
        ipw.ipw.slice_position = 0.5 * self.data.shape[
            self._axis_names[axis_name]]

        # Position the view for the scene
        views = dict(
            x=(0, 90),
            y=(90, 90),
            z=(0, 0),
        )
        scene.mlab.view(*views[axis_name])
        # 2D interaction: only pan and zoom
        scene.scene.interactor.interactor_style = tvtk.InteractorStyleImage()
        scene.scene.background = (0, 0, 0)

    @on_trait_change('scene_x.activated')
    def display_scene_x(self):
        return self.make_side_view('x')

    @on_trait_change('scene_y.activated')
    def display_scene_y(self):
        return self.make_side_view('y')

    @on_trait_change('scene_z.activated')
    def display_scene_z(self):
        return self.make_side_view('z')

    #---------------------------------------------------------------------------
    # The layout of the dialog created
    #---------------------------------------------------------------------------
    view = View(
        HGroup(
            Group(
                Item(
                    'scene_y',
                    editor=SceneEditor(scene_class=Scene),
                    height=250,
                    width=300),
                Item(
                    'scene_z',
                    editor=SceneEditor(scene_class=Scene),
                    height=250,
                    width=300),
                show_labels=False,
            ),
            Group(
                Item(
                    'scene_x',
                    editor=SceneEditor(scene_class=Scene),
                    height=250,
                    width=300),
                Item(
                    'scene3d',
                    editor=SceneEditor(scene_class=MayaviScene),
                    height=250,
                    width=300),
                show_labels=False,
            ),
        ),
        resizable=True,
        title='Volume Slicer',
    )


m = VolumeSlicer(data=data)
m.configure_traits()

#%% [markdown]
# ## Plot Submission
#
# Submit your plot for your assigned quantum numbers to your Chemistry instructors to get a point for this item.
#%% [markdown]
# ## eDimension Questions
#
# Answer the following questions on eDimension Week 10.
#%% [markdown]
# Q1: In the final function to calculate the hydrogen wave function, you are to use the other previous functions you have calculated. However, some of those functions rounds the result to 5 decimal places. The error on the final wave function magnitude is called ```____``` due to ```____```.
#
# 1. floating point error, rounding error.
# 1. propagation error, rounding error.
# 1. propagation error, floating point error.
# 1. rounding error, propagation error.
#
#%% [markdown]
# Q2: What is the effect when you increase the number of points $N_x$,$N_y$,$N_z$, while maintaining the values the other parameters?
#
# 1. increase of accuracy, decrease of computational time.
# 1. decrease of accuracy, increase of computational time.
# 1. increalse of accuracy, increase of computational time.
# 1. decrease of accuracy, decrease of computational time.
#
#%% [markdown]
# Q3: What is the effect of increasing the distance $r/a$, while maintaining the values
# of the other paramters?
#
# 1. increase of accuracy, no change in computational time.
# 1. decrease of accuracy, change in computational time.
# 1. increalse of accuracy, change in computational time.
# 1. decrease of accuracy, no change in computational time.
#%% [markdown]
# # Example Plots
#
#
# ![wave200.png](attachment:wave200.png)
# Magnitude plot for n = 2, l = 0, m = 0 using contour3d from mlab Mayavi package.
#
# ![wave211.png](attachment:wave211.png)
# Magnitude plot for n = 2, l = 1, m = 1 using contour3d from mlab Mayavi package.
#
# ![wave322.png](attachment:wave322.png)
# Magnitude plot for n = 3, l = 2, m = 2 using contour3d from mlab Mayavi package.
#
