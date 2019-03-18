#%% [markdown]
# # Schrödinger Assignment
# Solving the Schrödinger equation for hydrogen atom requires to code for the various mathematical components that constitute the wavefunction. Transfer the codes for required functions you submitted to vocareum corresponding to the weeks in this jupyter notebook.
#%% [markdown]
# ## Week 3:
#%% [markdown]
# Q1: Create two functions to convert degrees to radian and radian to degrees respectively. These functions should take 1 float argument and return the respective conversions each. Round to 5 decimal places.
#
# Hint: you can use Numpy trigonometric function by doing import numpy as np.

#%%
# Code:
import numpy as np
import scipy.constants as c


def deg_to_rad(deg):
    return round(deg / 180 * np.pi, 5)


def rad_to_deg(rad):
    return round(rad / np.pi * 180, 5)


# Test:
# assert statement will throw error if the result is wrong
# no output will be produced for correct results

eps = 1e-3
assert abs(deg_to_rad(90) - 1.5708) < eps
assert abs(deg_to_rad(180) - 3.14159) < eps
assert abs(deg_to_rad(270) - 4.71239) < eps
assert abs(rad_to_deg(3.14) - 179.90875) < eps
assert abs(rad_to_deg(3.14 / 2.0) - 89.95437) < eps
assert abs(rad_to_deg(3.14 * 3 / 4) - 134.93156) < eps

#%% [markdown]
# Q2: Create two functions to convert spherical to cartesian coordinates and cartesian to spherical coordinates. These functions should take 3 float arguments and return the 3 respective conversions. Round to 5 decimal places. (Pre-requisite: Part 2b) The convention is shown below.
# Hint: You can use Numpy trigonometric function by doing ```import numpy as np```.
#
# ![spherical.png](attachment:spherical.png)


#%%
###Code:
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


# Test:
# assert statement will throw error if the result is wrong
# no output will be produced for correct results

assert spherical_to_cartesian(3, 0, np.pi) == (-0.0, 0.0, 3.0)
assert spherical_to_cartesian(3, np.pi / 2.0, np.pi / 2.0) == (0.0, 3.0, 0.0)
assert spherical_to_cartesian(3, np.pi, 0) == (0.0, 0.0, -3.0)
assert cartesian_to_spherical(3, 0, 0) == (3.0, 1.5708, 0.0)
assert cartesian_to_spherical(0, 3, 0) == (3.0, 1.5708, 1.5708)
assert cartesian_to_spherical(0, 0, 3) == (3.0, 0.0, 0.0)
assert cartesian_to_spherical(0, -3, 0) == (3.0, 1.5708, -1.5708)
assert cartesian_to_spherical(0, 0, 3) == (3.0, 0.0, 0.0)

#%% [markdown]
# ## Week4:
#%% [markdown]
# Q1: Create a function that calculates the normalized angular solution. This function should take 4 float arguments and return the value of the normalized angular solution for the specific m, l, θ and Φ. The return value is a complex number rounded to 5 decimal places for both the real and the imaginary parts.
#
# Hint: You may want to use ```np.round()``` function to round the return value to 5 decimal places. You can also use ```np.complex()``` to ensure the output is a complex number.

#%%
### Code for angular wavefunction:
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


# Test for angular wavefunction:
# assert statement will throw error if the result is wrong
# no output will be produced for correct results

assert angular_wave_func(0, 0, 0, 0) == (0.28209 + 0j)
assert angular_wave_func(0, 1, c.pi, 0) == (-0.4886 + 0j)
assert angular_wave_func(1, 1, c.pi / 2, c.pi) == (0.34549 - 0j)
assert angular_wave_func(0, 2, c.pi, 0) == (0.63078 + 0j)

#%% [markdown]
# Q2: Create a function that calculates the normalized radial solution. This function should take 3 float arguments and return the value of the normalized radial solution for the specific n, l, and r.  The return value should be normalized to a^(-3/2), where a is the Bohr's radius, and np.rounded to 5 decimal places.
#
# Hint: You may want to use ```np.np.round()``` function to np.round the return value to 5 decimal places.

#%%
# Code for radial wavefunction:
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


# Test for radial wavefunction:
# assert statement will throw error if the result is wrong
# no output will be produced for correct results

assert radial_wave_func(1, 0, a) == 0.73576
assert radial_wave_func(1, 0, a) == 0.73576
assert radial_wave_func(2, 1, a) == 0.12381
assert radial_wave_func(2, 1, 2 * a) == 0.15019
assert radial_wave_func(3, 1, 2 * a) == 0.08281

#%% [markdown]
# ## Week 5:
#%% [markdown]
# Q1: Create a function called mgrid that takes in six arguments xstart, xend, xpoints, ystart, yend, ypoints. The first three arguments specifies the starting, ending and the number of points in the x axis. The last three arguments does the same for the y axis. The function should return a list of lists as described in numpy.mgrid when the step length is a complex number. You are not allowed to use numpy.mgrid or any other built-in function in your code, but you are strongly suggested to use numpy.mgrid to test your version of mgrid.


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


mgrid2d(0, 4, 5, 0, 4, 3)
#%%
# Test:
# assert statement will throw error if the result is wrong
# no output will be produced for correct results
import numpy as np

assert np.shape(mgrid2d(0, 4, 5, 0, 4,
                        5)) == np.shape(np.mgrid[0:4:5j, 0:4:5j])
assert np.allclose(mgrid2d(0, 4, 5, 0, 4, 5), np.mgrid[0:4:5j, 0:4:5j])

assert np.shape(mgrid2d(0, 5, 15, 0, 4,
                        10)) == np.shape(np.mgrid[0:5:15j, 0:4:10j])
assert np.allclose(mgrid2d(0, 5, 15, 0, 4, 10), np.mgrid[0:5:15j, 0:4:10j])

#%% [markdown]
# Q2. Create a function called mgrid that takes in nine arguments, three to specify each x, y, and z axis. The first three input arguments specifies the start (xstart), end (xend), and the number of points (xpoints) in the x axis. Similarly for the y and z axis. The function should return a list of lists as described in numpy.mgrid. You are not allowed to use numpy.mgrid or any other built-in function.
#
# However, you can use numpy.mgrid to test your own function and compare the result.


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


mgrid3d(1, 2, 3, 3, 4, 3, 5, 6, 4)

#%%
np.mgrid[1:2:3j, 3:4:3j, 5:6:4j]
#%%

# Test:
# assert statement will throw error if the result is wrong
# no output will be produced for correct results

assert np.shape(mgrid3d(0, 4, 5, 0, 4, 5, 0, 4,
                        5)) == np.shape(np.mgrid[0:4:5j, 0:4:5j, 0:4:5j])
assert np.allclose(
    mgrid3d(0, 4, 5, 0, 4, 5, 0, 4, 5), np.mgrid[0:4:5j, 0:4:5j, 0:4:5j])

assert np.shape(mgrid3d(0, 5, 15, 0, 4, 10, 1, 2,
                        3)) == np.shape(np.mgrid[0:5:15j, 0:4:10j, 1:2:3j])
assert np.allclose(
    mgrid3d(0, 5, 15, 0, 4, 10, 1, 2, 3), np.mgrid[0:5:15j, 0:4:10j, 1:2:3j])

#%% [markdown]
# ## Week 6:
#%% [markdown]
# Q1. Create a function that calculates the square of the magnitude of the real wave function. The function takes in several arguments:
# * ```n```: quantum number n
# * ```l```: quantum number l
# * ```m```: quantum number m
# * ```roa```: maximum distance to plot from the centre, normalized to Bohr radius, i.e. r/a.
# * ```Nx```: Number of points in the x axis.
# * ```Ny```: Number of points in the y axis.
# * ```Nz```: Number of points in the z axis.
#
# The function should return:
# * ```xx```: x location of all the points in a 3D Numpy array.
# * ```yy```: y location of all the points in a 3D Numpy array.
# * ```zz```: z location of all the points in a 3D Numpy array.
# * ```density```: The square of the magnitude of the real wave function, i.e. |Ψ|^2
#
# To obtain the real wave function from complex parts, i.e. m != 0:
#
# ![spherical_harmonics_real.jpg](attachment:spherical_harmonics_real.jpg)
#
# Hint: You may find the following functions to be useful:
# * ```fvec=numpy.vectorize(f)```: This function takes in a function and return its vectorized version of the function.
#
# * ```m=mag(c)```: This function takes in a complex number and returns its absolute value or its magnitude. Use your own function rather than numpy's built-in function.
#
# * ```ar=numpy.array(x)```: This function takes in a list and returns a numpy array. Numpy array is faster to process than Python's list.
#
# Hint: You will need to use all the previous functions you have done. Note that some of those functions may round the output to 5 decimal places and the final magnitude output from this function should also be rounded to 5 decimal places.

#%%
###Code:
import numpy as np


def hydrogen_wave_func(n, l, m, roa, Nx, Ny, Nz):
    xx, yy, zz = np.array(mgrid3d(-roa, roa, Nx, -roa, roa, Ny, -roa, roa, Nz))

    # '''vectorize cartesian to spherical function'''
    # cartesian_to_spherical_vectorized = np.vectorize(cartesian_to_spherical)
    # r, theta, phi = cartesian_to_spherical_vectorized(xx, yy, zz)
    r, theta, phi = np.vectorize(cartesian_to_spherical)(xx, yy, zz)

    ang = np.vectorize(angular_wave_func)(m, l, theta, phi)
    rad = np.vectorize(radial_wave_func)(
        n, l, r * a)  # radius multiply by a (Bohr radius)
    real_wave_func = np.absolute(ang * rad)**2
    return tuple(
        np.round([xx, yy, zz, real_wave_func],
                 5))  # use np's round function, convert nparray to tuple


print(hydrogen_wave_func(2, 1, 1, 5, 3, 4, 2))

#%%
###Test:
print('Test 1')
x, y, z, mag = hydrogen_wave_func(2, 1, 1, 8, 2, 2, 2)
print('x, y, z:')
print(x, y, z)
print('mag:')
print(mag)

print('\n')
print('Test 2')
x, y, z, mag = hydrogen_wave_func(2, 1, 1, 5, 3, 4, 2)
print('x, y, z:')
print(x, y, z)
print('mag:')
print(mag)

print('\n')
print('Test 3')
x, y, z, mag = hydrogen_wave_func(2, 0, 0, 3, 5, 4, 3)
print('x, y, z:')
print(x, y, z)
print('mag:')
print(mag)

#%% [markdown]
# Expected Output:
#
# Test 1
#
# x, y, z:
#
# [[[-8. -8.]
#
#   [-8. -8.]]
#
#  [[ 8.  8.]
#
#   [ 8.  8.]]] [[[-8. -8.]
#
#   [ 8.  8.]]
#
#  [[-8. -8.]
#
#   [ 8.  8.]]] [[[-8.  8.]
#
#   [-8.  8.]]
#
#  [[-8.  8.]
#
#   [-8.  8.]]]
#
# mag:
#
# [[[ 0.  0.]
#
#   [ 0.  0.]]
#
#  [[ 0.  0.]
#
#   [ 0.  0.]]]
#
#
# Test 2
#
# x, y, z:
#
# [[[-5. -5.]
#
#   [-5. -5.]
#
#   [-5. -5.]
#
#   [-5. -5.]]
#
#  [[ 0.  0.]
#
#   [ 0.  0.]
#
#   [ 0.  0.]
#
#   [ 0.  0.]]
#
#  [[ 5.  5.]
#
#   [ 5.  5.]
#
#   [ 5.  5.]
#
#   [ 5.  5.]]] [[[-5.      -5.     ]
#
#   [-1.66667 -1.66667]
#
#   [ 1.66667  1.66667]
#
#   [ 5.       5.     ]]
#
#  [[-5.      -5.     ]
#
#   [-1.66667 -1.66667]
#
#   [ 1.66667  1.66667]
#
#   [ 5.       5.     ]]
#
#  [[-5.      -5.     ]
#
#   [-1.66667 -1.66667]
#
#   [ 1.66667  1.66667]
#
#   [ 5.       5.     ]]] [[[-5.  5.]
#
#   [-5.  5.]
#
#   [-5.  5.]
#
#   [-5.  5.]]
#
#  [[-5.  5.]
#
#   [-5.  5.]
#
#   [-5.  5.]
#
#   [-5.  5.]]
#
#  [[-5.  5.]
#
#   [-5.  5.]
#
#   [-5.  5.]
#
#   [-5.  5.]]]
#
# mag:
#
# [[[  4.00000000e-05   4.00000000e-05]
#
#   [  1.70000000e-04   1.70000000e-04]
#
#   [  1.70000000e-04   1.70000000e-04]
#
#   [  4.00000000e-05   4.00000000e-05]]
#
#  [[  0.00000000e+00   0.00000000e+00]
#
#   [  0.00000000e+00   0.00000000e+00]
#
#   [  0.00000000e+00   0.00000000e+00]
#
#   [  0.00000000e+00   0.00000000e+00]]
#
#  [[  4.00000000e-05   4.00000000e-05]
#
#   [  1.70000000e-04   1.70000000e-04]
#
#   [  1.70000000e-04   1.70000000e-04]
#
#   [  4.00000000e-05   4.00000000e-05]]]
#
#
# Test 3
#
# x, y, z:
#
# [[[-3.  -3.  -3. ]
#
#   [-3.  -3.  -3. ]
#
#   [-3.  -3.  -3. ]
#
#   [-3.  -3.  -3. ]]
#
#  [[-1.5 -1.5 -1.5]
#
#   [-1.5 -1.5 -1.5]
#
#   [-1.5 -1.5 -1.5]
#
#   [-1.5 -1.5 -1.5]]
#
#  [[ 0.   0.   0. ]
#
#   [ 0.   0.   0. ]
#
#   [ 0.   0.   0. ]
#
#   [ 0.   0.   0. ]]
#
#  [[ 1.5  1.5  1.5]
#
#   [ 1.5  1.5  1.5]
#
#   [ 1.5  1.5  1.5]
#
#   [ 1.5  1.5  1.5]]
#
#  [[ 3.   3.   3. ]
#
#   [ 3.   3.   3. ]
#
#   [ 3.   3.   3. ]
#
#   [ 3.   3.   3. ]]] [[[-3. -3. -3.]
#
#   [-1. -1. -1.]
#
#   [ 1.  1.  1.]
#
#   [ 3.  3.  3.]]
#
#  [[-3. -3. -3.]
#
#   [-1. -1. -1.]
#
#   [ 1.  1.  1.]
#
#   [ 3.  3.  3.]]
#
#  [[-3. -3. -3.]
#
#   [-1. -1. -1.]
#
#   [ 1.  1.  1.]
#
#   [ 3.  3.  3.]]
#
#  [[-3. -3. -3.]
#
#   [-1. -1. -1.]
#
#   [ 1.  1.  1.]
#
#   [ 3.  3.  3.]]
#
#  [[-3. -3. -3.]
#
#   [-1. -1. -1.]
#
#   [ 1.  1.  1.]
#
#   [ 3.  3.  3.]]] [[[-3.  0.  3.]
#
#   [-3.  0.  3.]
#
#   [-3.  0.  3.]
#
#   [-3.  0.  3.]]
#
#  [[-3.  0.  3.]
#
#   [-3.  0.  3.]
#
#   [-3.  0.  3.]
#
#   [-3.  0.  3.]]
#
#  [[-3.  0.  3.]
#
#   [-3.  0.  3.]
#
#   [-3.  0.  3.]
#
#   [-3.  0.  3.]]
#
#  [[-3.  0.  3.]
#
#   [-3.  0.  3.]
#
#   [-3.  0.  3.]
#
#   [-3.  0.  3.]]
#
#  [[-3.  0.  3.]
#
#   [-3.  0.  3.]
#
#   [-3.  0.  3.]
#
#   [-3.  0.  3.]]]
#
# mag:
#
# [[[  5.60000000e-04   7.20000000e-04   5.60000000e-04]
#
#   [  7.10000000e-04   5.70000000e-04   7.10000000e-04]
#
#   [  7.10000000e-04   5.70000000e-04   7.10000000e-04]
#
#   [  5.60000000e-04   7.20000000e-04   5.60000000e-04]]
#
#  [[  6.90000000e-04   6.40000000e-04   6.90000000e-04]
#
#   [  6.80000000e-04   6.00000000e-05   6.80000000e-04]
#
#   [  6.80000000e-04   6.00000000e-05   6.80000000e-04]
#
#   [  6.90000000e-04   6.40000000e-04   6.90000000e-04]]
#
#
#  [[  7.20000000e-04   5.00000000e-04   7.20000000e-04]
#
#   [  5.70000000e-04   3.66000000e-03   5.70000000e-04]
#
#   [  5.70000000e-04   3.66000000e-03   5.70000000e-04]
#
#   [  7.20000000e-04   5.00000000e-04   7.20000000e-04]]
#
#  [[  6.90000000e-04   6.40000000e-04   6.90000000e-04]
#
#   [  6.80000000e-04   6.00000000e-05   6.80000000e-04]
#
#   [  6.80000000e-04   6.00000000e-05   6.80000000e-04]
#
#   [  6.90000000e-04   6.40000000e-04   6.90000000e-04]]
#
#  [[  5.60000000e-04   7.20000000e-04   5.60000000e-04]
#
#   [  7.10000000e-04   5.70000000e-04   7.10000000e-04]
#
#   [  7.10000000e-04   5.70000000e-04   7.10000000e-04]
#
#   [  5.60000000e-04   7.20000000e-04   5.60000000e-04]]]
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
# Mayavi code:

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
###Volume slicer code:
import numpy as np

from traits.api import HasTraits, Instance, Array, on_trait_change
from traitsui.api import View, Item, HGroup, Group

from tvtk.api import tvtk
from tvtk.pyface.scene import Scene

from mayavi import mlab
from mayavi.core.api import PipelineBase, Source
from mayavi.core.ui.api import SceneEditor, MayaviScene, MlabSceneModel

################################################################################
# Create some data
data = np.load('den_test.dat')


################################################################################
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
