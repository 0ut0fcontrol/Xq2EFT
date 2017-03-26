"""A neat code for tri-linear interpolation in numpy.

source: http://stackoverflow.com/questions/18170782/how-to-build-a-lookup-table-for-tri-linear-interpolation-in-numpy
point order:
        ^
        |z
        ________
       /|1     /|3
      / |     / |
     /__|_ _7/  |
    5|  |____|__|__y__>
     | /0    | /2
     |/      |/
     |_______|
   x/4        6
   v       
   
"""
import numpy as np
coords = np.array([47, 775, 41.3])
ndim = len(coords)
# You would get this with a call to:
# vertices = get_8_points(filename, *coords)
vertices = np.array([[  4.50e+01,   6.00e+02,   4.00e+01,   2.00e-01],
                     [  4.50e+01,   6.00e+02,   4.50e+01,   1.70e+00],
                     [  4.50e+01,   8.00e+02,   4.00e+01,   8.00e-01],
                     [  4.50e+01,   8.00e+02,   4.50e+01,   4.00e-01],
                     [  5.00e+01,   6.00e+02,   4.00e+01,   1.00e-01],
                     [  5.00e+01,   6.00e+02,   4.50e+01,   1.20e+00],
                     [  5.00e+01,   8.00e+02,   4.00e+01,   4.00e-01],
                     [  5.00e+01,   8.00e+02,   4.50e+01,   2.00e-01]])

for dim in range(ndim):
    vtx_delta = 2**(ndim - dim - 1)
    for vtx in range(vtx_delta):
        vertices[vtx, -1] += ((vertices[vtx + vtx_delta, -1] -
                               vertices[vtx, -1]) *
                              (coords[dim] -
                               vertices[vtx, dim]) /
                              (vertices[vtx + vtx_delta, dim] -
                               vertices[vtx, dim]))

print(vertices[0, -1])# prints 0.55075
