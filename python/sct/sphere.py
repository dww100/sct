#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Sphere model related functions used in the SCT suite of programs.
Includes functions to calculate scattering curve from sphere models
"""

# Copyright 2014 University College London

# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at

#    http://www.apache.org/licenses/LICENSE-2.0

# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import numpy as np
import scipy.optimize as optimize
import scipy.spatial.distance as dist
from bisect import bisect_right

def create_grid(atom_coords, x_axis, y_axis, z_axis):
    """Count atoms in grid boxes defined by passed axes intervals

    @type  atom_coords: list
    @param atom_coords: A list containing lists of x, y & z coordinates
                        (3 * floats)
    @type  x_axis:      numpy array
    @param x_axis:      Array of floats defining the intersection of grid boxes
                        along the x-axis
    @type  y_axis:      numpy array
    @param y_axis:      Array of floats defining the intersection of grid boxes
                        along the y-axis
    @type  z_axis:      numpy array
    @param z_axis:      Array of floats defining the intersection of grid boxes
                        along the z-axis
    @rtype:             numpy array
    @return:            3 dimensional array representing the number of atoms in
                        each grid box defined by the input axes divisions
    """

    # Initialize grid
    # dimension -1 as looking at intervals rather than axes points
    grid = np.zeros([len(x_axis) - 1, len(y_axis) - 1, len(z_axis) - 1])

    # Use a binary search to find which interval each atom coordinate fits in
    # then increment count in appropriate grid box
    for atom_coord in atom_coords:
        x = bisect_right(x_axis, atom_coord[0]) - 1
        y = bisect_right(y_axis, atom_coord[1]) - 1
        z = bisect_right(z_axis, atom_coord[2]) - 1
        grid[x, y, z] += 1

    return grid

def define_grid_axes(coords, box_side):
    """
    Create intervals of length box_side along x, y & z axes. The axes extend to
    include all input atomic coordinates.

    @type  coords:    list
    @param coords:    A list containing lists of x, y & z coordinates
                      (3 * floats)
    @type  box_side:  float
    @param box_side:  Length of the sides of the cubes used to create a grid to
                      coarse grain atomic structure into a sphere model
    @rtype:           numpy array, numpy array, numpy array
    @return:          Intervals of length box_side along the x, y and z axes.
                      Minimum and maximum values chosed to extend beyond the
                      input coordinates.
    """

    # Create limits of the axes which contain all coords
    coord_mins = np.floor(np.amin(coords, 0))
    coord_maxs = np.ceil(np.amax(coords, 0))

    # Create float valued arrays of intervals
    x_axis = np.arange(coord_mins[0] - 2.0 * box_side,
                       coord_maxs[0] + 2.0 * box_side, box_side)
    y_axis = np.arange(coord_mins[1] - 2.0 * box_side,
                       coord_maxs[1] + 2.0 * box_side, box_side)
    z_axis = np.arange(coord_mins[2] - 2.0 * box_side,
                       coord_maxs[2] + 2.0 * box_side, box_side)

    return x_axis, y_axis, z_axis

def grid_to_spheres(grid, radius, cutoff, x_axis, y_axis, z_axis):
    """
    Generate a list of sphere coordinates based on counts of atoms in the passed
    grid. A sphere is included in the model in each grid space with a content
    >= cutoff. Coordinates are based on the passed axes defining teh grid boxes.

    @type  grid:    numpy array
    @param grid:    3 dimensional array representing the number of atoms in
                    each grid box defined by the input axes divisions.
    @type  radius:  float
    @param radius:  Radius of spheres within the sphere model.
    @type  cutoff:  integer
    @param cutoff:  Cutoff for number of atoms within a grid square before a
                    a sphere is placed in that grid box.
    @type  x_axis:  numpy array
    @param x_axis:  Array of floats defining the intersection of grid boxes
                    along the x-axis.
    @type  y_axis:  numpy array
    @param y_axis:  Array of floats defining the intersection of grid boxes
                    along the y-axis.
    @type  z_axis:  numpy array
    @param z_axis:  Array of floats defining the intersection of grid boxes
                    along the z-axis.
    @rtype:         list
    @return:        A list containing lists of x, y & z coordinates
                    (3 * floats) of spheres within the sphere model.
    """

    # Find the locations in the grid containing more than cutoff atoms
    # A sphere will be placed at the centre of each of these cubes
    full = np.where(grid >= cutoff)

    sphere_coords = []

    no_spheres = len(full[0])

    for ndx in xrange(0, no_spheres):
        x = x_axis[full[0][ndx]] + radius
        y = y_axis[full[1][ndx]] + radius
        z = z_axis[full[2][ndx]] + radius
        sphere_coords.append([x,y,z])

    return sphere_coords

def create_sphere_model(atom_coords, cutoff, box_side, **kwargs):
    """
    Create a sphere model from the input atomic coordinates. First a grid is
    created with divisions of length box_side. Then the number of atoms in each
    grid box is counted. A sphere is then included in the final model placed at
    the centre of a grid box if more that cutoff atoms are within that box.

    @type  atom_coords: list
    @param atom_coords: A list containing lists of x, y & z coordinates
                        (3 * floats)
    @type  cutoff:      integer
    @param cutoff:      Cutoff for number of atoms within a grid square before a
                        a sphere is placed in that grid box.
    @type  box_side:    float
    @param box_side:    Length of the sides of the cubes used to create a grid
                        to coarse grain atomic structure into a sphere model.
    @keyword x_axis:    Array of floats defining the intersection of grid boxes
                        along the x-axis.
    @keyword y_axis:    Array of floats defining the intersection of grid boxes
                        along the y-axis.
    @keyword z_axis:    Array of floats defining the intersection of grid boxes
                        along the z-axis.
    @rtype:             list, numpy array, numpy array, numpy array
    @return:            
                        - A list containing lists of x, y & z coordinates
                        (3 * floats) of spheres within the sphere model.

                        - Array of floats defining the intersection of grid
                        boxes along the x-axis.

                        - Array of floats defining the intersection of grid
                        boxes along the y-axis.

                        - Array of floats defining the intersection of grid
                        boxes along the z-axis.
    """

    if ('xaxis' in kwargs) and ('yaxis' in kwargs) and ('zaxis' in kwargs):
        x_axis = kwargs['xaxis']
        y_axis = kwargs['yaxis']
        z_axis = kwargs['zaxis']
    else:
        # Get axes for the grid used to create sphere models
        # The grid we create contains all atoms and is divided into cubes with
        # dimension box_side
        x_axis, y_axis, z_axis = define_grid_axes(atom_coords, box_side)

    # Conversion to numpy array speeds up operations later
    atom_coords = np.array(atom_coords)

    # Create a grid containing the number of atoms in each cube defined
    # by the axes created above
    grid = create_grid(atom_coords, x_axis, y_axis, z_axis)

    # Convert grid to spheres
    # A sphere centre of a box is created if > cutoff atoms are within it
    sphere_coords = grid_to_spheres(grid, box_side/2.0, cutoff,
                                    x_axis, y_axis, z_axis)

    return sphere_coords, x_axis, y_axis, z_axis

def residual2_box(box_side, cutoff, atom_coords, targ_volume):
    """
    Compute squared residual of sphere model to the theroetical target volume

    @type  box_side:     float
    @param box_side:     Length of the sides of the cubes used to create a grid
                         to coarse grain atomic structure into a sphere model.
    @type  cutoff:       integer
    @param cutoff:       Cutoff for number of atoms within a grid square before
                         a sphere is placed in that grid box.
    @type  atom_coords:  list
    @param atom_coords:  A list containing lists of x, y & z coordinates
                         (3 * floats)
    @type  targ_volume:   float
    @param targ_volume:   The target volume for the dry model calculated from
                          the protein sequence.
    @rtype:               float
    @return:              Squared difference between sphere model volume and
                          target volume.
    """

    # Create sphere model
    sphere_coords, x_axis, y_axis, z_axis = create_sphere_model(atom_coords,
                                                               cutoff, box_side)
    # Compute volume of the model
    vol_spheres = len(sphere_coords) * box_side**3

    # Return squared residual
    return (vol_spheres - targ_volume)**2

def optimize_box_side(cutoff, coords, targ_vol, side_min, side_max, tolerance):
    """
    Get optimal box_side to reproduce target_volume (and deviation)

    @type  cutoff:     integer
    @param cutoff:     Cutoff for number of atoms within a grid square before
                       a sphere is placed in that grid box.
    @type  coords:     list
    @param coords:     A list containing lists of x, y & z coordinates
                       (3 * floats)
    @type  targ_vol:   float
    @param targ_vol:   The target volume for the dry model calculated from
                       the protein sequence.
    @type  side_min:   float
    @param side_min:   Minimum side length to consider in the optimization.
    @type  side_max:   float
    @param side_min:   Maximum side length to consider in the optimization.
    @type  tolerance:  float
    @param tolerance:  Percentage tolerance before optimization is said to
                       have converged.
    @rtype:            float
    @return:           Squared difference between sphere model volume and
                       target volume.
    """

    # Format the bounds within which to run the optimization
    side_bounds = (side_min, side_max)

    # Minimize the squared residuals between the sphere model ang target volume
    # Uses minimize_scalar from scipy.optimize
    opt = optimize.minimize_scalar(residual2_box, args=(cutoff, coords, targ_vol),
                                   bounds=side_bounds, method='bounded',
                                   options={'xtol' : tolerance})

    # Return the optimized box_side and residual
    return opt['x'], opt['fun']**0.5

def write_sphere_line(x, y, z, radius, out):
    """
    Write out a single sphere entry to a SCT sphere file format

    @type  x:       float
    @param x:       X coordinate of sphere
    @type  y:       float
    @param y:       Y coordinate of sphere
    @type  z:       float
    @param z:       Z coordinate of sphere
    @type  radius:  float
    @param radius:  radius of sphere
    @type  out:     file object
    @param out:     The file to which the output sphere line is written
    """

    # Add radius to coordinate to centre sphere in grid box

    xs = x + radius
    ys = y + radius
    zs = z + radius

    out.write("{0:10.2f}{1:10.2f}{2:10.2f}{3:10.2f}\n".format(xs, ys, zs, radius))

def write_spheres(coords, radius, out):
    """
    Write coordinates and radius to a SCT sphere file format for each sphere
    location in coords

    @type  coords:  list
    @param coords:  A list containing lists of x, y & z coordinates
                    (3 * floats) of spheres within the sphere model.
    @type  radius:  float
    @param radius:  Radius of spheres within the sphere model.
    @type  out:     file object
    @param out:     The file to which the output sphere line is written
    """
    
    for coord in coords:
        write_sphere_line(coord[0], coord[1], coord[2], radius, out)

def parse_sphere_line(line):
    """
    Parse line from sphere file, return dictionary with list of
    'coords' and 'radius'

    @type  line:  string
    @param line:  Line from a saved SCT format sphere model file
    @rtype:       dictionary
    @return:      Dictionary with:

                  - coords: A list containing x, y & z coordinates (3 * floats)
                  of a sphere within the sphere model.
                  - radius: Sphere radius (float)
    """

    data = {}

    # SCT sphere files have four columns of fixed width (3 coords + radius)

    data['coords'] = [float(line[0:10]), float(line[11:20]), float(line[21:30])]
    data['radius'] = float(line[31:40])

    return data

def read_mono_spheres(filename):
    """
    Read in SCT format sphere file and return coordinates and radius
    Note: assumed that we have only one radius of sphere here

    @type  filename:  string
    @param filename:  Path to the input SCT format sphere model file
    @rtype:           list, float
    @return:    
                    - A list containing lists of x, y & z coordinates
                    (3 * floats) of a sphere within the sphere model.
                    - Sphere radius (float)
    """

    spheres = []

    with open(filename) as f:
        for line in f:
            sphere = parse_sphere_line(line)
            spheres.append(sphere['coords'])
            radius = sphere['radius']

    return spheres, radius

def sphere_model_rg(coords, radius):
    """
    Calculate the radius of gyration of a set of spheres given their position
    and radius

    @type  coords:    list
    @param coords:    A list containing lists of x, y & z coordinates
                      (3 * floats)
    @type  radius:    float
    @param radius:    Sphere radius
    @rtype:           float
    @return:          Radius of gyration of input sphere distribution
    """

    # Radius of gyration of a cube = side length
    # Cube side length = 2 * radius of cube that fits inside it
    box_rg2 = (2 * radius)**2

    # R_g**2 = 1/N sum(r_k - r_mean)**2
    # where r_k is the position of each sphere and r_k their mean position

    no_spheres = len(coords)

    coords = np.array(coords)
    x_mean = np.mean(coords[ : , 0 ])
    y_mean = np.mean(coords[ : , 1 ])
    z_mean = np.mean(coords[ : , 2 ])

    rad_sq = 0.0

    for ii in range(0, no_spheres):
        d = (coords[ii,0] - x_mean)**2 + (coords[ii,1] - y_mean)**2 + (coords[ii,2] - z_mean)**2
        rad_sq += d

    return np.sqrt(rad_sq / no_spheres + box_rg2)

def hydrate_sphere(coords, hydration_pos, dist):
    """
    Add test hydration spheres to selected positions surrounding the input
    sphere coordinate

    @type  coords:         list
    @param coords:         List containing x, y & z coordinates (3 * floats) for
                           a sphere.
    @type  hydration_pos:  list
    @param hydration_pos:  List of positions to include in wet sphere model:
                           Position 1 = original sphere position,
                           2 to 27 positions on cube centred on original sphere.
    @type  dist:           float
    @param dist:           Distance between sphere centres ( 2 * sphere radius)
    @rtype:                list
    @return:               List containing coordinates of all spheres in
                           hydrated sphere model.
    """

    hyd_coords = []

    hyd_coords.append(coords)

    if 1 in hydration_pos:
        hyd_coords.append([coords[0] - dist, coords[1], coords[2]])

    if 2 in hydration_pos:
        hyd_coords.append([coords[0] + dist, coords[1], coords[2]])

    if 3 in hydration_pos:
        hyd_coords.append([coords[0], coords[1], coords[2] + dist])

    if 4 in hydration_pos:
        hyd_coords.append([coords[0], coords[1], coords[2] - dist])

    if 5 in hydration_pos:
        hyd_coords.append([coords[0], coords[1] + dist, coords[2]])

    if 6 in hydration_pos:
        hyd_coords.append([coords[0], coords[1] - dist, coords[2]])


    if 7 in hydration_pos:
        hyd_coords.append([coords[0] - dist, coords[1], coords[2] + dist])

    if 8 in hydration_pos:
        hyd_coords.append([coords[0] + dist, coords[1], coords[2] - dist])

    if 9 in hydration_pos:
        hyd_coords.append([coords[0], coords[1] + dist, coords[2] + dist])

    if 10 in hydration_pos:
        hyd_coords.append([coords[0], coords[1] + dist, coords[2] - dist])

    if 11 in hydration_pos:
        hyd_coords.append([coords[0] - dist, coords[1] - dist, coords[2]])

    if 12 in hydration_pos:
        hyd_coords.append([coords[0] + dist, coords[1] - dist, coords[2]])


    if 13 in hydration_pos:
        hyd_coords.append([coords[0] - dist, coords[1], coords[2] - dist])

    if 14 in hydration_pos:
        hyd_coords.append([coords[0] + dist, coords[1], coords[2] + dist])

    if 15 in hydration_pos:
        hyd_coords.append([coords[0], coords[1] - dist, coords[2] - dist])

    if 16 in hydration_pos:
        hyd_coords.append([coords[0], coords[1] - dist, coords[2] + dist])

    if 17 in hydration_pos:
        hyd_coords.append([coords[0] - dist, coords[1] + dist, coords[2]])

    if 18 in hydration_pos:
        hyd_coords.append([coords[0] + dist, coords[1] + dist, coords[2]])


    if 19 in hydration_pos:
        hyd_coords.append([coords[0] + dist, coords[1] + dist, coords[2] + dist])

    if 20 in hydration_pos:
        hyd_coords.append([coords[0] + dist, coords[1] - dist, coords[2] + dist])

    if 21 in hydration_pos:
        hyd_coords.append([coords[0] - dist, coords[1] + dist, coords[2] - dist])

    if 22 in hydration_pos:
        hyd_coords.append([coords[0] - dist, coords[1] - dist, coords[2] - dist])

    if 23 in hydration_pos:
        hyd_coords.append([coords[0] + dist, coords[1] + dist, coords[2] - dist])

    if 24 in hydration_pos:
        hyd_coords.append([coords[0] - dist, coords[1] - dist, coords[2] + dist])


    if 25 in hydration_pos:
        hyd_coords.append([coords[0] - dist, coords[1] + dist, coords[2] + dist])

    if 26 in hydration_pos:
        hyd_coords.append([coords[0] + dist, coords[1] - dist, coords[2] - dist])

    return hyd_coords

def hydrate_spheres(coords, hydration_pos, radius):
    """
    Add hydration layer spheres to the spheres in the input coordinate list

    @type  coords:         list
    @param coords:         List of lists containing x, y & z coordinates
                           (3 * floats) for each sphere in a sphere model.
    @type  hydration_pos:  list
    @param hydration_pos:  List of positions to include in wet sphere model:
                           Position 1 = original sphere position,
                           2 to 27 positions on cube centred on original sphere.
    @type  radius:         float
    @param radius:         Sphere radius.
    @rtype:                list
    @return:               A list containing lists of x, y & z coordinates
                           (3 * floats) for all test spheres in hydrated model.

                           B{Note}: This includes original spheres and
                           overlapping hydration spheres as no filtering 
                           has yet occured.
    """

    wet_spheres = []

    for coord in coords:
        # Separation between sphere centres is 2 * radius of a sphere
        wet_spheres += hydrate_sphere(coord, hydration_pos, 2.0 * radius)

    return wet_spheres

def hydrate_sphere_model(dry_spheres, hydration_no, box_side, water_cutoff, **kwargs):
    """
    Create hydrated sphere model in 4 steps (returning the appropriate sphere
    coordinates):

        1. Add spheres to the dry model in positions assigned through 
        hydration_no (Position 1 = original sphere position, 2 to 27 
        positions on cube centred on original sphere).
        2. Filter out excess spheres using create_sphere_model with water_cutoff
        (this cutoff should be high ~ 10-12 spheres per box)
        3. Add dry spheres back to the model. Some will have been removed 
        alongside the excess water in the last step
        4. Filter out overlapping spheres using create_sphere_model with cutoff
        set to 1.

    @type  dry_spheres:   list
    @param dry_spheres:   A list containing lists of x, y & z coordinates
                          (3 * floats) for all spheres in an non-hydrated sphere
                          model.
    @type  hydration_no:  integer
    @param hydration_no:  Number defining which positions on a cube surrounding
                          each input sphere should be used to place test
                          hydration spheres. Position 1 = original sphere
                          position, 2 to 27 positions on cube centred on
                          original sphere.
    @type  box_side:      float
    @param box_side:      Length of the sides of the cubes used to create a grid
                          to coarse grain atomic structure into a sphere model
                          (and here to place hydration spheres).
    @type  water_cutoff:  integer
    @param water_cutoff:  Cutoff used to remove overlapping and excess test
                          hydration spheres
    @keyword xaxis:       Array of floats defining the intersection of grid
                          boxes along the x-axis (from creation of the original
                          sphere model).
    @keyword yaxis:       Array of floats defining the intersection of grid
                          boxes along the y-axis (from creation of the original
                          sphere model).
    @keyword zaxis:       Array of floats defining the intersection of grid
                          boxes along the z-axis (from creation of the original
                          sphere model).
    @rtype:               list
    @return:              A list containing lists of x, y & z coordinates
                          (3 * floats) for all spheres in the hydrated sphere
                          model.
    """

    # Check to see if axes divisions from the creation of the dry model were
    # passed in.
    # If so use these when creating the sphere model, otherwise they will need
    # to be calculated
    if ('xaxis' in kwargs) and ('yaxis' in kwargs) and ('zaxis' in kwargs):
        old_axes = True
        x_axis = kwargs['xaxis']
        y_axis = kwargs['yaxis']
        z_axis = kwargs['zaxis']
    else:
        old_axes = False

    radius = box_side / 2.0

    # List of positions to include in wet sphere model:
    # Position 0 = original sphere position,
    # Positions 1 to 26 positions on cube centred on original sphere
    hydration_pos = xrange(0, hydration_no + 1)

    # Create hydrated model
    wet_spheres = hydrate_spheres(dry_spheres, hydration_pos, radius)

    if old_axes:
        wet_spheres, x_axis, y_axis, z_axis = create_sphere_model(wet_spheres,
                                                                  water_cutoff,
                                                                  box_side,
                                                                  xaxis = x_axis,
                                                                  yaxis = y_axis,
                                                                  zaxis = z_axis)
    else:
        wet_spheres, x_axis, y_axis, z_axis = create_sphere_model(wet_spheres,
                                                                  water_cutoff,
                                                                  box_side)

    wet_spheres = wet_spheres + dry_spheres

    # Remove un-necessary/overlapping spheres from model
    # Use grid system from pdb2sphere with a cutoff of 1 for the number of
    # 'atoms' per box required to add a sphere

    if old_axes:
        wet_spheres, x_axis, y_axis, z_axis = create_sphere_model(wet_spheres,
                                                                  1,
                                                                  box_side,
                                                                  xaxis = x_axis,
                                                                  yaxis = y_axis,
                                                                  zaxis = z_axis)

    else:
        wet_spheres, x_axis, y_axis, z_axis = create_sphere_model(wet_spheres,
                                                                  1,
                                                                  box_side)

    return wet_spheres

def residual2_cut(cutoff, box_side, dry_spheres, hydration_no, targ_volume):
    """
    Compute squared residual of sphere model to the theroetical target volume

    @type  cutoff:        integer
    @param cutoff:        Cutoff used to remove overlapping and excess test
                          hydration spheres
    @type  box_side:      float
    @param box_side:      Length of the sides of the cubes used to create a grid
                          to coarse grain atomic structure into a sphere model
                          (and place hydration spheres).
    @type  dry_spheres:   list
    @param dry_spheres:   A list containing lists of x, y & z coordinates
                          (3 * floats) for all spheres in an non-hydrated sphere
                          model.
    @type  hydration_no:  integer
    @param hydration_no:  Number defining which positions on a cube surrounding
                          each input sphere should be used to place test
                          hydration spheres. Position 1 = original sphere
                          position, 2 to 27 positions on cube centred on
                          original sphere.
    @type  targ_volume:   float
    @param targ_volume:   The target volume for the hydrated model calculated
                          from the protein sequence.
    @rtype:               float
    @return:              Squared difference between hydrated model volume and
                          the target volume.
    """

    # Create sphere model
    wet_model = hydrate_sphere_model(dry_spheres, hydration_no, box_side, cutoff)

    # Compute volume of the model
    vol_spheres = len(wet_model) * box_side**3

    # Return squared residual
    return (vol_spheres - targ_volume)**2

def optimize_watercut(box_side, coords, hydration_no, targ_vol, cut_min, cut_max, tolerance):
    """
    Get optimal cutoff to use to filter hydration spheres to reproduce
    target_volume (and deviation)

    @type  box_side:      float
    @param box_side:      Length of the sides of the cubes used to create a grid
                          to coarse grain atomic structure into a sphere model
                          (and place hydration spheres).
    @type  coords:        list
    @param coords:        A list containing lists of x, y & z coordinates
                          (3 * floats) for all spheres in an non-hydrated sphere
                          model.
    @type  hydration_no:  integer
    @param hydration_no:  Number defining which positions on a cube surrounding
                          each input sphere should be used to place test
                          hydration spheres. Position 1 = original sphere
                          position, 2 to 27 positions on cube centred on
                          original sphere.
    @type  targ_vol:      float
    @param targ_vol:      The target volume for the hydrated model calculated
                          from the protein sequence.
    @rtype:               integer, float
    @return:         
                          - Optimized cutoff value
                          - Deviation of the calculated volume from the target
                          volume

    """

    # Format the bounds within which to run the optimization
    cut_bounds = (cut_min, cut_max)

    # Minimize the squared residuals between the sphere model and target volumes
    # Uses minimize_scalar from scipy.optimize
    opt = optimize.minimize_scalar(residual2_cut, args=(box_side, coords, hydration_no, targ_vol),
                                   bounds=cut_bounds, method='bounded',
                                   options={'xtol' : tolerance})

    # Return the optimized cutoff and residual
    return int(round(opt['x'])), opt['fun']**0.5

def sphere_squared_form_factor(q, r):
    """
    Return the squared form factor for the magnitude of scattering vector, q,
    and sphere with given radius, r.

    @type  q:  float
    @param q:  Magnitude of scattering vector
    @type  r:  float
    @param r:  Sphere radius
    @rtype:    float
    @return:   Squared form factor
    """

    # Form factor^2 = (3 * [sin(q*r) - (q*r)*cos(q*r)])^2 / (q*r)^6

    qr = q * r
    qr6 = qr**6

    numerator = (3 * (np.sin(qr) - qr * np.cos(qr)))**2

    return numerator / qr6


def calc_r_hist(coords, no_bins):
    """
    Calculate histrogram of pair distances between the coordinates using
    no_bins bins. Return hist (counts), bin_edges and last non-zero bin no.

    @type  coords:   list
    @param coords:   A list containing lists of x, y & z coordinates
                     (3 * floats) of spheres
    @type  no_bins:  integer
    @param no_bins:  Number of bins to use when creating histogram of pair
                     distances between spheres.
    @rtype:          numpy array, numpy array, integer
    @return:   
                     1. Bins containing the counts of atoms in each division of
                     the histogram of pair distances between spheres.

                     2. The values at the edges of the bins in the histogram of
                     pair distances between spheres (corresponds to the counts
                     returned).

                     3. Number of the last bin with a non-zero count.
    """

    pair_dists = dist.pdist(coords, 'euclidean')
    hist, bin_edges = np.histogram(pair_dists, bins = no_bins)
    last_used = len([i for i in hist if i > 0])

    return hist, bin_edges, last_used

def spheres_to_sas_curve(coords, radius, q_max, n_points, **kwargs):
    """
    Calculate theoretical curve given sphere coordinates and radius using the
    Debye equation. Assumes all spheres are identical. A curve containing
    n_points is computed over a q range of 0 to q_max.

    @type  coords:    list
    @param coords:    A list containing lists of x, y & z coordinates
                      (3 * floats) of spheres.
    @type  radius:    float
    @param radius:    Sphere radius
    @type  q_max:     float
    @param q_max:     Maximum magnitude of the scattering vector, q, for which
                      to calculate the theoretical scattering intensity.
    @type  n_points:  integer
    @param n_points:  Number of points to calculate in the theoretical
                      scattering curve.
    @keyword rbins:   Number of bins to use when creating histogram of pair
                      distances between spheres (for use in the Debye equation).
                      Default = 400.
    @rtype:           numpy array
    @return:          Two dimensional array containing scattering vector
                      magnitude, q, and intensity, I, value pairs forming a
                      theoretical scattering curve.
    """

    # Number of bins to use when creating histogram of pair distances between
    # spheres
    n_bins = kwargs.get('rbins', 400)

    # Create set of q points for which the scattering curve is calculated
    q_delta = q_max / n_points

    q = []
    for n in range(1, n_points + 1):
        q.append(n * q_delta)

    natom = len(coords)

    # Debye equation:
    # I/Io = g * [n^-1 = 2 * n^-2 * sum(A_j * sin(q * r_j)/ q * r_j)]
    # Where g = squared form factor
    # n = number of spheres
    # r_j = the jth pair distance between spheres from a histogram of them (the
    # sum is over all such distances)
    # A_j = the number of distances r_j from the histogram

    # Create histogram of pair distances between spheres
    hist, bin_edges, last_used = calc_r_hist(coords, n_bins)

    # Create the sum of in teh Debye equation
    sigma = []

    for jj in range(0, n_points):
        sigma_tmp = 0.0
        for kk in range(0, n_bins):
            qr = q[jj] * bin_edges[kk + 1]
            sin_term = np.sin(qr)/ qr
            sigma_tmp += hist[kk] * sin_term
        sigma.append(sigma_tmp)

    inv_n = 1.0 / natom
    inv_n2 = 1.0 / natom**2

    scat = []

    for jj in range(0, n_points):
        scat.append([q[jj], sphere_squared_form_factor(q[jj], radius)
                    * (inv_n + 2.0 * inv_n2 * sigma[jj])])

    return np.array(scat)
