#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
PDB related functions used in the SCT suite of programs
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

import seq

def pdb_res_line_parse(line):
    """Parse a single line from a PDB file into a dictionary according to
    the standard PDB column definitions

    @type  line: string
    @param line: Line from a PDB file.
    @rtype:      dictionary
    @return:     A dictionary containing:
                 - atom_no:   integer, atom number
                 - atom_name: string, atom name
                 - res_id:    string, residue three letter code
                 - chain:     string, chain identifier
                 - res_no:    integer, residue number
                 - coords:    list of floats, the x, y & z coordinates
    """

    data = {}

    # Only interested in ATOM and HETATM lines not REMARKs etc.
    # The type of record on the line is determined by the first 6 characters
    record = line[0:6].strip()
    if record in ['ATOM','HETATM']:

        # Split data on the line according to column definitions for PDB files
        # Ignore residues that we can't handle in SCT
        # TODO: This should perhaps give a warning

        # res_id is a three letter residue code
        res_id = line[17:20].strip()
        if res_id in seq.all_residues:
            data['atom_no'] = int(line[6:11])
            data['atom_name'] = line[12:16].strip()
            data['res_id'] = res_id
            data['chain'] = line[21:22].strip()
            data['res_no'] = int(line[22:26])
            x = float(line[30:38])
            y = float(line[38:46])
            z = float(line[46:54])
            data['coords'] = [x, y, z]

    return data

def read_pdb_atom_data (filename):
    """
    Read PDB file and return residue ferquencies and atom coordinates

    @type  filename: string
    @param filename: Path to a PDB file.
    @rtype:   dictionary, list
    @return:  1. Dictionary for residue type frequency. Three letter residue
              code for the key and residue type frequency as the values.

              2. A list containing lists of x, y & z coordinates (3 * floats)
    """

    # Initialize a dictionary of all recognized residues (frequency = 0)
    res_freq = seq.residue_freq_dict()

    # Track the last residues seen
    last_res_no = 0

    atom_coords = []

    # Read in lines of the PDB
    # ATOM and HETATM records are interpretted - others ignored
    with open(filename) as f:
        for line in f:
            data = pdb_res_line_parse(line)

            if len(data) != 0:
                atom_coords.append(data['coords'])
                # If residue number has changed increment res_id tally
                if data['res_no'] != last_res_no:
                    res_freq[data['res_id']] += 1
                    last_res_no = data['res_no']

    return res_freq, atom_coords

def create_pdb_atom(res_no, res_id, atom_no, atom_type, coords, **kwargs):
    """
    Create ATOM line for output in a PDB file (output has a \n)

    @type  res_no:    integer
    @param res_no:    Residue number
    @type  res_id:    string
    @param res_id:    Residue name (3 letter code)
    @type  atom_no:   integer
    @param atom_no:   Atom number
    @type  atom_type: string
    @param atom_type: Atom type
    @type  coords:    list
    @param coords:    A list containing lists of x, y & z coordinates
                      (3 * floats)
    @keyword beta:    Value to put in the temperature factor (beta) column
    @keyword occ:     Value to put in the occupancy (occ) column
    @keyword chain:   Chain identifier (a single letter)
    @rtype:           string
    @return:          A line of text formatted as an ATOM line of a PDB
    """

    # Columns for temperature factor (beta) and occupancy (occ) often used to
    # store other data
    beta = kwargs.get('beta', 0.0)
    occ = kwargs.get('occ', 0.0)
    chain = kwargs.get('chain', 'A')


    line = "ATOM  {0:5d} {1:4s} {2:3s}".format(atom_no, atom_type, res_id)
    line += " {0:1s}{1:4d}    ".format(chain, res_no)
    line += "{0:8.3f}{1:8.3f}{2:8.3f}".format(coords[0],coords[1],coords[2])
    line += "{0:6.2f}{1:6.2f}\n".format(occ, beta)

    return line

def write_sphere_pdb(coords, radius, filename):
    """
    Writes out a PDb file containing entries for the spheres of a sphere model,
    using the beta column to store radius (all ATOMs are C1 SER)

    @type  coords:    list
    @param coords:    A list containing lists of x, y & z coordinates
                      (3 * floats)
    @type  radius:    float
    @param radius:    Sphere radius
    @type  filename:  string
    @param filename:  Filename for the output sphere model PDB
    """

    out_file = open(filename, 'w')

    count = 1

    for coord in coords:

        pdb_line = create_pdb_atom(count, 'SER', count, 'C1', coord, beta = radius)

        out_file.write(pdb_line)
        count += 1

    out_file.close()
