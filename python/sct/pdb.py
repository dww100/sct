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
import string
import sys

# Dictionary of CHARMM residue names to standard naming
charmm_resids = {
'ANE5AC': 'SIA',
'BNE5AC': 'SIA',
'AFUC': 'FUC',
'BFUC': 'FUC',
'AGAL': 'GAL',
'BGAL': 'GAL',
'AGLCNA': 'NAG',
'BGLCNA': 'NAG',
'AGALNA': 'NGA',
'BGALNA': 'NGA',
'AMAN': 'MAN',
'BMAN': 'MAN',
'ASH': 'ASP',
'CYM': 'CYS',
'CYX': 'CYS',
'GLH': 'GLU',
'HIP': 'HIS',
'HID': 'HIS',
'HIE': 'HIS',
'HSD': 'HIS',
'HSE': 'HIS',
'LYN': 'LYS',
'TYM': 'TYR'
}

charmm_residues = charmm_resids.keys()

# In the PDBs the CHARMM resids get truncated to three characters.
charmm_three_char = ['ANE','BNE', 'AFU', 'BFU', 'AGA', 'BGA', 'AGL', 'BGL', 'AMA', 'BMA']

PDB_ATOM_RECORD = (
    ('record', 0,  6, None),    # 0,  record type name
    ('atom_no', 6, 11, int),     # 1,  atom serial/index no
    ('atom_name', 12, 16, None),     # 2,  atom name
    ('alt', 16, None, None),    # 3,  altLoc
    ('res_id', 17, 20, None),  # 4,  residue name
    ('chain', 21, None, None),   # 5,  chain ID
    ('res_no', 22, 26, int),     # 6,  residue sequence number
    ('insert', 26, None, None), # 7,  insertion code
    ('x', 30, 38, float),       # 8,  x
    ('y', 38, 46, float),       # 9,  y
    ('z', 46, 54, float),       # 10, z
    ('occ', 54, 60, float),     # 11, occupancy
    ('b', 60, 66, float),       # 12, temperature factor
    ('segid', 72, 76, None),   # 13, segment ID
    ('element', 76, 78, None),  # 14, element
    ('charge', 78, 80, None),   # 15, charge
)

#  CHARMM Format strings - source/psffres.src
#
#  Format shown including CHEQ (which we are ignoring):
#  II,LSEGID,LRESID,LRES,TYPE(I),IAC(I),CG(I),AMASS(I),IMOVE(I),ECH(I),EHA(I)
#
# standard format:
# (I8,1X,A4,1X,A4,1X,A4,1X,A4,1X,I4,1X,2G14.6,I8,2G14.6)
# (I8,1X,A4,1X,A4,1X,A4,1X,A4,1X,A4,1X,2G14.6,I8,2G14.6)  XPLOR
# expanded format EXT:
# (I10,1X,A8,1X,A8,1X,A8,1X,A8,1X,I4,1X,2G14.6,I8,2G14.6)
# (I10,1X,A8,1X,A8,1X,A8,1X,A8,1X,A4,1X,2G14.6,I8,2G14.6) XPLOR
#
# Note: For some reason the extended format seems to like putting small
# charges in scientific format which leads to the length increasing by 2 chars
# Works on files I tested but YMMV

PSF_ATOM_RECORD = (
    ('atom_no', 0,  8, int),     # 0, atom serial/index no
    ('segid', 9, 13, None),     # 1, segment ID
    ('res_no', 14, 18, int),     # 2, residue sequence number
    ('res_id', 19, 23, None),  # 3, residue name
    ('atom_name', 24, 28, None),     # 4, atom name
    ('type', 29, 33, None),     # 5, atom type
    ('charge', 34, 48, float),  # 6, charge
    ('mass', 48, 62, float),    # 7, mass
    ('imove', 62, 70, int),     # 8, imove
#            ('ech', 70, 84, float),     # 9, ech
#            ('eha', 84, 100, float),    # 10, eha    
)

PSFEXT_ATOM_RECORD = (
    ('atom_no', 0,  10, int),    # 0, atom serial/index no
    ('segid', 11, 19, None),    # 1, segment ID
    ('res_no', 20, 28, int),     # 2, residue sequence number
    ('res_id', 29, 37, None),  # 3, residue name
    ('atom_name', 38, 46, None),     # 4, atom name
    ('type', 47, 52, None),     # 5, atom type
    ('charge', 53, 68, float),  # 6, charge
    ('mass', 68, 82, float),    # 7, mass
    ('imove', 82, 90, int),     # 8, imove
#            ('ech', 90, 104, float),     # 9, ech
#            ('eha', 104, 118, float),    # 10, eha  
)

# We will accept either standard or CHARMM residue naming
accept_resids = seq.all_residues + charmm_residues + charmm_three_char

def parse_line(line, schema, filename):
    """
    Parse an atom record line from a PDB or PSF file via schema.
    Based on the PDB parsing code written by Alisue (lambdalisue@hashnote.net):
    https://gist.github.com/lambdalisue/7209305
    
    Schema format:

        1) record name
        2) start index
        3) end index or None for single character
        4) processing function or None for `string.strip`

        so a truncated (but valid) schema for a PDB ATOM record cold be:

        ATOM_RECORD_SCHEMA = (
        ...    ('record', 0,  6,   None),   # 0,  record name
        ...    ('name', 6, 11,   int),      # 1,  serial
        )

    @type  line:   string
    @param line:   Line from a text file.
    @type  schema: list
    @param schema: List of 4-tuples. The format of the tuples is given above.
    @rtype:        dictionary
    @return:       Dictionary with keys specifiec as the record name in the 
                   schema.
    """
    
    vals = {}
    for record, start, end, fn in schema:
        if end is None:
            end = start + 1
        if fn is None:
            fn = string.strip
        try:
            vals[record] = fn(line[start:end])
        except ValueError, e:
            print "Error: File %s contains a line which could not be parsed according to the schema." % filename
            print "Python error message: %s" % e
            print "Line contents:", line
            raise IOError("Error parsing %s, ignoring this pdb file." % filename)
    return vals

def pdb_res_line_parse(line, filename):
    """Parse a single line from a PDB file into a dictionary according to
    the standard PDB column definitions

    @type  line: string
    @param line: Line from a PDB file.
    @rtype:      dictionary
    @return:     Contains:
                     - atom_no:   integer, atom number
                     - atom_name: string, atom name
                     - res_id:    string, residue three letter code
                     - chain:     string, chain identifier
                     - res_no:    integer, residue number
                     - coords:    list of floats, the x, y & z coordinates
    """

    # Only interested in ATOM and HETATM lines not REMARKs etc.
    # The type of record on the line is determined by the first 6 characters
    record = line[0:6].strip()
    
    if record in ['ATOM','HETATM']:

        data = parse_line(line, PDB_ATOM_RECORD, filename)
        if data['res_id'] in accept_resids:
        # Split data on the line according to column definitions for PDB files
        # Ignore residues that we can't handle in SCT
        # TODO: This should perhaps give a warning

            data['coords'] = [data['x'], data['y'], data['z']]
        else:
            raise IOError("Error: File %s contains unknown residue %s" % (filename, data['res_id']))
            data = {}
    else:
        data = {}
            
    return data

def read_pdb_atom_data (filename):
    """
    Read PDB file and return residue frequencies and atom coordinates

    @type  filename: string
    @param filename: Path to a PDB file.
    @rtype:   dictionary, list
    @return:  
              - Dictionary for residue type frequency. 
                Three letter residue code for the key and residue type
                frequency as the values. 
              - A list containing lists of x, y & z coordinates (3 * floats)
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
            data = pdb_res_line_parse(line, filename)

            if len(data) != 0:
                atom_coords.append(data['coords'])
                # If residue number has changed increment res_id tally

                if data['res_no'] != last_res_no:
                    # If a CHARMM residue name is used convert to standard naming
                    if data['res_id'] in charmm_resids:
                        res_name = charmm_resids[data['res_id']]
                    else:
                        res_name = data['res_id']
                    res_freq[res_name] += 1
                    last_res_no = data['res_no']

    return res_freq, atom_coords

def read_sphere_pdb(filename):
    """
    Read PDB file and return residue frequencies and atom coordinates

    @type  filename: string
    @param filename: Path to a PDB file.
    @rtype:          dictionary, list
    @return:   
              - A list containing lists of x, y & z coordinates (3 * floats)
              - float of the sphere radius
    """
    
    spheres = []
    radius = 0.0
    # Read in lines of the PDB
    # ATOM and HETATM records are interpretted - others ignored
    with open(filename) as f:
        for line in f:
            data = pdb_res_line_parse(line, filename)

            if len(data) != 0:
                spheres.append(data['coords'])
                radius = data['b']       
    
    return spheres, radius

def read_psf(filename):
    """
    Read in a PSF file and create a list of dictionaries containing the fields 
    defined in PSF_ATOM_RECORD/PSFEXT_ATOM_RECORD schemas.

    @type  filename: string
    @param filename: Path to a PSF file.
    @rtype:          list
    @return:         List of dictionaries with kets as defined in the
                     PSF_ATOM_RECORD/PSFEXT_ATOM_RECORD schemas.
    """    
    
    
    psf_file = open(filename,'r')
    # Process PSF header
    # Contains the flags that determine contents and format type 
    # (EXT or normal, not XPLOR or CHARMM)
    header = psf_file.readline().split()
    if header[0] != 'PSF':
        # This should be a proper exception! As should the other errors!
        print 'Boo not a PSF!'
        sys.exit()
        
    if 'EXT' in header:
        schema = PSFEXT_ATOM_RECORD
    else:
        schema = PSF_ATOM_RECORD
    # Skip blank line
    psf_file.readline()

    title = psf_file.readline().split()
    if title[1] != '!NTITLE':
        print 'Not a valid PSF (NTITLE) - sad times'
        sys.exit()

    # Skip remark lines (number given as the first entry on the NTITLE line)                
    for ii in range(int(title[0])):
        psf_file.readline()
        
    # Skip blank line
    psf_file.readline()

    atom_title = psf_file.readline().split()
    if atom_title[1] != '!NATOM':
        print 'Not a valid PSF (NATOM) - sad times'
        sys.exit()
    
    n_atoms = int(atom_title[0])
    
    atoms = []
    
    for ii in xrange(n_atoms):
        line = psf_file.readline()
        atoms.append(parse_line(line, schema, filename))

    psf_file.close()
    
    return atoms
    
def process_pdb_psf(psf_filename, pdb_filename):
    """
    Read in a PSF/PDB pair of files and create a list of dictionaries 
    containing the fields defined in PSF_ATOM_RECORD/PSFEXT_ATOM_RECORD 
    schemas plus 'coords' (from the PDB x,y,z coordinates) and 'chain' (from 
    the PDB chain column). Hydrogen atoms are filtered out.

    @type  filename: string
    @param filename: Path to a PSF file.
    @rtype:          list
    @return:         List of dictionaries with kets as defined in the
                     PSF_ATOM_RECORD/PSFEXT_ATOM_RECORD schemas plus 'coords' 
                     (from the PDB x,y,z coordinates) and 'chain' (from the PDB
                     chain column)..
    """
    
    psf_atoms = read_psf(psf_filename)
    n_atoms = len(psf_atoms)        
    
    if n_atoms == 0:
        print "No atoms foud in the PSF, expect pain"
    else:
        tmp_atoms = []
    
    # Read in lines of the PDB
    # ATOM and HETATM records are interpretted - others ignored
    with open(pdb_filename) as f:
        for line in f:
            if line[:6].strip() in ['ATOM', 'HETATM']:
                data = pdb_res_line_parse(line, pdb_filename)
                tmp_atoms.append(data)        
        
        if len(tmp_atoms) != n_atoms:
            print "Number of atoms in the PDB do not match the PSF"
        else:
            for ii in range(0, n_atoms - 1):
                #print psf_atoms[ii]['atom_name'] + '\t' + tmp_atoms[ii]['atom_name']
                psf_atoms[ii]['coords'] = tmp_atoms[ii]['coords']
                psf_atoms[ii]['chain'] = tmp_atoms[ii]['chain']
    
    
    atoms = []    
    
    for atom in psf_atoms:
        if atom['atom_name'][0] != 'H':
    
            if atom['res_id'] in charmm_residues:
                atom['res_id'] = charmm_resids[atom['res_id']]
                
            atoms.append(atom)                            
                
    return atoms

def create_pdb_atom(res_no, res_id, atom_no, atom_type, coords, **kwargs):
    """
    Create ATOM line for output in a PDB file (ends with a new line character)
    
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

def write_pdb(atoms, filename):
    """
    Writes out a PDB file based on the values in the input atoms

    @type  atoms:     list
    @param atoms:     A list of dictionaries containing values for the keys 
                      defined in PDB_ATOM_RECORD schema
    @type  filename:  string
    @param filename:  Filename for the output sphere model PDB
    """
    
    out_file = open(filename, 'w')
    
    for atom in atoms:
        pdb_line = create_pdb_atom(atom['res_no'], 
                                   atom['res_id'], 
                                   atom['atom_no'], 
                                   atom['atom_name'], 
                                   atom['coords'],
                                   chain = atom['chain'])
    
        out_file.write(pdb_line)

    out_file.close()

def write_sphere_pdb(coords, radius, filename):
    
    """
    Writes out a PDB file containing entries for the spheres of a sphere model,
    using the beta column to store radius (all ATOMs are C1 SER)

    @type  coords:    list
    @param coords:    A list containing lists of x, y & z coordinates (3 * floats)
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
