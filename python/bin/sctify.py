#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Sctify: converts a CHARMM PSF/PDF pair to a SCT compatible PDB
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

import argparse
import sct

def parse_arguments():
    """Parse command line arguments and ensure correct combinations present"""

    parser = argparse.ArgumentParser(
        description= 'Convert a CHARMM PSF/PDF pair to a SCT compatible PDB\n')

    parser.add_argument('-i','--input_pdb', nargs='?', type=str,
        dest='pdb_path', help = 'Path to the input PDB file',
        required=True)
    parser.add_argument('-p','--input_psf', nargs='?', type=str,
        dest='psf_path', help = 'Path to the input PSF file',
        required=True)
    parser.add_argument('-o','--output_pdb', nargs='?', type=str,
        dest='pdb_out', default=None,
        help = 'Path to the output PDB file')

    return parser.parse_args()
    
def main():

    args = parse_arguments()

    atoms = sct.pdb.process_pdb_psf(args.psf_path, args.pdb_path)
    sct.pdb.write_pdb(atoms, args.pdb_out)

if __name__ == "__main__":
    main()

