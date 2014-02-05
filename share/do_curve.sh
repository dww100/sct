#!/bin/bash
#
# do_curve
#
# Control script for calculating scattering curves for pdb models
# Initiates; BRKTOS97d, HYPRO, APS, SCT, SCTPL6g, RFAC700XN
# Versions of codes with _XXXX are edited to handle XXXX atom coordinates

#-------------------------------------------------------------------------------
# Copyright 1981-2014 University College London
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#    http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#
#-------------------------------------------------------------------------------




# Unless you know what you are doing copy this script to the directory
# containing your initial models as PDBs and edit the three variables below.
# To run: ./do_curve.sh

################################################################################
# Variables to be edited for each run
#
# PREFIX = the prefix for all initial model files (they should be named
# PREFIX1.pdb, PREFIX2.pdb etc.)
PREFIX='test'
# COUNT is the variable which increments to run through the numbered PDBs
# Initialize to the first numbered PDB you wish to analyse
COUNT=1
# LAST_NO is the number of teh last PDB to analyse
LAST_NO=1
################################################################################

mkdir {HYP,BRK,SCX,SCN,SUM}

while [ $COUNT -le $LAST_NO ]
do

    PREFIX_P=$PREFIX$COUNT
    PDB=$PREFIX_P'.pdb'
    BRK=$PREFIX_P'.brk'
    OUR_PATH='/usr0/usr3/people/davidw/sct'

#
################################################################################
#
##  Dry Sphere Modelling:
##    1. BRKTOS97b - note that now we use a more update version called
##  BRKTOS97d
##  To create dry sphere models of desired volume, I have found that it is
##  necessary to start from several initial models.
##  Position models at different positions in 'coordinate space'.  Experiment
##  with cube-side and cutoff values to produce 3-dimensional grid that
##  satisfies all of the models. Hopefully this will be robust enough so that a
##  decrease in dry model volume is due to steric overlap within the pdb model
##  and not because of its orientation within the grid!
#
#


    echo -e $PDB'\n5.57 4 1 0\n'$BRK'\n500 500 500\n-500 -500 -500\n' > tmp_in
    $OUR_PATH/brktos97d < tmp_in
    wc -l $BRK >> $PREFIX'dry_spheres.sum'

#
################################################################################
#
##  Wet Sphere Modelling:
##    1. HYPRO - add high number of spheres to dry model
##    2. BRKTOS97b - use high cutoff (>9) and reset grid to get wet sphere
##    model (1) of desired volume
##    3. HYPRO - repeat procedure for dry model to recreate HYPRO model but
##    without adding spheres (sphere=1)
##    4. BRKTOS97b - repeat with cutoff = 1 for unhyrated HYPRO model (3) to
##    align in same grid as hydrated model (2)
##    5. Cat sphere models together (2+4), and use sort and uniq commands to
##    remove duplicate sphere
##    This procedure produces a hydrated model without losing extended peptide
##    or carbohydrate structures which can happen if do not add initial
##    unhydrated sphere model back to hydrated model, however this means that
##    hydration is localised to the compact regions of the structure.
##    To get the desired volume for the hydrated models (5) it is necessary to
##    determine empirically the number of spheres used in (1) and the cutoff
##    used in (2). Note that we now use BRKTOS97d
#

    HYP=$PREFIX_P'.hyp'
    echo -e $BRK'\n'$HYP'1\n26\n' >tmp_in
    $OUR_PATH/hypro < tmp_in
    echo -e $HYP'1\n5.57 1 1 0\n'$HYP'2\n500 500 500\n-501 -501 -501\n' > tmp_in
    $OUR_PATH/brktos97d < tmp_in
    echo -e $BRK'\n'$HYP'3\n1\n' > tmp_in
    $OUR_PATH/hypro < tmp_in
    echo -e $HYP'3\n5.57 1 1 0\n'$HYP'4\n500 500 500\n-501 -501 -501\n' > tmp_in
    $OUR_PATH/brktos97d < tmp_in
    cat $HYP'2' > $HYP'5temp'
    cat $HYP'4' >> $HYP'5temp'
    sort -n $HYP'5temp' | uniq > $HYP'5'
    wc -l $HYP'5' >> $PREFIX'wet_spheres.sum'

#
################################################################################
#
#
##  APS Rg Calculation
##    1. Calculate Rg of dry sphere model
##    2. Calculate Rg of wet sphere model
#

    $OUR_PATH/aps < $BRK > tmp_out
    OUT_RG=`fgrep -i "radius of gyration   =" tmp_out`
    echo -e $BRK'  '$OUT_RG >> $PREFIX'Rg_dry_aps.sum'
    $OUR_PATH/aps < $HYP'5' > tmp_out
    OUT_RG=`fgrep -i "radius of gyration   =" tmp_out`
    echo -e $HYP'5 '$OUT_RG >> $PREFIX'Rg_wet_aps.sum'

################################################################################
#
##  SCT calculation of scattering curves
##    1. Calculate neutron curve from dry model, using smearing values
##    of 0.1 6 0.016 for RAL
##    2. Calculate X-ray curve from wet model, unsmeared

    N=$PREFIX_P'.scn'
    X=$PREFIX_P'.scx'
    echo -e $BRK'\n'$N'\nx\n\n0\n1 1\n0.1 6 0.016\n0.92\n0.25 2.78\n' > tmp_in
    $OUR_PATH/sct_6000 < tmp_in
    echo -e `grep -i 'VALUE OF ' x`'   '$N >> $PREFIX'NEUT_sct.sum'
    echo -e `grep -i ' TOTAL OF CROSS-TERMS OF MIXED SPHERES IS' x`'   '$N >> $PREFIX'NEUT_sct.sum'
    echo -e '------------------------------------------------------------------' >> $PREFIX'NEUT_sct.sum'
    echo -e $HYP'5\n'$X'\nx\n\n0\n1 0\n0.92\n0.25 3.47\n' > tmp_in
    $OUR_PATH/sct_6000 < tmp_in
    echo -e `grep -i 'VALUE OF ' x`'   '$X >> $PREFIX'XRAY_sct.sum'
    echo -e `grep -i ' TOTAL OF CROSS-TERMS OF MIXED SPHERES IS' x`'   '$X >> $PREFIX'XRAY_sct.sum'
    echo -e '------------------------------------------------------------------' >> $PREFIX'XRAY_sct.sum'
    cp $N neutron
    cp $X XRAY

################################################################################
##  SCTPL6 Calculation of Rg
##    1. Calculate Rg of neutron curve
##    2. Calculate Rg of X-ray curve


    echo -e 'm\n'$N' '$PREFIX_P'\n1\nneutron\nn\n2\n0 0.07\n0.016 0.05\n\n12\n' > tmp_in
    $OUR_PATH/sctpl6g < tmp_in
    cat Sctpl6_Summary >> $PREFIX'NEUTsctpl_rg.sum'
    rm Sctpl6_Summary

    echo -e 'm\n'$X' '$PREFIX_P'\n1\nXRAY\nn\n2\n0 0.07\n0.016 0.05\n\n12\n' >tmp_in
    $OUR_PATH/sctpl6g < tmp_in
    cat Sctpl6_Summary >> $PREFIX'XRAYsctpl_rg.sum'
    rm Sctpl6_Summary

##  SCTPL6 Calculation of Rxs
##    1. Calculate Rxs of neutron curve
##    2. Calculate Rxs of X-ray curve

    echo -e 'm\n'$N' '$PREFIX_P'\n1\nneutron\nn\n1\n0.01 0.15\n0.055 0.1\n\n12\n' > tmp_in
    $OUR_PATH/sctpl6g < tmp_in
    cat Sctpl6_Summary >> $PREFIX'NEUTsctpl_rxs.sum'
    rm Sctpl6_Summary

    echo -e 'm\n'$X' '$PREFIX_P'\n1\nXRAY\nn\n1\n0.01 0.15\n0.055 0.1\n\n12\n' > tmp_in

    $OUR_PATH/sctpl6g < tmp_in
    cat Sctpl6_Summary >> $PREFIX'XRAYsctpl_rxs.sum'
    rm Sctpl6_Summary

################################################################################
#
##  RFACXN Calculation of R-factor
##    1. Calculate R-factor for neutron4 curve against r___.q, Q xx --> 0.2Ang
##    2. Calculate R-factor for X-ray curve against e__.3, Q xx --> 0.2Ang


    echo -e $N'\nneut.dat\nRfac.out\n0.2\n' > tmp_in
    $OUR_PATH/rfacXN < tmp_in
    cat Rfac.out >> $PREFIX'Neutron_rfac.sum'
    rm Rfac.out


    #echo -e $X'\nFH_68_100cAV\nRfac.out\n0.2\n' > tmp_in
    #$OUR_PATH/rfac700XN < tmp_in
    #cat Rfac.out >> $PREFIX'Xray_rfac.sum'
    #rm Rfac.out


##    Clean Up Files

    mv $HYP'5' HYP/
    mv $X SCX/
    mv $N SCN/
    mv $BRK BRK/

    rm $HYP'1'
    rm $HYP'2'
    rm $HYP'3'
    rm $HYP'4'
    rm $HYP'5temp'
    rm neutron
    rm XRAY

    COUNT=$[$COUNT+1]

done

mv *sum SUM/
