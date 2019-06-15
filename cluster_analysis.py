#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 12 14:30:58 2018

@author: Maxim Makeev
"""
#
import os
import sys
import math
#
#
#---Function checks the foces on molecules in identified clusters-
def force_on_mol( ll1, ll2, fxx, fyy, fzz ):
    s_fx = float( 0.0 )
    s_fy = float( 0.0 )
    s_fz = float( 0.0 )
    for inn in range( 0, ll2-ll1 ):
        s_fx = s_fx + fxx[inn]
        s_fy = s_fy + fyy[inn] 
        s_fz = s_fz + fzz[inn]
    f_const = float(1.0/16.0)*float( 0.043363 )
    res_f = min( s_fx*f_const, s_fy*f_const, s_fz*f_const ) 
    return res_f
#
#
#
def main():
    n_mol_types = []
    n_atom_mol  = []
#---------------------------------------------------------------
    CON_2_ATTY = []
    CON_2_ELEM = []
#---------------------------------------------------------------
#Containers for atomic IDs, types and Cartesian coordinates-----
    ATID = []
    ATTY = []
    MOTY = []
#
    ATXX = []
    ATYY = []
    ATZZ = []
#
    FFXX = []
    FFYY = []
    FFZZ = []
#----------------------------------------------------------------
    currwd = os.getcwd( )
    print("Current directory: code and input file:\n", currwd)
    clusterinput = open( "input_cluster.dat","r" )
#
    AI = clusterinput.readlines()
#
    clusterinput.close()
#----------------------------------------------------------------
    tmp01 = AI[1].split()
    rcut = float( tmp01[0] )
    print( rcut )
    tmp02 = AI[3].split()
    ntypes = int( tmp02[0] )
    print( ntypes )
    tmp03 = AI[5].split()
    tmp04 = AI[6].split()
    for kk in range( 0,ntypes ):
        CON_2_ATTY.append( int( tmp03[kk] ) )
        CON_2_ELEM.append( str( tmp04[kk] ) )
    for i in range( 0,ntypes): 
        print( CON_2_ATTY[i], end=" " )
        print('')
    for i in range( 0,ntypes): 
        print( CON_2_ELEM[i], end=" " )
        print('')
    tmp05 = AI[8].split()
    ion_type = int( tmp05[0] )
    print( ion_type )
    tmp06 = AI[10].split()
    mol_types = int( tmp06[0] )
    print( mol_types )
    tmp07 = AI[12].split()
    for kl in range( 0, mol_types):
        n_mol_types.append( int( tmp07[kl] ) )
    for i in range( 0,mol_types): print( n_mol_types[i] )
    tmp08 = AI[14].split()
    for km in range( 0,mol_types ):
        n_atom_mol.append( int( tmp08[km] ) )
    for i in range( 0,mol_types): print( n_atom_mol[i], end=" " )
    print('')
    tmp09 = AI[16].split()
    force_max = float( tmp09[0] )
    print( force_max )
    tmp10 = AI[18].split()
    filename = str( tmp10[0] )
    print( filename )
    tmp11 = AI[20].split()
    workdir0 = str( tmp11[0] )
    tmp12 = AI[22].split()
    struct = tmp12[0]
    frame_num  = int( tmp12[1] )
    print( struct, frame_num)
    tmp12 = AI[24].split()
    mode = str( tmp12[0] ) 
    print("mode=:", mode)
#-----------------------------------------------------------------    
    os.chdir( workdir0 )
    pathdir = os.getcwd()
    print( pathdir )
    coorfile = open( filename,"r" )
    #
    A = coorfile.readlines()
    #
    coorfile.close()
#
#-Create arrays of molecular ranges------------------------------
#
    s_index = [0]*int( mol_types )
    e_index = [0]*int( mol_types ) 
    n_start = int( 0 )
    for im in range( 0,mol_types ):
        s_index[im] = n_start 
        e_index[im] = n_start + n_mol_types[im]*int( n_atom_mol[im] ) 
    n_start = e_index[im]
#   print( s_index[im],  e_index[im] )
#
    num_atoms_total = int( 0 )
    for na in range( 0,int(mol_types) ):
        print( n_mol_types[ na ] )
        num_atoms_total = num_atoms_total + n_mol_types[ na ]*\
        n_atom_mol[ na ]
    if( mode == 'verbose' ):
        print( "Total number of atoms:", num_atoms_total  )
#-----------------------------------------------------------------
#-----------------------------------------------------------------
#---Read a frame from the file------------------------------------
    start_line = (num_atoms_total + int(9) )*frame_num
    end_line   = (num_atoms_total + int(9) )*(frame_num + int(1) ) 
#--Read atomic ids, types, and coordinates from the first frame---
#--The code uses the LAMMPS dunp file structure:------------------
#--See LAMMPS _dump_ file structure: ITEM: ATOMS id...------------
#
    lco = [0]*9
    natom = 0
    key1 = 'ATOMS'
    for i in range( start_line,end_line ):
        tmp = A[i].split()
        j = i - start_line
        if j == 8 and len(tmp) > 1 and tmp[1] == key1:
            m = len( tmp )
            for n in range( 0,m ):
                if tmp[n] == 'id':
                    lco[0] = (n-2)
                if tmp[n] == 'mol':
                    lco[1] = (n-2)
                if tmp[n] == 'type':
                    lco[2] = (n-2)
                if tmp[n] == 'x':
                    lco[3] = (n-2)
                if tmp[n] == 'y':
                    lco[4] = (n-2)
                if tmp[n] == 'z':
                    lco[5] = (n-2)
                if tmp[n] == 'fx':
                    lco[6] = (n-2)
                if tmp[n] == 'fy':
                    lco[7] = (n-2)
                if tmp[n] == 'fz':
                    lco[8] = (n-2)
        if( j == 3): NA = int( tmp[0] )
        if( j == 3): print( 'Number of atoms in each frame: {0:d}'.format(NA) )
        if( j == 5): Lx = float( tmp[1] ) - float( tmp[0] )
        if( j == 6): Ly = float( tmp[1] ) - float( tmp[0] )
        if( j == 7): Lz = float( tmp[1] ) - float( tmp[0] )
        if( j >= 9 and natom < NA):
            natom += 1
            ATID.append( int(  tmp[lco[0]] ) )
            ATTY.append( int(  tmp[lco[2]] ) )
            MOTY.append( int(  tmp[lco[1]] ) )
#
            ATXX.append( float(tmp[lco[3]] ) )
            ATYY.append( float(tmp[lco[4]] ) )
            ATZZ.append( float(tmp[lco[5]] ) )
#
            FFXX.append( float(tmp[lco[6]] ) )
            FFYY.append( float(tmp[lco[7]] ) )
            FFZZ.append( float(tmp[lco[8]] ) )
#------------------------------------------------------------------
    if( NA != natom):
        print( "Major consistency check failed:" )
        print( "Configuration was not read correctly." )
        print( "NA=:", NA, "natoms=:", natom)
        sys.exit()
#----Ordering Data by the atomic IDs and creating a list of ions---
#
    ION_ID = [ ]
#------------------------------------------------------------------    
    ORTY = [int( 0   ) for ii in range( 0,natom ) ]
    ORMO = [int( 0   ) for ii in range( 0,natom ) ]
    ORXX = [float(0.0) for ii in range( 0,natom ) ]
    ORYY = [float(0.0) for ii in range( 0,natom ) ]
    ORZZ = [float(0.0) for ii in range( 0,natom ) ]
    ORFX = [float(0.0) for ii in range( 0,natom ) ]
    ORFY = [float(0.0) for ii in range( 0,natom ) ]
    ORFZ = [float(0.0) for ii in range( 0,natom ) ]
    for ik in range( 0,natom ):
        c_id = int( ATID[ik] )
        ORTY[c_id-1] = ATTY[ik]
        ORMO[c_id-1] = MOTY[ik]
        ORXX[c_id-1] = ATXX[ik]
        ORYY[c_id-1] = ATYY[ik]
        ORZZ[c_id-1] = ATZZ[ik]
        ORFX[c_id-1] = FFXX[ik]
        ORFY[c_id-1] = FFYY[ik]
        ORFZ[c_id-1] = FFZZ[ik]
    for ll in range( 0,natom ):
        if( ORTY[ll] == ion_type):
            ION_ID.append( ll )
#
    n_of_ion = len( ION_ID ) 
#
#-----------------------------------------------------------------------
#
    IDCLUS = [ [ ] ]
    TYCLUS = [ [ ] ]
    MOCLUS = [ [ ] ]
    XXCLUS = [ [ ] ]
    YYCLUS = [ [ ] ]
    ZZCLUS = [ [ ] ]
    FXCLUS = [ [ ] ]
    FYCLUS = [ [ ] ]
    FZCLUS = [ [ ] ]
#
#---Building clusters: find molecules arouns ions of type ion_type------ 
#---This part counts molecules around ions of type ion_type-------------
#
    num_of_clust = int( 0 )
    red_flag = 0
    for ii in range( 0,n_of_ion ):
        ion_curr = ION_ID[ii]
        IDCLUS.append( [ ] )   
        TYCLUS.append( [ ] ) 
        MOCLUS.append( [ ] )
        XXCLUS.append( [ ] )
        YYCLUS.append( [ ] ) 
        ZZCLUS.append( [ ] )
        FXCLUS.append( [ ] )
        FYCLUS.append( [ ] ) 
        FZCLUS.append( [ ] )
#-----------------------------------------------------------------------
        IDCLUS[ ii ].append( int(  ion_curr ) )   
        TYCLUS[ ii ].append( ORTY[ ion_curr ] )
        MOCLUS[ ii ].append( ORMO[ ion_curr ] )
        XXCLUS[ ii ].append( ORXX[ ion_curr ] )
        YYCLUS[ ii ].append( ORYY[ ion_curr ] ) 
        ZZCLUS[ ii ].append( ORZZ[ ion_curr ] )
        FXCLUS[ ii ].append( ORFX[ ion_curr ] )
        FYCLUS[ ii ].append( ORFY[ ion_curr ] ) 
        FZCLUS[ ii ].append( ORFZ[ ion_curr ] )
#    print( ii, IDCLUS[ii][0], TYCLUS[ii][0],  MOCLUS[ii][0] )
        if( mode == 'verbose' ): 
            print("Processing ion number:", ii )
        for kx in range( 0, mol_types):
            for ky in range( s_index[kx], e_index[kx], n_atom_mol[kx] ):
                ml1 = ky
                ml2 = ky + n_atom_mol[ kx ]
                red_flag = 0
                for mj in range( ml1,ml2 ):
                        dx = float( ORXX[mj] - ORXX[ion_curr] )
                        dy = float( ORYY[mj] - ORYY[ion_curr] )
                        dz = float( ORZZ[mj] - ORZZ[ion_curr] )
                        if( dx >  float(0.5)*Lx ): dx = dx - Lx
                        if( dy >  float(0.5)*Ly ): dy = dy - Ly
                        if( dz >  float(0.5)*Lz ): dz = dz - Lz
                        if( dx < -float(0.5)*Lx ): dx = dx + Lx
                        if( dy < -float(0.5)*Ly ): dy = dy + Ly
                        if( dz < -float(0.5)*Lz ): dz = dz + Lz
                        rsq = math.sqrt( dx*dx + dy*dy + dz*dz )
                        if red_flag == 1:
                            break
                        if( rsq < rcut and mj != ion_curr ): 
#include molecule in the cluster ion_curr
                            red_flag = 1
                            m_num = int ( ( mj - s_index[kx] )/( n_atom_mol[kx] ) )
                            l1 = s_index[kx] + m_num*n_atom_mol[kx]
                            l2 = l1 + n_atom_mol[kx]
                            tem_fx = []
                            tem_fy = []
                            tem_fz = []
                            for kf in range( l1, l2):
                                tem_fx.append( float( ORFX[kk] ) )
                                tem_fy.append( float( ORFY[kk] ) )
                                tem_fz.append( float( ORFZ[kk] ) )
                            m_ff = \
                            force_on_mol(l1, l2, tem_fx, tem_fy, tem_fz)
                            if( m_ff < force_max):
                                for kk in range( l1,l2 ):
                                    IDCLUS[ ii ].append (int(kk) )
                                    TYCLUS[ ii ].append( int( ORTY[ kk ] ) )
                                    MOCLUS[ ii ].append( int( ORMO[ kk ] ) )
                                    XXCLUS[ ii ].append( float( ORXX[kk] ) )
                                    YYCLUS[ ii ].append( float( ORYY[kk] ) )
                                    ZZCLUS[ ii ].append( float( ORZZ[kk] ) )
                                    FXCLUS[ ii ].append( float( ORFX[kk] ) )
                                    FYCLUS[ ii ].append( float( ORFY[kk] ) )
                                    FZCLUS[ ii ].append( float( ORFZ[kk] ) )
#-----------atom belongs to the molecule:------------------------------
#------Remove effects due to periodic boundary conditions--------------
    for ixx in range( 0,n_of_ion ):
        for iyy in range( 0,len( IDCLUS[ixx] )  ):
            dx = float( XXCLUS[ ixx ][iyy] - XXCLUS[ ixx ][0] )
            dy = float( YYCLUS[ ixx ][iyy] - YYCLUS[ ixx ][0] )
            dz = float( ZZCLUS[ ixx ][iyy] - ZZCLUS[ ixx ][0] )
            if( dx >  float(0.5)*Lx ): 
                XXCLUS[ixx][iyy] = XXCLUS[ixx][iyy] - Lx
            if( dy >  float(0.5)*Ly ): 
                YYCLUS[ixx][iyy] = YYCLUS[ixx][iyy] - Ly
            if( dz >  float(0.5)*Lz ): 
                ZZCLUS[ixx][iyy] = ZZCLUS[ixx][iyy] - Lz
            if( dx < -float(0.5)*Lx ): 
                XXCLUS[ixx][iyy] = XXCLUS[ixx][iyy] + Lx
            if( dy < -float(0.5)*Ly ): 
                YYCLUS[ixx][iyy] = YYCLUS[ixx][iyy] + Ly
            if( dz < -float(0.5)*Lz ): 
                ZZCLUS[ixx][iyy] = ZZCLUS[ixx][iyy] + Lz           
#
    for jj in range( 0,n_of_ion ):
        if ( len(IDCLUS[jj]) > 1  ): num_of_clust += 1
    if( mode == 'verbose' ):
        print( "Total number of clusters=:", num_of_clust )
#  
#----Output-*.xyz file------------------------------------------------
    for ifile in range( 0,n_of_ion ):
        if ifile < 10: 
            fn1 = 0
            fn2 = ifile
            filename_p = "Cluster_"+str(fn1)+str(fn2)+str(".xyz")
        else:
            fn11 = ifile
            filename_p = "Cluster_"+str(fn11)+str(".xyz")
        outf1 = open(filename_p,"w+")
#    
        number_atoms = int( len( IDCLUS[ifile]  ) )
        if( number_atoms > int(1) ):
            outf1 = open(filename_p,"w+")
            outf1.write("%d\n\n"% ( number_atoms ) )
            for i in range( 0, number_atoms ):
                n_type_w = TYCLUS[ifile][i]
                atom_type = CON_2_ELEM[ n_type_w - 1]
                outf1.write("%2s %15.10f %15.10f %15.10f\n"%( \
                    atom_type, XXCLUS[ifile][i], \
                    YYCLUS[ifile][i], ZZCLUS[ifile][i] ) )
            outf1.close()
#
    if( mode == 'verbose' ):
        print("The cluster configurations are written to *.xyz files") 
#
main()



