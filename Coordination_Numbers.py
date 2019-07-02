#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Jun  2 20:12:03 2019

@author: maximmakeev
"""

import os
import sys
import math
#
#
#-------------------------------------------------------------------
def main():
#---Partial RDFs are computed between types: nrdf[n] and nrdf[m]----
# 
    nrdf0 = []
    nrdf1 = []
    mass  = []
#-------------------------------------------------------------------
#Containers for atomic ID, type and Cartesian coordinates
#-------------------------------------------------------------------
    currwd = os.getcwd( )
    print("Current directory is:\n", currwd)
    rdfinput = open( "i_coord.dat","r" )
#
    AI = rdfinput.readlines()
#
    rdfinput.close()
#----Parameters-loaded from the input file-------------------------
#
    tmp01 = AI[1].split()
    rcut = float( tmp01[0] )
    tmp02 = AI[3].split()
    ddr = float( tmp02[0] )
    tmp03 = AI[5].split()
    NA = int ( tmp03[0] )
    tmp04 = AI[7].split()
    ntypes = int ( tmp04[0] )
    tmp05 = AI[10].split()
    for kk in range( 0,ntypes ):
        mass.append( tmp05[kk] )
    tmp06 = AI[12].split()
    nrdfs = int( tmp06[0] )
    tmp07 = AI[15].split()
    tmp08 = AI[16].split()
    for ll in range( 0,nrdfs ):
        nrdf0.append( tmp07[ll] )
        nrdf1.append( tmp08[ll] )
    tmp09 = AI[19].split()
    input_mode = tmp09[0]
    filename   = tmp09[1]
    tmp10 = AI[21].split()
    workdir0 = tmp10[0]
#
    os.chdir( workdir0 )
    pathdir = os.getcwd()
    print("Working directory is:")
    print( pathdir )
#
#----Find the MD trajectory files to be used for calculation-------
#
#---For configurations in multiple files, count the number of files
    if (input_mode == 'multi' ):
        listFiles = [ ]
        listoffiles = os.listdir()
        for entry in listoffiles:
            if entry.startswith( 'Mg_2TFSI_G1.lammpstrj.' ):
                listFiles.append(entry)
                num_files = len( listFiles )
        print("# of files/frames (num_files) read:", num_files )
#
#---For configurations in a single file, count the number of frames
#
    num_frame_tot = int( 0 )
    if ( input_mode == 'single' ):
        filenameONE = filename
        rdffile = open( filenameONE,"r" )
        B = rdffile.readlines()
        rdffile.close()
        file_size = len( B )
        print( file_size )
        for ii in range( 0,file_size ):
            tmp = B[ii].split()
            if ( len( tmp ) > 1 and tmp[1] == "TIMESTEP" ): num_frame_tot += 1
            n_per_frame = int( file_size/num_frame_tot )
            num_files = int( num_frame_tot )
        print( "num_files/frames (num_files) read:=:", num_files )
#-----------------------------------------------------------------
#
#--Define containers to be used to sum configurations------------- 
    CON_P_SUM = [0.0]*len(nrdf0)
#
#----Loop over the coordinate files-------------------------------
#----Loop over trajectory files-----------------------------------
    for i_tr in range( 0,int( 1 ) ):
        if (input_mode == 'multi' ):
            filenameT = listFiles[ i_tr ]
            print( "Processing file:", filenameT )
            rdffile = open( filenameT,"r" )
            A = rdffile.readlines()   
            rdffile.close()
        if (input_mode == 'single' ):
            filenameT = filenameONE
            A = []
            print( "Processing frame:", i_tr )
            start_line = int(n_per_frame)*(i_tr)
            end_line   = int(n_per_frame)*int(i_tr + int(1) ) 
            for jjj in range( start_line, end_line ):
                A.append( B[ jjj ] ) 
        nlines = len( A )
        print( nlines)
        lco = [0]*5
        natom = 0
        key1 = 'ATOMS'
        ATID = []
        ATTY = []
        ATXX = []
        ATYY = []
        ATZZ = []
        for i in range( 0,nlines ):
            tmp = A[i].split()
            if len(tmp) > 1 and tmp[1] == key1:
                m = len( tmp )
                for n in range( 0,m ):
                    if tmp[n] == 'id':
                        lco[0] = (n-2)
                    if tmp[n] == 'type':
                        lco[1] = (n-2)
                    if tmp[n] == 'x':
                        lco[2] = (n-2)
                    if tmp[n] == 'y':
                        lco[3] = (n-2)
                    if tmp[n] == 'z':
                        lco[4] = (n-2)
            if( i == 3): NA = int( tmp[0] )
            if( i == 3): print( 'Number of atoms=: {0:d}'.format(NA) )
            if( i == 5): Lx = float( tmp[1] ) - float( tmp[0] )
            if( i == 6): Ly = float( tmp[1] ) - float( tmp[0] )
            if( i == 7): Lz = float( tmp[1] ) - float( tmp[0] )
            if( i >= 9 ):
                natom += 1
                ATID.append( int(  tmp[lco[0]] ) )
                ATTY.append( int(  tmp[lco[1]] ) )
                ATXX.append( float(tmp[lco[2]] ) )
                ATYY.append( float(tmp[lco[3]] ) )
                ATZZ.append( float(tmp[lco[4]] ) )
#
        if( NA != natom):
            print( "Major consistency check failed:" )
            print( "Configuration was not read correctly." )
            print( "NA=:", NA, "natoms=:", natom)
            sys.exit()
        n_a_pairs = len( nrdf0 )
        setID = {*()}
        atomtypes = [0]*( ntypes+1 )
#
        for jj in range( 0,natom ):
            npp = ATTY[jj]
            setID.add( npp )
            atomtypes[npp] += 1
            nset = len( setID )
#------------------------------------------------------------------------
        if( ntypes != nset):
            print( "Consistency check failed:" )
            print( "Number of atomic types in the config file is \
                  different from the corresponding value in input file" )
            print( "ntypes=:", ntypes, "nset=:", nset)
#------------------------------------------------------------------------
#                    
        CON_P = [0.0]*len( nrdf0 )
# 
#------Coordination Numbers Calculations---------------------------------
        for i in range( 0,natom ):
            if( i%1000 == 0): 
                print("Processing atom #:", i )
            for j in range(i+1, natom):
                dx = (ATXX[i] - ATXX[j])
                dy = (ATYY[i] - ATYY[j])
                dz = (ATZZ[i] - ATZZ[j])
                if( dx >  float(0.5)*Lx ): dx = dx - Lx
                if( dy >  float(0.5)*Ly ): dy = dy - Ly
                if( dz >  float(0.5)*Lz ): dz = dz - Lz
                if( dx < -float(0.5)*Lx ): dx = dx + Lx
                if( dy < -float(0.5)*Ly ): dy = dy + Ly
                if( dz < -float(0.5)*Lz ): dz = dz + Lz
                rsq = math.sqrt( dx*dx + dy*dy + dz*dz )
                if( rsq < rcut ):
                    for kl in range( 0,n_a_pairs ):
                        nta1 = int( nrdf0[kl] ) 
                        nta2 = int( nrdf1[kl] )
                        if( int(ATTY[i]) == nta1 and int(ATTY[j]) == nta2 ):
                            CON_P[ kl ] = CON_P[ kl ] + 1.0
                        if( int(ATTY[j]) == nta1 and int(ATTY[i]) == nta2 ):
                            CON_P[ kl ] = CON_P[ kl ] + 1.0
#          
        print( "Done computing COOR_NUM:for Mg_2TFSI_G1.lammpstrj", i_tr )
#---Normalization Procedure for the full RDF and partical RDFs
#
        for kl in range( 0, n_a_pairs ):
            npp = int( nrdf0[ kl ] )
            print (CON_P[kl])
            CON_P[ kl ] = CON_P[ kl ]/atomtypes[ npp ]                     
#    
#----------------------------------------------------------------------      
        for kl in range( 0,n_a_pairs ):
           CON_P_SUM[ kl ] =  CON_P_SUM[ kl ] + CON_P[kl]  
#
    for kl in range( 0,n_a_pairs ):
        CON_P_SUM[ kl ] =  CON_P_SUM[ kl ]/float( 1.0 ) 
#
#---------------------------------------------------------------------
#
    outf = open("coord_number.dat","w+")
    for kk in range( 0, n_a_pairs ):
        p1 = nrdf0[ kk ]
        p2 = nrdf1[ kk ]
        outf.write("%d %d %25.20f\n" % ( int(p1), int(p2),\
                                        CON_P_SUM[ kk ] ) )
    outf.close()
#---------------------------------------------------------------------
#
    print("Coordination nambers are written in coord_number.dat file")     
  
main()   
    
    
    
    
    
    
    
    
    
    
    
    

