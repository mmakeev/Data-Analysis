#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 15:05:14 2019

@author: maximmakeev
"""

import sys
import os
import platform
#import math
#
#
#
print ( sys.version )
#
#
ppath = sys.argv[0] 
print ( "Number of arguments\n", ppath )
#-----------------------------------------------------
print( os.name )
print( platform.system() )
print( platform.release() )
#
#-----------------------------------------------------
currwd = os.getcwd( )
print("Current directory is:\n", currwd)
fileinput = open( "input.file","r" )
#
AI = fileinput.readlines()
#
fileinput.close()
#----Parameters-loaded from the input file------------
#
tmp01 = AI[1].split()
cf_name = str( tmp01[0] )
tmp02 = AI[3].split()
mypath_set = str( tmp02[0] )
mypath = str( mypath_set )
os.chdir( mypath )
workdir = os.getcwd() 
print( workdir )
#-----------------------------------------------------
filelist = []
#-----------------------------------------------------
##A = 16350
NA = 130800
IDD = [ 0.0 for ii in range( 0,NA ) ]
QAT = [ 0.0 for ii in range( 0,NA ) ]
RXX = [ 0.0 for ii in range( 0,NA ) ]
RYY = [ 0.0 for ii in range( 0,NA ) ]
RZZ = [ 0.0 for ii in range( 0,NA)  ]
CLIST = []
#-----------------------------------------------------
#--s^{2} C^{2}/(m^{3} kg) [J/K] [m^3]
#-- ( 10^{-20}*10^{-38} ) / (10^{-30}*10^{-12}*10^{-23}*10^{5} )
con = pow(1.60217662,2)*100.0/( 8.85*1.3806485279*9.0 )
#-----------------------------------------------------
mdirs = os.scandir( workdir )
for entry in mdirs:
    if entry.is_file():
        cfile = entry.name
        if cfile.startswith( cf_name ):
            filelist.append( cfile )
print( len(filelist) )
nfile = len( filelist )
filelist.sort(key = len)
for entry in range(nfile):
    print( filelist[entry] )
PAVX = [0.0 for ii in range( 0,nfile) ]
PAVY = [0.0 for ii in range( 0,nfile) ]
PAVZ = [0.0 for ii in range( 0,nfile) ]
P2AV = [0.0 for ii in range( 0,nfile) ]
#-----------------------------------------------------
for nfi in range( nfile ):
#    print("file", nfi )
    filename = filelist[ nfi ]
#    print( filename )
    fileseq = open( filename, "r" )
    A = fileseq.readlines()
    fileseq.close()
    totnum = len( A )
    if nfi == 0:
        l1 = A[5].split() 
        l2 = A[6].split() 
        l3 = A[7].split()
        Lx = float( l1[1] ) - float( l1[0] )/10.0
        Ly = float( l2[1] ) - float( l2[0] )/10.0
        Lz = float( l3[1] ) - float( l3[0] )/10.0
        VV = Lx*Ly*Lz
        natom = int( A[3].split()[0] )
        print( natom )
        numa = 0
    else:
        numa = 0
        for nn in range(9, totnum):
            tmp = A[nn].split() 
            IDD[numa] =   int( tmp[1] )
            QAT[numa] = float( tmp[2] )
            RXX[numa] = float( tmp[3] ) - Lx*0.25
            RYY[numa] = float( tmp[4] ) - Lx*0.25
            RZZ[numa] = float( tmp[5] ) - Lx*0.25
            numa += 1
    summ = float( 0.0 )
    sumX = float( 0.0 )
    sumY = float( 0.0 )
    sumZ = float( 0.0 )
    QTOT = float( 0.0 )
    for nat in range( 0, NA ):
            QTOT = QTOT + QAT[nat]
            qcoo1 = QAT[nat]*RXX[nat] 
            qcoo2 = QAT[nat]*RYY[nat] 
            qcoo3 = QAT[nat]*RZZ[nat] 
            sumX = sumX + qcoo1
            sumY = sumY + qcoo2
            sumZ = sumZ + qcoo3
    summ = sumX*sumX + sumY*sumY + sumZ*sumZ
    PAVX[nfi] = sumX*(con/VV)
    PAVY[nfi] = sumY*(con/VV)
    PAVZ[nfi] = sumZ*(con/VV)
    P2AV[nfi] = summ*(con/VV)
    print( nfi, P2AV[nfi], QTOT, "\n"  )
fsumX = float( 0.0 )
fsumY = float( 0.0 )
fsumZ = float( 0.0 )
fsum2 = float( 0.0 )
for nf in range( 1, nfile ):
    fsumX = fsumX + PAVX[nf]
    fsumY = fsumY + PAVY[nf]
    fsumZ = fsumZ + PAVZ[nf]
    fsum2 = fsum2 + P2AV[nf]
fsumX = fsumX/float( nfile-1 )
fsumY = fsumY/float( nfile-1 )
fsumZ = fsumZ/float( nfile-1 )
fsum1 = fsumX*fsumX + fsumY*fsumY + fsumZ*fsumZ
diel  = fsum2/float( nfile-1 ) - fsum1 + 1.0
print( "dielectric constant =", diel, QTOT )


        
        
    

