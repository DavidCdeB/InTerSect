# 
# This script extracts ET, TS EL E0 from many outputs and prints out:
# V vs F @ each Temperature.
# These files are needed for the last step of the QHA methodology,
# where we need tables of F vs V 
# in order to then compute P=dF/dV and then G(P,T)

import re
import os
import glob
from itertools import islice
import numpy as np
import sys

n_volume = []
path='./'
template = os.path.join(path, '*.out')

# Setting the number of formula units as a raw_input:
n_F_u = raw_input("""
Please type as an integer the number of formula units in the primitive cell. 
For example, Calcite I contains 2 formula units in the primitive (rombohedral) cell and 6 formula units in the crystallographic (hexagonal) cell. Thus, the number to be introduced is:   2 <and press ENTER>
""")

n_F_u = float(n_F_u)
n_F_u = int(float(n_F_u))

# Extracting each thermodynamic variable:
ET = []
TS = []
EL = []
E0 = []
VOLUME_EACH = []
T = []

for fname in glob.glob(template):
  print fname
  f = open(fname, 'r')
  real_part = False

  for line in f:

        if re.match(r"^ ET            :", line):
         start = line.find(':') + 8
         end = line.find(':') + 22
         result_ET = line[start:end]
         ET.append(result_ET)

        if re.match(r"^ TS            :", line):
         start = line.find(':') + 8
         end = line.find(':') + 22
         result_TS = line[start:end]
         TS.append(result_TS)

        if re.match(r"^ EL            :", line):
         start = line.find(':') + 4
         end = line.find(':') + 22
         result_EL = line[start:end]
         EL.append(result_EL)

        if re.match(r"^ E0            :", line):
         start = line.find(':') + 8
         end = line.find(':') + 22
         result_E0 = line[start:end]
         E0.append(result_E0)

        if re.match(r"^ AT \(T =", line):
         start = line.find('T =') + 4
         end = line.find('K')
         result_Temperatures = line[start:end]
         T.append(result_Temperatures)

        if 'LATTICE PARAMETERS  (ANGSTROMS AND DEGREES) - PRIMITIVE CELL' in line:  
                  print "line 1 = ", line
                  f.next()
                  each_volume_times_4 = []
                  each_volume_times_100 = []
                  
                  parameters = (''.join(islice(f, 1)))
                  columns = parameters.split()
                  each_volume = columns[6]
                  print 'each_volume = ', each_volume

                  VOLUME_EACH.append(each_volume)

# Transform each element of the list from <str> to <float64>:
VOLUME_EACH = [float(i) for i in VOLUME_EACH]
EL = [float(i) for i in EL]
E0 = [float(i) for i in E0]
ET = [float(i) for i in ET]
TS = [float(i) for i in TS]
T = [float(i) for i in T]

# Transform each element of the list to a numpy array:
VOLUME_EACH = np.array(VOLUME_EACH)
EL = np.array(EL)
E0 = np.array(E0)
ET = np.array(ET)
TS = np.array(TS)
T = np.array(T)

# Divide per F.U.:
VOLUME_EACH = VOLUME_EACH/n_F_u
EL = EL/n_F_u
E0 = E0/n_F_u
ET = ET/n_F_u
TS = TS/n_F_u

print 'VOLUME_EACH = ', VOLUME_EACH
print 'len(VOLUME_EACH) =  ' , len(VOLUME_EACH)

print 'EL = ', EL
print 'len(EL) =  ' , len(EL)

print 'E0 = ', E0
print 'len(E0) =  ' , len(E0)

print 'T = ', T
print 'len(T) =  ' , len(T)

print 'ET = ', ET
print 'len(ET) =  ' , len(ET)

print 'TS = ', TS
print 'len(TS) =  ' , len(TS)

#T_chosen = [10.00, 30.10, 70.30]

#ITEMS = []
#for t in T_chosen:
#  itemindex = np.argwhere(T==t)
## itemindex = np.where(T==t)
#  print itemindex
#  print type(itemindex)
#  itemindex = int(itemindex)
#  ITEMS.append(itemindex)
#print ITEMS

#for  i in ITEMS:
#    print ET[i]

#n_T = 4
#n_V = 2
#T = T[:n_T]

print 'ET.shape = ', ET.shape 

n_volume = len(VOLUME_EACH)
print 'n_volume = ', n_volume

n_T = len(T) / n_volume
print 'n_T = ', n_T

ET = np.reshape(ET, (n_volume, n_T))
TS = np.reshape(TS, (n_volume, n_T))
T  = np.reshape(T, (n_volume, n_T))

print 'ET reshaped = ', ET.shape 

print 'ET reduced = ', ET
print 'TS reduced = ', TS
print 'T reduced = ', T

rows = ET.shape[0]
cols = ET.shape[1]

print 'rows = ', rows
print 'cols = ', cols

rows = ET.shape[0]
cols = ET.shape[1]

F_all = []
for x, indx_EL, indx_E0 in zip(range(0, rows), range(len(EL)), range(len(E0))):
    print 'EL[indx_EL] =  ', EL[indx_EL]
    print 'E0[indx_E0] =  ', E0[indx_E0]
    print 'ET[x] =  ', ET[x]
    print 'TS[x] =  ', TS[x]
    aux = []    
    for y in range(0, cols):
        F = EL[indx_EL] + E0[indx_E0] + ET[x,y] - TS[x,y]
        print F
        aux.append(F)
    F_all.append(aux)

print ' F_all = ', F_all
print ' T = ', T
print ' VOLUME_EACH = ', VOLUME_EACH

for i in F_all:
 print F_all

# Transform to a np.array:
F_all = np.array(F_all)

print ' F_all = ', F_all
for i in F_all:
 print type(F_all)

print ' F_all[:,0] = ', F_all[:, 0]
print ' F_all[:,1] = ', F_all[:, 1]

cols_T = T.shape[1]
rows_T = T.shape[0]

print 'cols_T = ', cols_T
print 'rows = ', rows
print 'coles = ', cols

F_all_each_V_at_cte_T = []
for indx, t  in zip(range(0, cols), range(0, cols_T) ):
   aux_T = T[:,t] 
   aux_F = F_all[:,indx]
   print ' aux_F = ', aux_F
   print ' aux_T[0] = ', aux_T[0]

   output_array = np.vstack((VOLUME_EACH, aux_F)).T
   np.savetxt('F_vs_V_%0.2fK.dat'  %aux_T[0], output_array, header="Volume           F at %0.2fK" %aux_T[0], fmt="%0.13f")


