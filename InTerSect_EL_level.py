#

import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
from matplotlib.font_manager import FontProperties
import sys
from sympy import *
import sympy as sym
import os
from itertools import chain


# Intial candidates for fit, per FU: - thus, the E vs V input data has to be per FU
E0_init = -941.510817926696  # -1882.50963222/2.0 
V0_init = 63.54960592453 #125.8532/2.0 
B0_init = 76.3746233515232 #74.49 
B0_prime_init = 4.05340727164527 #4.15

#def BM(V, E0, V0, B0, B0_prime):
#        return  E0+ (2.293710449E+17)*(1E-21)*( (9.0/16.0)*(V0*B0) * (  (((V0/V)**(2.0/3.0)-1.0)**3.0)*B0_prime  + ((V0/V)**(2.0/3.0)-1)**2  *  (6.0-4.0*(V0/V)**(2.0/3.0))  ))

def BM(x, a, b, c, d):
         return  (2.293710449E+17)*(1E-21)* (a + b*x + c*x**2 + d*x**3 )

#def P(V, V0, B0, B0_prime):
#     f0=(3.0/2.0)*B0
#     f1=((V0/V)**(7.0/3.0))-((V0/V)**(5.0/3.0))
#     f2=((V0/V)**(2.0/3.0))-1
#     pressure= f0*f1*(1+(3.0/4.0)*(B0_prime-4)*f2)
#     return pressure

def P(x, b, c, d):
     return -b - 2*c*x - 3 *d*x**2

def H(x, a, b, c, d):
     return  a + b*x + c*x**2 + d*x**3

#filefolder_for_E_and_P_vs_V = '/home/david/Trabajo/structures/Calcite_I_and_II/PBE-D3__SHRINK_8_8__bipolar_18_18__TOLINTEG_8_18__XXLGRID_TOLDEE_8'

# Calcite I (Red triangles): 
V_not_p_f_unit_C_I, E_not_p_f_unit_C_I = np.loadtxt('./calcite_I.dat', skiprows = 1).T

# 14 (Empty grey triangles):
V_14_not_p_f_unit, E_14_not_p_f_unit = np.loadtxt('./calcite_II.dat', skiprows = 1).T

V_161_not_p_f_unit, E_161_not_p_f_unit = np.loadtxt('./calcite_I_son_SG_161.dat', skiprows = 1).T

# If the data is not per f unit, do this:
nFU_C_I = 2.0
nFU_C_II = 4.0
E_C_I = E_not_p_f_unit_C_I/nFU_C_I
V_C_I = V_not_p_f_unit_C_I/nFU_C_I

E_14 = E_14_not_p_f_unit/nFU_C_II
V_14 = V_14_not_p_f_unit/nFU_C_II

E_161 = E_161_not_p_f_unit/nFU_C_I
V_161 = V_161_not_p_f_unit/nFU_C_I

init_vals = [E0_init, V0_init, B0_init, B0_prime_init]

popt_C_I, pcov_C_I = curve_fit(BM, V_C_I, E_C_I, p0=init_vals)
popt_14, pcov_14 = curve_fit(BM, V_14, E_14, p0=init_vals)
popt_161, pcov_161 = curve_fit(BM, V_161, E_161, p0=init_vals)

# Linspace for plotting the fitting curves:
V_C_I_lin = np.linspace(V_C_I[0], V_C_I[-1], 10000)
V_14_lin = np.linspace(V_14[0], V_14[-1], 10000)
V_161_lin = np.linspace(V_161[0], V_161[-1], 10000)

# Plotting the fitting curves:
p2, = plt.plot(V_C_I_lin, BM(V_C_I_lin, *popt_C_I), color='black', label='Cubic fit Calcite I' )
p6, = plt.plot(V_14_lin, BM(V_14_lin, *popt_14), 'b', label='Cubic fit Calcite II')
p161, = plt.plot(V_161_lin, BM(V_161_lin, *popt_161), 'brown', label='Cubic fit 161')

# Plotting the scattered points: 
p1 = plt.scatter(V_C_I, E_C_I, color='red', marker="^", label='Calcite I', s=100)
p5 = plt.scatter(V_14, E_14, color='grey', marker="^", facecolors='none', label='Calcite II', s=100)

#p7 = plt.scatter(V_7, E_7, color='magenta', marker="^", facecolors='none', label='S.G. 7', s=100)

p8 = plt.scatter(V_161, E_161, color='green', marker="^", facecolors='none', label='S.G. 161', s=100)

fontP = FontProperties()
#fontP.set_size('small')
fontP.set_size('15')

plt.legend((p1, p2, p5, p6, p8, p161), ("Calcite I", "Cubic fit Calcite I", "Calcite II", 'Cubic fit Calcite II', "S.G. 161", "Cubic fit S.G. 161"), prop=fontP)

global V0, B0, B0_prime
E0   =     popt_C_I[0] 
V0   =     popt_C_I[1]
B0   =     popt_C_I[2]
B0_prime = popt_C_I[3]

pressures_per_F_unit_C_I = P(V_C_I, V0, B0, B0_prime)
print 'popt_C_I = ', popt_C_I
print 'popt_C_I[1:] = ', popt_C_I[1:]
output_array_2 = np.vstack((E_not_p_f_unit_C_I, V_not_p_f_unit_C_I, E_C_I, V_C_I, pressures_per_F_unit_C_I)).T
np.savetxt('Volumes_and_pressures_C_I.dat', output_array_2, header="Energy (a.u.) \t Volume (A^3) \t Energy / FU (a.u.) \t Volume / FU (A^3) \t Pressures (GPa)", fmt="%0.13f")

global V0_14, B0_14, B0_prime_14
E0_14   =     popt_14[0] 
V0_14   =     popt_14[1]
B0_14   =     popt_14[2]
B0_prime_14 = popt_14[3]


pressures_per_F_unit_14 = P(V_14, V0_14, B0_14, B0_prime_14)
output_array_2 = np.vstack((E_14_not_p_f_unit, V_14_not_p_f_unit, E_14, V_14, pressures_per_F_unit_14)).T
np.savetxt('Volumes_and_pressures_14.dat', output_array_2, header="Energy (a.u.) \t Volume (A^3) \t Energy / FU (a.u.) \t Volume / FU (A^3) \t Pressures (GPa)", fmt="%0.13f")

global V0_161, B0_161, B0_prime_161
E0_161   =     popt_161[0] 
V0_161   =     popt_161[1]
B0_161   =     popt_161[2]
B0_prime_161 = popt_161[3]


pressures_per_F_unit_161 = P(V_161, V0_161, B0_161, B0_prime_161)
output_array_3 = np.vstack((E_161_not_p_f_unit, V_161_not_p_f_unit, E_161, V_161, pressures_per_F_unit_161)).T
np.savetxt('Volumes_and_pressures_161.dat', output_array_3, header="Energy (a.u.) \t Volume (A^3) \t Energy / FU (a.u.) \t Volume / FU (A^3) \t Pressures (GPa)", fmt="%0.13f")

plt.xlabel('V / F.U. (Angstrom$^{3}$)', fontsize=20)
plt.ylabel(r'$E$ / F.U. (a.u.)', fontsize=20)
plt.suptitle("PBE-D3, pob-TZVP, SHRINK 8 8, Bipolar 18 18, TOLINTEG 8 18, XXLGRID, TOLDEE 8")
plt.title("(0.87 - 1.08)$V_{eq}$, 60 volumes", fontsize=10)
plt.ticklabel_format(useOffset=False)

plt.savefig('calcite_I_and_II_all_2_summary_better_plot.pdf', bbox_inches='tight')
plt.figure()


# Plotting the fitting curves:
p2, = plt.plot(V_C_I_lin, BM(V_C_I_lin, *popt_C_I), color='black', label='Cubic fit Calcite I' )

# Plotting the scattered points: 
p1 = plt.scatter(V_C_I, E_C_I, color='red', marker="^", label='Calcite I', s=100)
fontP = FontProperties()
#fontP.set_size('small')
fontP.set_size('15')

plt.legend((p1, p2), ("Calcite I", "Cubic fit Calcite I"), prop=fontP)

plt.xlabel('V / F.U. (Angstrom$^{3}$)', fontsize=20)
plt.ylabel(r'$E$ / F.U. (a.u.)', fontsize=20)
plt.suptitle("PBE-D3, pob-TZVP, SHRINK 8 8, Bipolar 18 18, TOLINTEG 8 18, XXLGRID, TOLDEE 8")
plt.title("(0.87 - 1.08)$V_{eq}$, 60 volumes", fontsize=10)
plt.ticklabel_format(useOffset=False)
plt.savefig('calcite_I_summary_better_plot.pdf', bbox_inches='tight')


# Plotting P vs V:
fig = plt.figure()

p2, = plt.plot(V_C_I_lin, P(V_C_I_lin, V0, B0, B0_prime), color='black', label='Cubic fit Calcite I' )
p6, = plt.plot(V_14_lin, P(V_14_lin, V0_14, B0_14, B0_prime_14), 'b', label='Cubic fit Calcite II')

p161, = plt.plot(V_161_lin, P(V_161_lin, V0_161, B0_161, B0_prime_161), 'brown', label='Cubic fit 161')


# Plotting the scattered points: 
p1 = plt.scatter(V_C_I, pressures_per_F_unit_C_I, color='red', marker="^", label='Calcite I', s=100)
p5 = plt.scatter(V_14, pressures_per_F_unit_14, color='grey', marker="^", facecolors='none', label='Calcite II', s=100)
#p7 = plt.scatter(V_7, E_7, color='magenta', marker="^", facecolors='none', label='S.G. 7', s=100)

p8 = plt.scatter(V_161, pressures_per_F_unit_161, color='green', marker="^", facecolors='none', label='S.G. 161', s=100)

fontP = FontProperties()
#fontP.set_size('small')
fontP.set_size('13')

plt.legend((p1, p2, p5, p6, p8, p161), ("Calcite I", "Cubic fit Calcite I", "Calcite II", 'Cubic fit Calcite II', "S.G. 161", "Cubic fit S.G. 161"), prop=fontP)


plt.xlabel('V / F.U. (Angstrom$^{3}$)', fontsize=20)
plt.ylabel(r'$P = -\frac{\partial E}{\partial V}$ (GPa)', fontsize=20)
plt.suptitle("PBE-D3, pob-TZVP, SHRINK 8 8, Bipolar 18 18, TOLINTEG 8 18, XXLGRID, TOLDEE 8")
plt.title("(0.87 - 1.08)$V_{eq}$, 60 volumes", fontsize=10)
plt.ticklabel_format(useOffset=False)

plt.savefig('calcite_I_and_II_all_2_summary_better_plot_P_vs_V.pdf', bbox_inches='tight')


#000000000000000000

H_C_I = E_C_I + pressures_per_F_unit_C_I * V_C_I * (2.293710449E+17)*(1E-21) 
H_14 = E_14 + pressures_per_F_unit_14 * V_14 * (2.293710449E+17)*(1E-21)
H_161 = E_161 + pressures_per_F_unit_161 * V_161 * (2.293710449E+17)*(1E-21)

output_array_3 = np.vstack((E_C_I, V_C_I, pressures_per_F_unit_C_I, H_C_I)).T
np.savetxt('E_V_P_H__C_I.dat', output_array_3, header="Energy / FU (a.u.) \t Volume / FU (A^3) \t Pressure / F.U. (GPa) \t Enthalpy (a.u.)", fmt="%0.13f") 

output_array_4 = np.vstack((E_14, V_14, pressures_per_F_unit_14, H_14)).T
np.savetxt('E_V_P_H__14.dat', output_array_4, header="Energy / FU (a.u.) \t Volume / FU (A^3) \t Pressure / F.U. (GPa) \t Enthalpy (a.u.)", fmt="%0.13f") 
#sys.exit()
output_array_5 = np.vstack((E_161, V_161, pressures_per_F_unit_161, H_161)).T
np.savetxt('E_V_P_H__161.dat', output_array_5, header="Energy / FU (a.u.) \t Volume / FU (A^3) \t Pressure / F.U. (GPa) \t Enthalpy (a.u.)", fmt="%0.13f") 

# Saving into variables:

P_lin_C_I = P(V_C_I_lin, V0, B0, B0_prime)
H_lin_C_I = BM(V_C_I_lin, *popt_C_I) +  P(V_C_I_lin, V0, B0, B0_prime) * V_C_I_lin * (2.293710449E+17)*(1E-21)

P_lin_14 = P(V_14_lin, V0_14, B0_14, B0_prime_14)
H_lin_14 = BM(V_14_lin, *popt_14) + P(V_14_lin, V0_14, B0_14, B0_prime_14) * V_14_lin * (2.293710449E+17)*(1E-21) 

P_lin_161 = P(V_161_lin, V0_161, B0_161, B0_prime_161)
H_lin_161 = BM(V_161_lin, *popt_161) + P(V_161_lin, V0_161, B0_161, B0_prime_161) * V_161_lin * (2.293710449E+17)*(1E-21) 


print ' P_lin_C_I = ', P_lin_C_I
print ' type(P_lin_C_I) = ', type(P_lin_C_I)
print ' H_lin_C_I = ', H_lin_C_I
print ' P_lin_14  = ', P_lin_14
print ' H_lin_14  = ', H_lin_14

output_array_1 = np.vstack((P_lin_C_I, H_lin_C_I)).T
np.savetxt('P_lin_C_I__H_lin_C_I.dat', output_array_1, header="P(GPa) \t   H per F unit (a.u)", fmt="%0.13f")

output_array_2 = np.vstack((P_lin_14, H_lin_14)).T
np.savetxt('P_lin_14__H_lin_14.dat', output_array_2, header="P(GPa) \t    H per F unit (a.u)", fmt="%0.13f")

output_array_3 = np.vstack((P_lin_161, H_lin_161)).T
np.savetxt('P_lin_14__H_lin_161.dat', output_array_3, header="P(GPa) \t    H per F unit (a.u)", fmt="%0.13f")


fig = plt.figure()

#  Obtaining a cubic expression for H(P):

# Reminder:
#H_C_I = E_C_I + pressures_per_F_unit_C_I * V_C_I * (2.293710449E+17)*(1E-21) 
#H_14 = E_14 + pressures_per_F_unit_14 * V_14 * (2.293710449E+17)*(1E-21)

init_vals = [E0_init, V0_init, B0_init, B0_prime_init]

popt_HofP_C_I, pcov_HofP_C_I = curve_fit(H, pressures_per_F_unit_C_I, H_C_I, p0=init_vals)
popt_HofP_14, pcov_HofP_14 = curve_fit(H, pressures_per_F_unit_14, H_14, p0=init_vals)
popt_HofP_161, pcov_HofP_161 = curve_fit(H, pressures_per_F_unit_161, H_161, p0=init_vals)


print " %%%%%%%%%%%%%%%%%%%%%% pressures_per_F_unit_161 = ", pressures_per_F_unit_161

pressures_per_F_unit_161_sorted = np.sort(pressures_per_F_unit_161)
print " %%%%%%%%%%%%%%%%%%%%%% np.sort(pressures_per_F_unit_161) ", pressures_per_F_unit_161_sorted
#print " %%%%%%%%%%%%%%%%%%%%%% np.sort(pressures_per_F_unit_161)[::-1] ", pressures_per_F_unit_161_sorted[::-1]
pressures_per_F_unit_161_lin = pressures_per_F_unit_161_sorted #[::-1]


# Linspace for plotting the fitting curves:
P_C_I_lin = np.linspace(pressures_per_F_unit_C_I[0], pressures_per_F_unit_C_I[-1], 10000)
P_14_lin = np.linspace(pressures_per_F_unit_14[0], pressures_per_F_unit_14[-1], 10000)
P_161_lin = np.linspace(pressures_per_F_unit_161_lin[0], pressures_per_F_unit_161_lin[-1], 10000)

# Plotting the fitting curves:
p2, = plt.plot(P_C_I_lin, H(P_C_I_lin, *popt_HofP_C_I), color='black', label='Cubic fit Calcite I' )
p6, = plt.plot(P_14_lin, H(P_14_lin, *popt_HofP_14), 'b', label='Cubic fit Calcite II')
p7, = plt.plot(P_161_lin, H(P_161_lin, *popt_HofP_161), 'brown', label='Cubic fit S.G. 161')

# Plotting the scattered points: 
p1 = plt.scatter(pressures_per_F_unit_C_I, H_C_I, color='red', marker="^", label='Calcite I', s=100)
p5 = plt.scatter(pressures_per_F_unit_14, H_14, color='grey', marker="^", facecolors='none', label='Calcite II', s=100)
p161 = plt.scatter(pressures_per_F_unit_161, H_161, color='green', marker="^", facecolors='none', label='Calcite II', s=100)

fontP = FontProperties()
#fontP.set_size('small')
fontP.set_size('13')

plt.legend((p1, p2, p5, p6, p161, p7), ("Calcite I", "Cubic fit Calcite I", "Calcite II", 'Cubic fit Calcite II', 'S.G. 161', 'Cubic fit S.G. 161'), prop=fontP)

global a0, a1, a2, a3
a0     =     popt_HofP_C_I[0] 
a1     =     popt_HofP_C_I[1]
a2     =     popt_HofP_C_I[2]
a3     =     popt_HofP_C_I[3]

global a0_s2, a1_s2, a2_s2, a3_s2
a0_s2        =     popt_HofP_14[0] 
a1_s2        =     popt_HofP_14[1]
a2_s2        =     popt_HofP_14[2]
a3_s2        =     popt_HofP_14[3]

print 'a0 = ', a0
print 'a1 = ', a1
print 'a2 = ', a2
print 'a3 = ', a3

print 'a0_s2 = ', a0_s2
print 'a1_s2 = ', a1_s2
print 'a2_s2 = ', a2_s2
print 'a3_s2 = ', a3_s2


print """ 
The equations are the following:
G_I (P) = a0 + a1*P + a2*P**2 + a3*P**3 
G_II (P) = a0_s2 + a1_s2*P + a2_s2*P**2 + a3_s2*P**3 
"""
print('G_I (P) = ({a0}) + ({a1})*P + ({a2})*P**2 + ({a3})*P**3 '.format(a0 = a0, a1 = a1, a2 = a2, a3 = a3, ))

print """
"""
print('G_II (P) = ({a0_s2}) + ({a1_s2})*P + ({a2_s2})*P**2 + ({a3_s2})*P**3 '.format(a0_s2 = a0_s2, a1_s2 = a1_s2, a2_s2 = a2_s2, a3_s2 = a3_s2))

print """
"""

print """
G_I (P) = G_II (P)
"""
# Set the boundaries for P here:
#P_C_I_lin
#P_14_lin

#z_fit = a0 + a1*P_C_I_lin + a2*P_C_I_lin**2 + a3*P_C_I_lin**3 		
#z_fit_2 = a0_s2 + a1_s2*P_14_lin + a2_s2*P_14_lin**2 + a3_s2*P_14_lin**3  

# Setting "P" to be symbolic:
P = sym.symbols('P') #, real=True)

def z_I(P):
        return   a0 + a1*P + a2*P**2 + a3*P**3 

def z_II(P):
        return   a0_s2 + a1_s2*P + a2_s2*P**2 + a3_s2*P**3 

#sol = sym.solve(z_I(P) - z_II(P) , P)
#print 'sol_ z_I(P) - z_II(P)  =', sol
#
#print 'sol[1] = ', sol[1]
#candidate = sol[2]
#P_real_intersection = re(candidate)
#print 'sol_ H_I(P) - H_II(P)[1]  =', P_real_intersection
#H_real_intersection = z_I(P_real_intersection)

# Crude intersection:
sol = sym.solve(z_I(P) - z_II(P) , P)
print 'sol_ H_I(P) - H_II(P)  =', sol

# Use of evalf to obtain better precision:
evalf_result = [x.evalf() for x in sol]
print '[x.evalf() for x in sol] = ', evalf_result

# Now, let's grab the real part of the evalf_result:
real_roots = []
for x in evalf_result:
  each_real_root = re(x)
  real_roots.append(each_real_root)

print 'real_roots = ', real_roots
for i in real_roots:
 print type(i)

# Transform each element of the list from <sympy.core.numbers.Float> to <float64>:
real_roots = [float(i) for i in real_roots]

for i in real_roots:
 print type(i)

# Transform each element of the list to a numpy array:
real_roots = np.array(real_roots)

for i in real_roots:
 print type(i)

print real_roots

# Let's grab the root located between 0.1GPa and 4GPa (true for CI-CII phase trans.)
real_roots_zero_to_four = real_roots[(real_roots >= 0.1) & (real_roots <= 4.0)]
print 'real_roots_zero_to_four = ', real_roots_zero_to_four

P_real_intersection = real_roots_zero_to_four
H_real_intersection = z_I(real_roots_zero_to_four)


plt.xlabel(r'$P$ (GPa)', fontsize=20)
plt.ylabel(r'$(H = E + PV)$ / F.U. (a.u.)', fontsize=15)
plt.suptitle("PBE-D3, pob-TZVP, SHRINK 8 8, Bipolar 18 18, TOLINTEG 8 18, XXLGRID, TOLDEE 8")
plt.title("(0.87 - 1.08)$V_{eq}$, 60 volumes.", fontsize=10)
plt.ticklabel_format(useOffset=False)
ax = fig.add_subplot(111)
ax.annotate('Analytic\nIntersection\nP= %g GPa\nG = %g a.u.' %(P_real_intersection, H_real_intersection), xy=(P_real_intersection, H_real_intersection), xytext=(P_real_intersection+2.5767, H_real_intersection-0.05), fontsize=15,
            arrowprops=dict(arrowstyle="->", connectionstyle="arc3", color='purple'),
            )
plt.savefig('calcite_I_and_II_all_2_summary_better_plot_delta_H_exact_expression_with_intersection.pdf', bbox_inches='tight')

plt.show()

