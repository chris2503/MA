import numpy as np
import matplotlib.pyplot as plt
import scipy.constants as const
from scipy import optimize
import re
import os
import sys


m_det_or = np.array([34.21, 35.15, 34.63, 33.91, 35.50, 35.50, 35.50, 35.50, 35.50]) * 1e-3	#kg
m_detpaint_glyptal_or = np.array([0.06, 0.07, 0.06, 0.07]) * 1e-3							#kg
m_detpaint_epoxy_or= np.array([0.14, 0.13, 0.14, 0.12, 0.16]) * 1e-3								#kg

#m_det = np.mean(m_det_or)
#m_detpaint_glyptal = np.mean(m_detpaint_glyptal_or)
#m_detpaint_epoxy = np.mean(m_detpaint_epoxy_or)

rho_glyptal= 1.441*1e3		# kg/m^3
rho_epoxy = 1.135*1e3

x_D = 10.2*1e-3				# m													
z_D = 16.0*1e-3				# m	


O = x_D**2 + z_D*x_D*4
V_glyptal = []
d_glyptal = []
V_epoxy = []
d_epoxy = []
for i in range(len(m_detpaint_glyptal_or)):
	V_glyptal.append(m_detpaint_glyptal_or[i]/rho_glyptal)
	d_glyptal.append(V_glyptal[i]/O)

for i in range(len(m_detpaint_epoxy_or)):
	V_epoxy.append(m_detpaint_epoxy_or[i]/rho_epoxy)
	d_epoxy.append(V_epoxy[i]/O)
print("\nGlyptal-Lackdicke: ", d_glyptal)
print("Epoxy-Lackdicke: ", d_epoxy)

print()
print("Means:")
print("Glyptal-Lackdicke: %f" %(np.mean(d_glyptal)))
print("Epoxy-Lackdicke: %f" %(np.mean(d_epoxy)))