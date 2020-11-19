"""A Simple Plotting Script for Project 2 of CSE 701.
This plot shows a 2D surface and its divergence as a check for my AD."""
# This import registers the 3D projection, but is otherwise unused.
from mpl_toolkits.mplot3d import Axes3D

import matplotlib.pyplot as plt
import numpy as np

FIG = plt.figure()
# Will be for f(xi)
AX = FIG.add_subplot(121, projection='3d')
# Will be for f'(xi)
AXP = FIG.add_subplot(122, projection='3d')
# Load the data, skip the first/header row
DATA = np.loadtxt('file1.out', skiprows=1)
# Assign our data to a vector corresponding to the data
X = DATA.T[0]
Y = DATA.T[1]
F = DATA.T[2]
FP = DATA.T[3]

N = 30

# Setup linspace to span our independent variables
# NOTE: This is just going to produce a line that crosses the diagonal of an R^N space rather
# than array of whole space to compare with data; see todo
# TODO: Should *really* do a grid based on size of data and then compare that to the result

X_REAL = np.linspace(-2, 2, N+1)
Y_REAL = np.linspace(-2, 2, N+1)
XV, YV = np.meshgrid(X_REAL, Y_REAL)
FV = XV*np.sin(YV)
FPV = np.sin(YV)+XV*np.cos(YV)
# f_real = np.cos(x_real)*np.sin(y_real)
# fp_real = np.cos(x_real)*np.cos(y_real) - np.sin(x_real)*np.sin(y_real)
# F_REAL = X_REAL*np.sin(Y_REAL)
# FP_REAL = np.sin(Y_REAL)+X_REAL*np.cos(Y_REAL)

for m, zlow, zhigh in [('o', -50, -25)]:
    # Set f(xi) plot
    AX.scatter(X, Y, F, marker=m, label='AD F(x,y)')
    AX.scatter(XV, YV, FV, marker=m, label='numpy F(x,y)')
    # Set f'(xi) plot
    AXP.scatter(X, Y, FP, marker=m, label='AD F\'(x,y)')
    AXP.scatter(XV, YV, FPV, marker=m, label='numpy F\'(x,y)')

# print(F[::31]-FV)
# print(FP[::31]-FPV)

# Set labels for f' plot
AXP.legend()
AXP.set_title('F\' = sin(y) + xcos(y)')
AXP.set_xlabel('X')
AXP.set_ylabel('Y')
AXP.set_zlabel('F\'(X, Y)')

# Set labels for f plot
AX.legend()
AX.set_title('F = xsin(y)')
AX.set_xlabel('X')
AX.set_ylabel('Y')
AX.set_zlabel('F(X, Y)')

plt.show()
