

Get Outlook for iOS
________________________________________
From: Cameron Smith
Sent: Monday, September 1, 2025 1:35:28 PM
To: Oliver Chu <och42@uclive.ac.nz>
Subject: 302 assessment 
 
Hey mate,
 
Here is all the code that has been used. The report is also attached. I think I have all the data needed and diagrams. Still need to add in the overall drawings and assembly matrices. The writing in there is purely ChatGPT that we can use for inspiration as to what needs to be written. 
 
 302 - Ollie and Cam.docx
 
Frame code
"""
Created on Thu Jul 31 14:16:22 2025
 
@author: camsmith
"""
import numpy as np 
import matplotlib.pyplot as plt
 
def local_frame(E, A, L, I): #This matrix stays the same always
    stiff = (E * I) / L**3 #Sets the stiffness for the element
    beta = A * L**2 / I
    ke = stiff * np.array([
    [ beta,     0,       0,     -beta,     0,       0     ],
    [ 0,       12,     6*L,       0,     -12,     6*L    ],
    [ 0,      6*L,   4*L**2,      0,    -6*L,   2*L**2   ],
    [-beta,     0,       0,      beta,     0,       0     ],
    [ 0,     -12,    -6*L,        0,      12,    -6*L    ],
    [ 0,      6*L,   2*L**2,      0,    -6*L,   4*L**2   ]
])
    return ke
    
def upplambda(theta): #Uppercase lambda converts from local to global coord system
    c = np.cos(theta)
    s = np.sin(theta)
    lam = np.array([[c, s, 0], [-s, c, 0], [0, 0, 1]])
    upplambda = np.block([
        [lam,               np.zeros((3,3))],
        [np.zeros((3,3)),   lam]
    ])
    
    return upplambda
 
def global_frame(ke_local, upplambda):
    return upplambda.T @ ke_local @ upplambda #Returning the stiffness matrix in global coord
 
I = 13628778.666666666e-12  # = 1.3628778667e-5 m^4
A = 4.496e-3
 
"""Element 1"""
E1, A1, L1, I1 = 200e9, A, 2.5, I #Material Props
theta_deg1 = 90 #Angle roatated
theta_rad1 = np.radians(theta_deg1)
ke1 = local_frame(E1, A1, L1, I1)
khat1 = global_frame(ke1, upplambda(theta_rad1))
             #D#1  2  3  4  5  6
A1 = np.array([[0, 0, 0, 1, 0, 0], #q1
               [0, 0, 0, 0, 1, 0], #q2
               [0, 0, 0, 0, 0, 1], #q3
               [0, 0, 0, 0, 0, 0], #q4
               [0, 0, 0, 0, 0, 0], #q5
               [0, 0, 0, 0, 0, 0], #q6
               [0, 0, 0, 0, 0, 0], #q7
               [0, 0, 0, 0, 0, 0], #q8
               [0, 0, 0, 0, 0, 0], #q9
               [0, 0, 0, 0, 0, 0], #q10
               [0, 0, 0, 0, 0, 0], #q11
               [0, 0, 0, 0, 0, 0], #q12
               [0, 0, 0, 0, 0, 0], #q13
               [0, 0, 0, 0, 0, 0], #q14
               [0, 0, 0, 0, 0, 0]]) #q15
         
"""Element 2"""
E2, A2, L2, I2 = 200e9, A, 2.5, I #Material Props
theta_deg2 = 90 #Angle roatated
theta_rad2 = np.radians(theta_deg2)
ke2 = local_frame(E2, A2, L2, I2)
khat2 = global_frame(ke2, upplambda(theta_rad2))
             #D#1  2  3  4  5  6
A2 = np.array([[1, 0, 0, 0, 0, 0], #q1
               [0, 1, 0, 0, 0, 0], #q2
               [0, 0, 1, 0, 0, 0], #q3
               [0, 0, 0, 1, 0, 0], #q4
               [0, 0, 0, 0, 1, 0], #q5
               [0, 0, 0, 0, 0, 1], #q6
               [0, 0, 0, 0, 0, 0], #q7
               [0, 0, 0, 0, 0, 0], #q8
               [0, 0, 0, 0, 0, 0], #q9
               [0, 0, 0, 0, 0, 0], #q10
               [0, 0, 0, 0, 0, 0], #q11
               [0, 0, 0, 0, 0, 0], #q12
               [0, 0, 0, 0, 0, 0], #q13
               [0, 0, 0, 0, 0, 0], #q14
               [0, 0, 0, 0, 0, 0]]) #q15
         
               
"""Element 3"""
E3, A3, L3, I3 = 200e9, A, 2.5, I #Material Props
theta_deg3 = 90 #Angle roatated
theta_rad3 = np.radians(theta_deg3)
ke3 = local_frame(E3, A3, L3, I3)
khat3 = global_frame(ke3, upplambda(theta_rad3))
             #D#1  2  3  4  5  6
A3 = np.array([[0, 0, 0, 0, 0, 0], #q1
               [0, 0, 0, 0, 0, 0], #q2
               [0, 0, 0, 0, 0, 0], #q3
               [1, 0, 0, 0, 0, 0], #q4
               [0, 1, 0, 0, 0, 0], #q5
               [0, 0, 1, 0, 0, 0], #q6
               [0, 0, 0, 1, 0, 0], #q7
               [0, 0, 0, 0, 1, 0], #q8
               [0, 0, 0, 0, 0, 1], #q9
               [0, 0, 0, 0, 0, 0], #q10
               [0, 0, 0, 0, 0, 0], #q11
               [0, 0, 0, 0, 0, 0], #q12
               [0, 0, 0, 0, 0, 0], #q13
               [0, 0, 0, 0, 0, 0], #q14
               [0, 0, 0, 0, 0, 0]]) #q15
         
 
"""Element 4"""
E4, A4, L4, I4 = 200e9, A, np.sqrt(2.5**2 + 1**2), I #Material Props
theta_deg4 = -68.2 #Angle roatated
theta_rad4 = np.radians(theta_deg4)
ke4 = local_frame(E4, A4, L4, I4)
khat4 = global_frame(ke4, upplambda(theta_rad4))
             #D#1  2  3  4  5  6
A4 = np.array([[1, 0, 0, 0, 0, 0], #q1
               [0, 1, 0, 0, 0, 0], #q2
               [0, 0, 1, 0, 0, 0], #q3
               [0, 0, 0, 0, 0, 0], #q4
               [0, 0, 0, 0, 0, 0], #q5
               [0, 0, 0, 0, 0, 0], #q6
               [0, 0, 0, 0, 0, 0], #q7
               [0, 0, 0, 0, 0, 0], #q8
               [0, 0, 0, 0, 0, 0], #q9
               [0, 0, 0, 0, 0, 0], #q10
               [0, 0, 0, 0, 0, 0], #q11
               [0, 0, 0, 0, 0, 0], #q12
               [0, 0, 0, 0, 0, 0], #q13
               [0, 0, 0, 0, 0, 0], #q14
               [0, 0, 0, 0, 0, 0]]) #q15
         
 
"""Element 5"""
E5, A5, L5, I5 = 200e9, A, 1, I #Material Props
theta_deg5 = 0 #Angle roatated
theta_rad5 = np.radians(theta_deg5)
ke5 = local_frame(E5, A5, L5, I5)
khat5 = global_frame(ke5, upplambda(theta_rad5))
             #D#1  2  3  4  5  6
A5 = np.array([[1, 0, 0, 0, 0, 0], #q1
               [0, 1, 0, 0, 0, 0], #q2
               [0, 0, 1, 0, 0, 0], #q3
               [0, 0, 0, 0, 0, 0], #q4
               [0, 0, 0, 0, 0, 0], #q5
               [0, 0, 0, 0, 0, 0], #q6
               [0, 0, 0, 0, 0, 0], #q7
               [0, 0, 0, 0, 0, 0], #q8
               [0, 0, 0, 0, 0, 0], #q9
               [0, 0, 0, 0, 0, 0], #q10
               [0, 0, 0, 0, 0, 0], #q11
               [0, 0, 0, 0, 0, 0], #q12
               [0, 0, 0, 1, 0, 0], #q13
               [0, 0, 0, 0, 1, 0], #q14
               [0, 0, 0, 0, 0, 1]]) #q15
         
 
"""Element 6"""
E6, A6, L6, I6 = 200e9, A, np.sqrt(2.5**2 + 1**2), I #Material Props
theta_deg6 = 68.2 #Angle roatated
theta_rad6 = np.radians(theta_deg6)
ke6 = local_frame(E6, A6, L6, I6)
khat6 = global_frame(ke6, upplambda(theta_rad6))
             #D#1  2  3  4  5  6
A6 = np.array([[1, 0, 0, 0, 0, 0], #q1
               [0, 1, 0, 0, 0, 0], #q2
               [0, 0, 1, 0, 0, 0], #q3
               [0, 0, 0, 0, 0, 0], #q4
               [0, 0, 0, 0, 0, 0], #q5
               [0, 0, 0, 0, 0, 0], #q6
               [0, 0, 0, 0, 0, 0], #q7
               [0, 0, 0, 0, 0, 0], #q8
               [0, 0, 0, 0, 0, 0], #q9
               [0, 0, 0, 1, 0, 0], #q10
               [0, 0, 0, 0, 1, 0], #q11
               [0, 0, 0, 0, 0, 1], #q12
               [0, 0, 0, 0, 0, 0], #q13
               [0, 0, 0, 0, 0, 0], #q14
               [0, 0, 0, 0, 0, 0]]) #q15
         
 
"""Element 7"""
E7, A7, L7, I7 = 200e9, A, 1, I #Material Props
theta_deg7 = 0 #Angle roatated
theta_rad7 = np.radians(theta_deg7)
ke7 = local_frame(E7, A7, L7, I7)
khat7 = global_frame(ke7, upplambda(theta_rad7))
             #D#1  2  3  4  5  6
A7 = np.array([[0, 0, 0, 0, 0, 0], #q1
               [0, 0, 0, 0, 0, 0], #q2
               [0, 0, 0, 0, 0, 0], #q3
               [1, 0, 0, 0, 0, 0], #q4
               [0, 1, 0, 0, 0, 0], #q5
               [0, 0, 1, 0, 0, 0], #q6
               [0, 0, 0, 0, 0, 0], #q7
               [0, 0, 0, 0, 0, 0], #q8
               [0, 0, 0, 0, 0, 0], #q9
               [0, 0, 0, 1, 0, 0], #q10
               [0, 0, 0, 0, 1, 0], #q11
               [0, 0, 0, 0, 0, 1], #q12
               [0, 0, 0, 0, 0, 0], #q13
               [0, 0, 0, 0, 0, 0], #q14
               [0, 0, 0, 0, 0, 0]]) #q15
         
 
"""Element 8"""
E8, A8, L8, I8 = 200e9, A, 3, I #Material Props
theta_deg8 = -68.2 #Angle roatated
theta_rad8 = np.radians(theta_deg8)
ke8 = local_frame(E8, A8, L8, I8)
khat8 = global_frame(ke8, upplambda(theta_rad8))
             #D#1  2  3  4  5  6
A8 = np.array([[0, 0, 0, 0, 0, 0], #q1
               [0, 0, 0, 0, 0, 0], #q2
               [0, 0, 0, 0, 0, 0], #q3
               [0, 0, 0, 0, 0, 0], #q4
               [0, 0, 0, 0, 0, 0], #q5
               [0, 0, 0, 0, 0, 0], #q6
               [1, 0, 0, 0, 0, 0], #q7
               [0, 1, 0, 0, 0, 0], #q8
               [0, 0, 1, 0, 0, 0], #q9
               [0, 0, 0, 1, 0, 0], #q10
               [0, 0, 0, 0, 1, 0], #q11
               [0, 0, 0, 0, 0, 1], #q12
               [0, 0, 0, 0, 0, 0], #q13
               [0, 0, 0, 0, 0, 0], #q14
               [0, 0, 0, 0, 0, 0]]) #q15
         
 
"""Element 9"""
E9, A9, L9, I9 = 200e9, A, 2.5, I #Material Props
theta_deg9 = 90 #Angle roatated
theta_rad9 = np.radians(theta_deg9)
ke9 = local_frame(E9, A9, L9, I9)
khat9 = global_frame(ke9, upplambda(theta_rad9))
             #D#1  2  3  4  5  6
A9 = np.array([[0, 0, 0, 0, 0, 0], #q1
               [0, 0, 0, 0, 0, 0], #q2
               [0, 0, 0, 0, 0, 0], #q3
               [0, 0, 0, 0, 0, 0], #q4
               [0, 0, 0, 0, 0, 0], #q5
               [0, 0, 0, 0, 0, 0], #q6
               [0, 0, 0, 0, 0, 0], #q7
               [0, 0, 0, 0, 0, 0], #q8
               [0, 0, 0, 0, 0, 0], #q9
               [0, 0, 0, 0, 0, 0], #q10
               [0, 0, 0, 0, 0, 0], #q11
               [0, 0, 0, 0, 0, 0], #q12
               [0, 0, 0, 1, 0, 0], #q13
               [0, 0, 0, 0, 1, 0], #q14
               [0, 0, 0, 0, 0, 1]]) #q15
         
 
"""Element 10"""
E10, A10, L10, I10 = 200e9, A, 2.5, I #Material Props
theta_deg10 = 90 #Angle roatated
theta_rad10 = np.radians(theta_deg10)
ke10 = local_frame(E10, A10, L10, I10)
khat10 = global_frame(ke10, upplambda(theta_rad10))
             #D#1  2  3  4  5  6
A10 = np.array([[0, 0, 0, 0, 0, 0], #q1
                [0, 0, 0, 0, 0, 0], #q2
                [0, 0, 0, 0, 0, 0], #q3
                [0, 0, 0, 0, 0, 0], #q4
                [0, 0, 0, 0, 0, 0], #q5
                [0, 0, 0, 0, 0, 0], #q6
                [0, 0, 0, 0, 0, 0], #q7
                [0, 0, 0, 0, 0, 0], #q8
                [0, 0, 0, 0, 0, 0], #q9
                [0, 0, 0, 1, 0, 0], #q10
                [0, 0, 0, 0, 1, 0], #q11
                [0, 0, 0, 0, 0, 1], #q12
                [1, 0, 0, 0, 0, 0], #q13
                [0, 1, 0, 0, 0, 0], #q14
                [0, 0, 1, 0, 0, 0]]) #q15
 
"""Indivdual KG Matrixs"""
KG1 = A1 @ khat1 @ A1.T
KG2 = A2 @ khat2 @ A2.T
KG3 = A3 @ khat3 @ A3.T
KG4 = A4 @ khat4 @ A4.T
KG5 = A5 @ khat5 @ A5.T
KG6 = A6 @ khat6 @ A6.T
KG7 = A7 @ khat7 @ A7.T
KG8 = A8 @ khat8 @ A8.T
KG9 = A9 @ khat9 @ A9.T
KG10 = A10 @ khat10 @ A10.T
 
"""Total"""
KG = KG1 + KG2 + KG3 + KG4 + KG5 + KG6 + KG7 + KG8 + KG9 + KG10
 
 
"""Q Loads"""
# Need to change the relative loads, length, assembly matrix and angle
# All values that need to be changed are in th top 4 lines of each load
 
 
"""Nodal Q Vector""" #Enter in Newtons
Q = np.array([[15000],
              [0], 
              [0], 
              [30000], 
              [0],
              [0],
              [15000],
              [0], 
              [0], 
              [0], 
              [0],
              [0],
              [0],
              [0],
              [0]])
 
 
"""Q loads from structure"""
Qtot = Q #Q_udl #+ Q_pl + Q_lv + Q_lvl + Q_pla + Q_dax + Q_ax #add in required loads
q = np.linalg.solve(KG,  Qtot)
 
"""Forces"""
F1 = khat1 @ A1.T @ q
F2 = khat2 @ A2.T @ q
F3 = khat3 @ A3.T @ q
F4 = khat4 @ A4.T @ q
F5 = khat5 @ A5.T @ q
F6 = khat6 @ A6.T @ q
F7 = khat7 @ A7.T @ q
F8 = khat8 @ A8.T @ q
F9 = khat9 @ A9.T @ q
F10 = khat10 @ A10.T @ q
 
"""Deflections"""
d1 = upplambda(theta_rad1) @ A1.T @ q
d2 = upplambda(theta_rad2) @ A2.T @ q
d3 = upplambda(theta_rad3) @ A3.T @ q
d4 = upplambda(theta_rad4) @ A4.T @ q
d5 = upplambda(theta_rad5) @ A5.T @ q
d6 = upplambda(theta_rad6) @ A6.T @ q
d7 = upplambda(theta_rad7) @ A7.T @ q
d8 = upplambda(theta_rad8) @ A8.T @ q
d9 = upplambda(theta_rad9) @ A9.T @ q
d10 = upplambda(theta_rad10) @ A10.T @ q
 
"""Internal forces"""
f1 = ke1 @ d1
f2 = ke2 @ d2
f3 = ke3 @ d3
f4 = ke4 @ d4
f5 = ke5 @ d5
f6 = ke6 @ d6
f7 = ke7 @ d7
f8 = ke8 @ d8
f9 = ke9 @ d9
f10 = ke10 @ d10
 
"""Tricky plot"""
def plot_deflected_shape(node1XG, node1YG, node2XG, node2YG, d_e, Disp_mag, N_points, L, alpha):
    """plot full deflected shape"""
    x_e = np.linspace(0, L, N_points)
    
    """Axial shape functions"""
    psi1 = 1 - x_e / L
    psi2 = x_e / L
    
    """Tansverse shape functions"""
    N1 = 1 - (3*x_e**2)/L**2 + (2*x_e**3)/L**3
    N2 = x_e**3/L**2 - (2*x_e**2)/L + x_e
    N3 = (3*x_e**2)/L**2 - (2*x_e**3)/L**3
    N4 =x_e**3/L**2 - x_e**2/L
    
    """Displacements"""
    u_x = psi1*d_e[0,0] + psi2*d_e[3,0]
    v_x = N1*d_e[1,0] + N2*d_e[2,0] + N3*d_e[4,0] + N4*d_e[5,0]
    
    """Delfections"""
    Deflections_XG = u_x*np.cos(alpha) - v_x*np.sin(alpha)
    Deflections_YG = u_x*np.sin(alpha) + v_x*np.cos(alpha)
    
    """Undeflected baseline of elements"""
    Undeflected_baseline_XG = np.linspace(node1XG, node2XG, N_points)
    Undeflected_baseline_YG = np.linspace(node1YG, node2YG, N_points)  
    
    """Deflected position of all points within element"""
    Deflected_XG = Undeflected_baseline_XG + Disp_mag*Deflections_XG
    Deflected_YG = Undeflected_baseline_YG + Disp_mag*Deflections_YG
    
    """Plot"""
    plt.plot(Undeflected_baseline_XG, Undeflected_baseline_YG, 'b.-')
    plt.plot(Deflected_XG,Deflected_YG, 'r.-')
    
Disp_mag = 10
N_points = 10
 
element_1 = plot_deflected_shape(0, 0, 0, 2.5, d1, Disp_mag, N_points, L1, theta_rad1)
element_2 = plot_deflected_shape(0, 2.5, 0, 5, d2, Disp_mag, N_points, L2, theta_rad2)
element_3 = plot_deflected_shape(0, 5, 0, 7.5, d3, Disp_mag, N_points, L3, theta_rad3)
element_4 = plot_deflected_shape(0, 2.5, 1, 0, d4, Disp_mag, N_points, L4, theta_rad4)
element_5 = plot_deflected_shape(0, 2.5, 1, 2.5, d5, Disp_mag, N_points, L5, theta_rad5)
element_6 = plot_deflected_shape(0, 2.5, 1, 5, d6, Disp_mag, N_points, L6, theta_rad6)
element_7 = plot_deflected_shape(0, 5, 1, 5, d7, Disp_mag, N_points, L7, theta_rad7)
element_8 = plot_deflected_shape(0, 7.5, 1, 5, d8, Disp_mag, N_points, L8, theta_rad8)
element_9 = plot_deflected_shape(1, 0, 1, 2.5, d9, Disp_mag, N_points, L9, theta_rad9)
element_10 = plot_deflected_shape(1, 2.5, 1, 5, d10, Disp_mag, N_points, L10, theta_rad10)
 
plt.grid(True)
plt.xlim(-2, 3)      # extend x-axis to 3
plt.ylim(0, 8)      # set y-axis limit to fit your structure
plt.show()
 
"""Deflection at the tip"""
print("Deflection in y and x direction at the tip")
print(f"x = {float(q[6,0]) * 1e3:.3f} mm, y = {float(q[7,0]) * 1e3:.3f} mm, rotation = {float(q[8,0]) * 1e3:.3f}mrad\n")
 
"""Reaction forces"""
# LEFT support at (0,0): only element 1, node i (local DOFs 0..2)
R_left_x = float(F1[0,0])
R_left_y = float(F1[1,0])
R_left_M = float(F1[2,0])
 
# RIGHT support at (1,0): element 4 node j (DOFs 3..5) + element 9 node i (DOFs 0..2)
R_right_x = float(F4[3,0] + F9[0,0])
R_right_y = float(F4[4,0] + F9[1,0])
R_right_M = float(F4[5,0] + F9[2,0])
 
print("Left support")
print(f"x = {R_left_x:.0f} N, y = {R_left_y:.0f} N, M = {R_left_M:.0f} Nm\n")
 
print("Right support")
print(f"x = {R_right_x:.0f} N, y = {R_right_y:.0f} N, M = {R_right_M:.0f} Nm\n")
 
 
print(f"The axial force of element 1: {float(abs(f1[0,0])):.3} N")
print(f"The axial force of element 2: {float(abs(f2[0,0])):.3} N")
print(f"The axial force of element 3: {float(abs(f3[0,0])):.3} N")
print(f"The axial force of element 4: {float(abs(f4[0,0])):.3} N")
print(f"The axial force of element 5: {float(abs(f5[0,0])):.3} N")
print(f"The axial force of element 6: {float(abs(f6[0,0])):.3} N")
print(f"The axial force of element 7: {float(abs(f7[0,0])):.3} N")
print(f"The axial force of element 8: {float(abs(f8[0,0])):.3} N")
print(f"The axial force of element 9: {float(abs(f9[0,0])):.3} N")
print(f"The axial force of element 10: {float(abs(f10[0,0])):.3} N\n")
 
 
"""Bending + Axial stress for element 1"""
c = 70e-3      # [m] distance to outer fibre   
 
Me1 = max(abs(float(f1[2,0])), abs(float(f1[5,0])))
Be1 = (Me1 * c) / I1
Af1 = abs(float(f1[0,0])) / A
sigma_tot1 = (Be1 + Af1) / 1e6
 
print("The axial + bending stress for each element")
print(f"Element 1 combined stress: {sigma_tot1:.3f} MPa")
 
 
"""Bending + Axial stress for element 2"""
Me2 = max(abs(float(f2[2,0])), abs(float(f2[5,0])))
Be2 = (Me2 * c) / I2
Af2 = abs(float(f2[0,0])) / A
sigma_tot2 = (Be2 + Af2) / 1e6
print(f"Element 2 combined stress: {sigma_tot2:.3f} MPa")
 
 
"""Bending + Axial stress for element 3"""
Me3 = max(abs(float(f3[2,0])), abs(float(f3[5,0])))
Be3 = (Me3 * c) / I3
Af3 = abs(float(f3[0,0])) / A
sigma_tot3 = (Be3 + Af3) / 1e6
print(f"Element 3 combined stress: {sigma_tot3:.3f} MPa")
 
 
"""Bending + Axial stress for element 4"""
Me4 = max(abs(float(f4[2,0])), abs(float(f4[5,0])))
Be4 = (Me4 * c) / I4
Af4 = abs(float(f4[0,0])) / A
sigma_tot4 = (Be4 + Af4) / 1e6
print(f"Element 4 combined stress: {sigma_tot4:.3f} MPa")
 
 
"""Bending + Axial stress for element 5"""
Me5 = max(abs(float(f5[2,0])), abs(float(f5[5,0])))
Be5 = (Me5 * c) / I5
Af5 = abs(float(f5[0,0])) / A
sigma_tot5 = (Be5 + Af5) / 1e6
print(f"Element 5 combined stress: {sigma_tot5:.3f} MPa")
 
 
"""Bending + Axial stress for element 6"""
Me6 = max(abs(float(f6[2,0])), abs(float(f6[5,0])))
Be6 = (Me6 * c) / I6
Af6 = abs(float(f6[0,0])) / A
sigma_tot6 = (Be6 + Af6) / 1e6
print(f"Element 6 combined stress: {sigma_tot6:.3f} MPa")
 
 
"""Bending + Axial stress for element 7"""
Me7 = max(abs(float(f7[2,0])), abs(float(f7[5,0])))
Be7 = (Me7 * c) / I7
Af7 = abs(float(f7[0,0])) / A
sigma_tot7 = (Be7 + Af7) / 1e6
print(f"Element 7 combined stress: {sigma_tot7:.3f} MPa")
 
 
"""Bending + Axial stress for element 8"""
Me8 = max(abs(float(f8[2,0])), abs(float(f8[5,0])))
Be8 = (Me8 * c) / I8
Af8 = abs(float(f8[0,0])) / A
sigma_tot8 = (Be8 + Af8) / 1e6
print(f"Element 8 combined stress: {sigma_tot8:.3f} MPa")
 
 
"""Bending + Axial stress for element 9"""
Me9 = max(abs(float(f9[2,0])), abs(float(f9[5,0])))
Be9 = (Me9 * c) / I9
Af9 = abs(float(f9[0,0])) / A
sigma_tot9 = (Be9 + Af9) / 1e6
print(f"Element 9 combined stress: {sigma_tot9:.3f} MPa")
 
 
"""Bending + Axial stress for element 10"""
Me10 = max(abs(float(f10[2,0])), abs(float(f10[5,0])))
Be10 = (Me10 * c) / I10
Af10 = abs(float(f10[0,0])) / A
sigma_tot10 = (Be10 + Af10) / 1e6
print(f"Element 10 combined stress: {sigma_tot10:.3f} MPa")
 
 
Bar code
 
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Aug 20 21:46:42 2025
 
@author: oliverchu
"""
 
 
import numpy as np
import matplotlib.pyplot as plt
 
 
def local_bar(E,A,L):
    """generate element stiffness matrix"""
    return (E*A)/L * np.array([[1, -1],[-1, 1]])
 
 
def global_bar(K, alpha):
    """generate element stiffness matirx in global coords"""
    c = np.cos(alpha)
    s = np.sin(alpha)
    lambda_matrix = np.array([[c,s,0,0],[0,0,c,s]])
    trans_lambda = np.transpose(lambda_matrix)
    return (trans_lambda @ K @ lambda_matrix), lambda_matrix
 
I = 13628778.666666666e-12
 
"""Element 1"""
L1, I1, A1, E1, alpha1 = 2.5, I, 4.496e-3, 200e9, np.radians(90)
K1 = local_bar(E1, A1, L1)
K1_hat, lambda1 = global_bar(K1, alpha1)
assem1 = np.array([
    [0,0,1,0],
    [0,0,0,1],
    [0,0,0,0],
    [0,0,0,0],
    [0,0,0,0],
    [0,0,0,0],
    [0,0,0,0],
    [0,0,0,0],
    [0,0,0,0],
    [0,0,0,0]])
KG1 = assem1 @ K1_hat @ np.transpose(assem1)
 
 
"""Element 2"""
L2, I2, A2, E2 ,alpha2 = 2.5, I, 4.496e-3, 200e9, np.radians(90)
K2 = local_bar(E2, A2, L2)
K2_hat, lambda2 = global_bar(K2, alpha2)
assem2 = np.array([
    [1,0,0,0],
    [0,1,0,0],
    [0,0,1,0],
    [0,0,0,1],
    [0,0,0,0],
    [0,0,0,0],
    [0,0,0,0],
    [0,0,0,0],
    [0,0,0,0],
    [0,0,0,0]])
KG2 = assem2 @ K2_hat @ np.transpose(assem2)
 
 
"""Element 3"""
L3, I3, A3, E3, alpha3 = 2.5, I, 4.496e-3, 200e9, np.radians(90)
K3 = local_bar(E3, A3, L3)
K3_hat, lambda3 = global_bar(K3, alpha3)
assem3 = np.array([
    [0,0,0,0],
    [0,0,0,0],
    [1,0,0,0],
    [0,1,0,0],
    [0,0,1,0],
    [0,0,0,1],
    [0,0,0,0],
    [0,0,0,0],
    [0,0,0,0],
    [0,0,0,0]])
KG3 = assem3 @ K3_hat @ np.transpose(assem3)
 
 
"""Element 4"""
L4, I4, A4, E4 ,alpha4 = np.sqrt(2.5**2+1**2), I, 4.496e-3, 200e9, np.radians(-68.2)
K4 = local_bar(E4, A4, L4)
K4_hat, lambda4 = global_bar(K4, alpha4)
assem4 = np.array([
    [1,0,0,0],
    [0,1,0,0],
    [0,0,0,0],
    [0,0,0,0],
    [0,0,0,0],
    [0,0,0,0],
    [0,0,0,0],
    [0,0,0,0],
    [0,0,0,0],
    [0,0,0,0]])
KG4 = assem4 @ K4_hat @ np.transpose(assem4)
 
 
"""Element 5"""
L5, I5, A5, E5, alpha5 = 1, I, 4.496e-3, 200e9, np.radians(0)
K5 = local_bar(E5, A5, L5)
K5_hat, lambda5 = global_bar(K5, alpha5)
assem5 = np.array([
    [1,0,0,0],
    [0,1,0,0],
    [0,0,0,0],
    [0,0,0,0],
    [0,0,0,0],
    [0,0,0,0],
    [0,0,0,0],
    [0,0,0,0],
    [0,0,1,0],
    [0,0,0,1]])
KG5 = assem5 @ K5_hat @ np.transpose(assem5)
 
 
"""Element 6"""
L6, I6, A6, E6 ,alpha6 = np.sqrt(2.5**2+1**2), I, 4.496e-3, 200e9, np.radians(68.2)
K6 = local_bar(E6, A6, L6)
K6_hat, lambda6 = global_bar(K6, alpha6)
assem6 = np.array([
    [1,0,0,0],
    [0,1,0,0],
    [0,0,0,0],
    [0,0,0,0],
    [0,0,0,0],
    [0,0,0,0],
    [0,0,1,0],
    [0,0,0,1],
    [0,0,0,0],
    [0,0,0,0]])
KG6 = assem6 @ K6_hat @ np.transpose(assem6)
 
 
"""Element 7"""
L7, I7, A7, E7, alpha7 = 1, I, 4.496e-3, 200e9, np.radians(0)
K7 = local_bar(E7, A7, L7)
K7_hat, lambda7 = global_bar(K7, alpha7)
assem7 = np.array([
    [0,0,0,0],
    [0,0,0,0],
    [1,0,0,0],
    [0,1,0,0],
    [0,0,0,0],
    [0,0,0,0],
    [0,0,1,0],
    [0,0,0,1],
    [0,0,0,0],
    [0,0,0,0]])
KG7 = assem7 @ K7_hat @ np.transpose(assem7)
 
 
"""Element 8"""
L8, I8, A8, E8 ,alpha8 = np.sqrt(2.5**2+1**2), I, 4.496e-3, 200e9, np.radians(-68.2)
K8 = local_bar(E8, A8, L8)
K8_hat, lambda8 = global_bar(K8, alpha8)
assem8 = np.array([
    [0,0,0,0],
    [0,0,0,0],
    [0,0,0,0],
    [0,0,0,0],
    [1,0,0,0],
    [0,1,0,0],
    [0,0,1,0],
    [0,0,0,1],
    [0,0,0,0],
    [0,0,0,0]])
KG8 = assem8 @ K8_hat @ np.transpose(assem8)
 
 
"""Element 9"""
L9, I9, A9, E9, alpha9 = 2.5, I, 4.496e-3, 200e9, np.radians(90)
K9 = local_bar(E9, A9, L9)
K9_hat, lambda9 = global_bar(K9, alpha9)
assem9 = np.array([
    [0,0,0,0],
    [0,0,0,0],
    [0,0,0,0],
    [0,0,0,0],
    [0,0,0,0],
    [0,0,0,0],
    [0,0,0,0],
    [0,0,0,0],
    [0,0,1,0],
    [0,0,0,1]])
KG9 = assem9 @ K9_hat @ np.transpose(assem9)
 
 
"""Element 10"""
L10, I10, A10, E10 ,alpha10 = 2.5, I, 4.496e-3, 200e9, np.radians(90)
K10 = local_bar(E10, A10, L10)
K10_hat, lambda10 = global_bar(K10, alpha10)
assem10 = np.array([
    [0,0,0,0],
    [0,0,0,0],
    [0,0,0,0],
    [0,0,0,0],
    [0,0,0,0],
    [0,0,0,0],
    [0,0,1,0],
    [0,0,0,1],
    [1,0,0,0],
    [0,1,0,0]])
KG10 = assem10 @ K10_hat @ np.transpose(assem10)
 
 
"""Total KG"""
KG = KG1 + KG2 + KG3 + KG4 +KG5 + KG6 + KG7 + KG8 + KG9 + KG10
 
 
"""Global forcing vector"""
Q = np.array([[15e3], [0], [30e3], [0], [15e3], [0], [0], [0], [0], [0]])
 
 
"""Q Total"""
Qtot = Q
 
 
"""Solving overall structure design"""
q = np.linalg.solve(KG, Qtot)
 
 
"""Deflections for each element"""
D1 = np.transpose(assem1) @ q
D2 = np.transpose(assem2) @ q
D3 = np.transpose(assem3) @ q
D4 = np.transpose(assem4) @ q
D5 = np.transpose(assem5) @ q
D6 = np.transpose(assem6) @ q
D7 = np.transpose(assem7) @ q
D8 = np.transpose(assem8) @ q
D9 = np.transpose(assem9) @ q
D10 = np.transpose(assem10) @ q
 
 
 
"""Element deflection vectors in global coords"""
d1 = lambda1 @ D1
d2 = lambda2 @ D2
d3 = lambda3 @ D3
d4 = lambda4 @ D4
d5 = lambda5 @ D5
d6 = lambda6 @ D6
d7 = lambda7 @ D7
d8 = lambda8 @ D8
d9 = lambda9 @ D9
d10 = lambda10 @ D10
 
 
"""Element nodal forcing vectors in local coords"""
f1 = K1 @ d1
f2 = K2 @ d2
f3 = K3 @ d3
f4 = K4 @ d4
f5 = K5 @ d5
f6 = K6 @ d6
f7 = K7 @ d7
f8 = K8 @ d8
f9 = K9 @ d9
f10 = K10 @ d10
 
 
"""Element nodal forcing vectors in global coords"""
F1 = K1_hat @ D1
F2 = K2_hat @ D2
F3 = K3_hat @ D3
F4 = K4_hat @ D4
F5 = K5_hat @ D5
F6 = K6_hat @ D6
F7 = K7_hat @ D7
F8 = K8_hat @ D8
F9 = K9_hat @ D9
F10 = K10_hat @ D10
 
 
def plot_truss_element(xi, yi, xj, yj, D_e, disp_mag=1.0, n_pts=10):
    """
    Plot undeflected and deflected shape of a TRUSS (bar) element.
    (xi, yi) -> (xj, yj) are node coordinates in metres.
    D_e is 4x1 or (4,) global DOFs [u_i, v_i, u_j, v_j].
    disp_mag is a visual scale factor only.
    """
    ui, vi, uj, vj = np.asarray(D_e).reshape(-1)  # [u_i, v_i, u_j, v_j]
    # baseline (straight line)
    X0 = np.linspace(xi, xj, n_pts)
    Y0 = np.linspace(yi, yj, n_pts)
    # linear interpolation of nodal displacements along the bar
    U = np.linspace(ui, uj, n_pts)
    V = np.linspace(vi, vj, n_pts)
    # deflected coordinates
    Xd = X0 + disp_mag * U
    Yd = Y0 + disp_mag * V
    # plot
    plt.plot(X0, Y0, 'b.-', linewidth=1, markersize=3)
    plt.plot(Xd, Yd, 'r.-', linewidth=1, markersize=3)
 
    
    
Disp_mag = 20
N_points = 10
 
plot_truss_element(0, 0,   0, 2.5, D1,  Disp_mag, N_points)
plot_truss_element(0, 2.5, 0, 5,   D2,  Disp_mag, N_points)
plot_truss_element(0, 5,   0, 7.5, D3,  Disp_mag, N_points)
plot_truss_element(0, 2.5, 1, 0,   D4,  Disp_mag, N_points)
plot_truss_element(0, 2.5, 1, 2.5, D5,  Disp_mag, N_points)
plot_truss_element(0, 2.5, 1, 5,   D6,  Disp_mag, N_points)
plot_truss_element(0, 5,   1, 5,   D7,  Disp_mag, N_points)
plot_truss_element(0, 7.5, 1, 5,   D8,  Disp_mag, N_points)
plot_truss_element(1, 0,   1, 2.5, D9,  Disp_mag, N_points)
plot_truss_element(1, 2.5, 1, 5,   D10, Disp_mag, N_points)
 
plt.grid(True)
plt.xlim(-2, 3)      # extend x-axis to 3
plt.ylim(0, 8)      # set y-axis limit to fit your structure
plt.show()
 
"""Deflection at the tip"""
print("Deflection in y and x direction at the tip")
print(f"x = {float(q[4,0]) * 1e3:.3f} mm, y = {float(q[5,0]) * 1e3:.3f} mm\n")
 
 
A = 4.496e-3
"""Reaction forces"""
# LEFT support at (0,0): only element 1, node i (local DOFs 0..2)
R_left_x = float(F1[0,0])
R_left_y = float(F1[1,0])
 
 
# RIGHT support at (1,0): element 4 node j (DOFs 3..5) + element 9 node i (DOFs 0..2)
R_right_x = float(F4[2,0] + F9[0,0])
R_right_y = float(F4[3,0] + F9[1,0])
 
 
print("Left support")
print(f"x = {R_left_x:.0f} N, y = {R_left_y:.0f} N\n")
 
print("Right support")
print(f"x = {R_right_x:.0f} N, y = {R_right_y:.0f} N")
 
 
print(f"The axial force of element 1: {float(abs(f1[0,0])):.3} N")
print(f"The axial force of element 2: {float(abs(f2[0,0])):.3} N")
print(f"The axial force of element 3: {float(abs(f3[0,0])):.3} N")
print(f"The axial force of element 4: {float(abs(f4[0,0])):.3} N")
print(f"The axial force of element 5: {float(abs(f5[0,0])):.3} N")
print(f"The axial force of element 6: {float(abs(f6[0,0])):.3} N")
print(f"The axial force of element 7: {float(abs(f7[0,0])):.3} N")
print(f"The axial force of element 8: {float(abs(f8[0,0])):.3} N")
print(f"The axial force of element 9: {float(abs(f9[0,0])):.3} N")
print(f"The axial force of element 10: {float(abs(f10[0,0])):.3} N\n")
 
 
"""Axial stress for element 1"""
Af1 = abs(float(f1[0,0])) / A
sigma_tot1 = (Af1) / 1e6
 
print("The axial for each element")
print(f"Element 1 combined stress: {sigma_tot1:.3f} MPa")
 
"""Axial stress for element 2"""
Af2 = abs(float(f2[0,0])) / A
sigma_tot2 = (Af2) / 1e6
print(f"Element 2 combined stress: {sigma_tot2:.3f} MPa")
 
"""Axial stress for element 3"""
Af3 = abs(float(f3[0,0])) / A
sigma_tot3 = (Af3) / 1e6
print(f"Element 3 combined stress: {sigma_tot3:.3f} MPa")
 
"""Axial stress for element 4"""
Af4 = abs(float(f4[0,0])) / A
sigma_tot4 = (Af4) / 1e6
print(f"Element 4 combined stress: {sigma_tot4:.3f} MPa")
 
"""Axial stress for element 5"""
Af5 = abs(float(f5[0,0])) / A
sigma_tot5 = (Af5) / 1e6
print(f"Element 5 combined stress: {sigma_tot5:.3f} MPa")
 
"""Axial stress for element 6"""
Af6 = abs(float(f6[0,0])) / A
sigma_tot6 = (Af6) / 1e6
print(f"Element 6 combined stress: {sigma_tot6:.3f} MPa")
 
"""Axial stress for element 7"""
Af7 = abs(float(f7[0,0])) / A
sigma_tot7 = (Af7) / 1e6
print(f"Element 7 combined stress: {sigma_tot7:.3f} MPa")
 
"""Axial stress for element 8"""
Af8 = abs(float(f8[0,0])) / A
sigma_tot8 = (Af8) / 1e6
print(f"Element 8 combined stress: {sigma_tot8:.3f} MPa")
 
"""Axial stress for element 9"""
Af9 = abs(float(f9[0,0])) / A
sigma_tot9 = (Af9) / 1e6
print(f"Element 9 combined stress: {sigma_tot9:.3f} MPa")
 
"""Axial stress for element 10"""
Af10 = abs(float(f10[0,0])) / A
sigma_tot10 = (Af10) / 1e6
print(f"Element 10 combined stress: {sigma_tot10:.3f} MPa")
 
 
# Low wind code 
 
"""
Created on Thu Jul 31 14:16:22 2025
 
@author: camsmith
"""
import numpy as np 
import matplotlib.pyplot as plt
 
def local_frame(E, A, L, I): #This matrix stays the same always
    stiff = (E * I) / L**3 #Sets the stiffness for the element
    beta = A * L**2 / I
    ke = stiff * np.array([
    [ beta,     0,       0,     -beta,     0,       0     ],
    [ 0,       12,     6*L,       0,     -12,     6*L    ],
    [ 0,      6*L,   4*L**2,      0,    -6*L,   2*L**2   ],
    [-beta,     0,       0,      beta,     0,       0     ],
    [ 0,     -12,    -6*L,        0,      12,    -6*L    ],
    [ 0,      6*L,   2*L**2,      0,    -6*L,   4*L**2   ]
])
    return ke
    
def upplambda(theta): #Uppercase lambda converts from local to global coord system
    c = np.cos(theta)
    s = np.sin(theta)
    lam = np.array([[c, s, 0], [-s, c, 0], [0, 0, 1]])
    upplambda = np.block([
        [lam,               np.zeros((3,3))],
        [np.zeros((3,3)),   lam]
    ])
    
    return upplambda
 
def global_frame(ke_local, upplambda):
    return upplambda.T @ ke_local @ upplambda #Returning the stiffness matrix in global coord
 
"""Element 1"""
E1, A1, L1, I1 = 200e9, 4.496e-3, 2.5, 13628778.666666666e-12 #Material Props
theta_deg1 = 90 #Angle roatated
theta_rad1 = np.radians(theta_deg1)
ke1 = local_frame(E1, A1, L1, I1)
khat1 = global_frame(ke1, upplambda(theta_rad1))
             #D#1  2  3  4  5  6
A1 = np.array([[0, 0, 0, 1, 0, 0], #q1
               [0, 0, 0, 0, 1, 0], #q2
               [0, 0, 0, 0, 0, 1], #q3
               [0, 0, 0, 0, 0, 0], #q4
               [0, 0, 0, 0, 0, 0], #q5
               [0, 0, 0, 0, 0, 0], #q6
               [0, 0, 0, 0, 0, 0], #q7
               [0, 0, 0, 0, 0, 0], #q8
               [0, 0, 0, 0, 0, 0], #q9
               [0, 0, 0, 0, 0, 0], #q10
               [0, 0, 0, 0, 0, 0], #q11
               [0, 0, 0, 0, 0, 0], #q12
               [0, 0, 0, 0, 0, 0], #q13
               [0, 0, 0, 0, 0, 0], #q14
               [0, 0, 0, 0, 0, 0]]) #q15
         
"""Element 2"""
E2, A2, L2, I2 = 200e9, 4.496e-3, 2.5, 13628778.666666666e-12 #Material Props
theta_deg2 = 90 #Angle roatated
theta_rad2 = np.radians(theta_deg2)
ke2 = local_frame(E2, A2, L2, I2)
khat2 = global_frame(ke2, upplambda(theta_rad2))
             #D#1  2  3  4  5  6
A2 = np.array([[1, 0, 0, 0, 0, 0], #q1
               [0, 1, 0, 0, 0, 0], #q2
               [0, 0, 1, 0, 0, 0], #q3
               [0, 0, 0, 1, 0, 0], #q4
               [0, 0, 0, 0, 1, 0], #q5
               [0, 0, 0, 0, 0, 1], #q6
               [0, 0, 0, 0, 0, 0], #q7
               [0, 0, 0, 0, 0, 0], #q8
               [0, 0, 0, 0, 0, 0], #q9
               [0, 0, 0, 0, 0, 0], #q10
               [0, 0, 0, 0, 0, 0], #q11
               [0, 0, 0, 0, 0, 0], #q12
               [0, 0, 0, 0, 0, 0], #q13
               [0, 0, 0, 0, 0, 0], #q14
               [0, 0, 0, 0, 0, 0]]) #q15
         
               
"""Element 3"""
E3, A3, L3, I3 = 200e9, 4.496e-3, 2.5, 13628778.666666666e-12 #Material Props
theta_deg3 = 90 #Angle roatated
theta_rad3 = np.radians(theta_deg3)
ke3 = local_frame(E3, A3, L3, I3)
khat3 = global_frame(ke3, upplambda(theta_rad3))
             #D#1  2  3  4  5  6
A3 = np.array([[0, 0, 0, 0, 0, 0], #q1
               [0, 0, 0, 0, 0, 0], #q2
               [0, 0, 0, 0, 0, 0], #q3
               [1, 0, 0, 0, 0, 0], #q4
               [0, 1, 0, 0, 0, 0], #q5
               [0, 0, 1, 0, 0, 0], #q6
               [0, 0, 0, 1, 0, 0], #q7
               [0, 0, 0, 0, 1, 0], #q8
               [0, 0, 0, 0, 0, 1], #q9
               [0, 0, 0, 0, 0, 0], #q10
               [0, 0, 0, 0, 0, 0], #q11
               [0, 0, 0, 0, 0, 0], #q12
               [0, 0, 0, 0, 0, 0], #q13
               [0, 0, 0, 0, 0, 0], #q14
               [0, 0, 0, 0, 0, 0]]) #q15
         
 
"""Element 4"""
E4, A4, L4, I4 = 200e9, 4.496e-3, np.sqrt(2.5**2 + 1**2), 13628778.666666666e-12 #Material Props
theta_deg4 = -68.2 #Angle roatated
theta_rad4 = np.radians(theta_deg4)
ke4 = local_frame(E4, A4, L4, I4)
khat4 = global_frame(ke4, upplambda(theta_rad4))
             #D#1  2  3  4  5  6
A4 = np.array([[1, 0, 0, 0, 0, 0], #q1
               [0, 1, 0, 0, 0, 0], #q2
               [0, 0, 1, 0, 0, 0], #q3
               [0, 0, 0, 0, 0, 0], #q4
               [0, 0, 0, 0, 0, 0], #q5
               [0, 0, 0, 0, 0, 0], #q6
               [0, 0, 0, 0, 0, 0], #q7
               [0, 0, 0, 0, 0, 0], #q8
               [0, 0, 0, 0, 0, 0], #q9
               [0, 0, 0, 0, 0, 0], #q10
               [0, 0, 0, 0, 0, 0], #q11
               [0, 0, 0, 0, 0, 0], #q12
               [0, 0, 0, 0, 0, 0], #q13
               [0, 0, 0, 0, 0, 0], #q14
               [0, 0, 0, 0, 0, 0]]) #q15
         
 
"""Element 5"""
E5, A5, L5, I5 = 200e9, 4.496e-3, 1, 13628778.666666666e-12 #Material Props
theta_deg5 = 0 #Angle roatated
theta_rad5 = np.radians(theta_deg5)
ke5 = local_frame(E5, A5, L5, I5)
khat5 = global_frame(ke5, upplambda(theta_rad5))
             #D#1  2  3  4  5  6
A5 = np.array([[1, 0, 0, 0, 0, 0], #q1
               [0, 1, 0, 0, 0, 0], #q2
               [0, 0, 1, 0, 0, 0], #q3
               [0, 0, 0, 0, 0, 0], #q4
               [0, 0, 0, 0, 0, 0], #q5
               [0, 0, 0, 0, 0, 0], #q6
               [0, 0, 0, 0, 0, 0], #q7
               [0, 0, 0, 0, 0, 0], #q8
               [0, 0, 0, 0, 0, 0], #q9
               [0, 0, 0, 0, 0, 0], #q10
               [0, 0, 0, 0, 0, 0], #q11
               [0, 0, 0, 0, 0, 0], #q12
               [0, 0, 0, 1, 0, 0], #q13
               [0, 0, 0, 0, 1, 0], #q14
               [0, 0, 0, 0, 0, 1]]) #q15
         
 
"""Element 6"""
E6, A6, L6, I6 = 200e9, 4.496e-3, np.sqrt(2.5**2 + 1**2), 13628778.666666666e-12 #Material Props
theta_deg6 = 68.2 #Angle roatated
theta_rad6 = np.radians(theta_deg6)
ke6 = local_frame(E6, A6, L6, I6)
khat6 = global_frame(ke6, upplambda(theta_rad6))
             #D#1  2  3  4  5  6
A6 = np.array([[1, 0, 0, 0, 0, 0], #q1
               [0, 1, 0, 0, 0, 0], #q2
               [0, 0, 1, 0, 0, 0], #q3
               [0, 0, 0, 0, 0, 0], #q4
               [0, 0, 0, 0, 0, 0], #q5
               [0, 0, 0, 0, 0, 0], #q6
               [0, 0, 0, 0, 0, 0], #q7
               [0, 0, 0, 0, 0, 0], #q8
               [0, 0, 0, 0, 0, 0], #q9
               [0, 0, 0, 1, 0, 0], #q10
               [0, 0, 0, 0, 1, 0], #q11
               [0, 0, 0, 0, 0, 1], #q12
               [0, 0, 0, 0, 0, 0], #q13
               [0, 0, 0, 0, 0, 0], #q14
               [0, 0, 0, 0, 0, 0]]) #q15
         
 
"""Element 7"""
E7, A7, L7, I7 = 200e9, 4.496e-3, 1, 13628778.666666666e-12 #Material Props
theta_deg7 = 0 #Angle roatated
theta_rad7 = np.radians(theta_deg7)
ke7 = local_frame(E7, A7, L7, I7)
khat7 = global_frame(ke7, upplambda(theta_rad7))
             #D#1  2  3  4  5  6
A7 = np.array([[0, 0, 0, 0, 0, 0], #q1
               [0, 0, 0, 0, 0, 0], #q2
               [0, 0, 0, 0, 0, 0], #q3
               [1, 0, 0, 0, 0, 0], #q4
               [0, 1, 0, 0, 0, 0], #q5
               [0, 0, 1, 0, 0, 0], #q6
               [0, 0, 0, 0, 0, 0], #q7
               [0, 0, 0, 0, 0, 0], #q8
               [0, 0, 0, 0, 0, 0], #q9
               [0, 0, 0, 1, 0, 0], #q10
               [0, 0, 0, 0, 1, 0], #q11
               [0, 0, 0, 0, 0, 1], #q12
               [0, 0, 0, 0, 0, 0], #q13
               [0, 0, 0, 0, 0, 0], #q14
               [0, 0, 0, 0, 0, 0]]) #q15
         
 
"""Element 8"""
E8, A8, L8, I8 = 200e9, 4.496e-3, 3, 13628778.666666666e-12#Material Props
theta_deg8 = -68.2 #Angle roatated
theta_rad8 = np.radians(theta_deg8)
ke8 = local_frame(E8, A8, L8, I8)
khat8 = global_frame(ke8, upplambda(theta_rad8))
             #D#1  2  3  4  5  6
A8 = np.array([[0, 0, 0, 0, 0, 0], #q1
               [0, 0, 0, 0, 0, 0], #q2
               [0, 0, 0, 0, 0, 0], #q3
               [0, 0, 0, 0, 0, 0], #q4
               [0, 0, 0, 0, 0, 0], #q5
               [0, 0, 0, 0, 0, 0], #q6
               [1, 0, 0, 0, 0, 0], #q7
               [0, 1, 0, 0, 0, 0], #q8
               [0, 0, 1, 0, 0, 0], #q9
               [0, 0, 0, 1, 0, 0], #q10
               [0, 0, 0, 0, 1, 0], #q11
               [0, 0, 0, 0, 0, 1], #q12
               [0, 0, 0, 0, 0, 0], #q13
               [0, 0, 0, 0, 0, 0], #q14
               [0, 0, 0, 0, 0, 0]]) #q15
         
 
"""Element 9"""
E9, A9, L9, I9 = 200e9, 4.496e-3, 2.5, 13628778.666666666e-12 #Material Props
theta_deg9 = 90 #Angle roatated
theta_rad9 = np.radians(theta_deg9)
ke9 = local_frame(E9, A9, L9, I9)
khat9 = global_frame(ke9, upplambda(theta_rad9))
             #D#1  2  3  4  5  6
A9 = np.array([[0, 0, 0, 0, 0, 0], #q1
               [0, 0, 0, 0, 0, 0], #q2
               [0, 0, 0, 0, 0, 0], #q3
               [0, 0, 0, 0, 0, 0], #q4
               [0, 0, 0, 0, 0, 0], #q5
               [0, 0, 0, 0, 0, 0], #q6
               [0, 0, 0, 0, 0, 0], #q7
               [0, 0, 0, 0, 0, 0], #q8
               [0, 0, 0, 0, 0, 0], #q9
               [0, 0, 0, 0, 0, 0], #q10
               [0, 0, 0, 0, 0, 0], #q11
               [0, 0, 0, 0, 0, 0], #q12
               [0, 0, 0, 1, 0, 0], #q13
               [0, 0, 0, 0, 1, 0], #q14
               [0, 0, 0, 0, 0, 1]]) #q15
         
 
"""Element 10"""
E10, A10, L10, I10 = 200e9, 4.496e-3, 2.5, 13628778.666666666e-12 #Material Props
theta_deg10 = 90 #Angle roatated
theta_rad10 = np.radians(theta_deg10)
ke10 = local_frame(E10, A10, L10, I10)
khat10 = global_frame(ke10, upplambda(theta_rad10))
             #D#1  2  3  4  5  6
A10 = np.array([[0, 0, 0, 0, 0, 0], #q1
                [0, 0, 0, 0, 0, 0], #q2
                [0, 0, 0, 0, 0, 0], #q3
                [0, 0, 0, 0, 0, 0], #q4
                [0, 0, 0, 0, 0, 0], #q5
                [0, 0, 0, 0, 0, 0], #q6
                [0, 0, 0, 0, 0, 0], #q7
                [0, 0, 0, 0, 0, 0], #q8
                [0, 0, 0, 0, 0, 0], #q9
                [0, 0, 0, 1, 0, 0], #q10
                [0, 0, 0, 0, 1, 0], #q11
                [0, 0, 0, 0, 0, 1], #q12
                [1, 0, 0, 0, 0, 0], #q13
                [0, 1, 0, 0, 0, 0], #q14
                [0, 0, 1, 0, 0, 0]]) #q15
 
"""Indivdual KG Matrixs"""
KG1 = A1 @ khat1 @ A1.T
KG2 = A2 @ khat2 @ A2.T
KG3 = A3 @ khat3 @ A3.T
KG4 = A4 @ khat4 @ A4.T
KG5 = A5 @ khat5 @ A5.T
KG6 = A6 @ khat6 @ A6.T
KG7 = A7 @ khat7 @ A7.T
KG8 = A8 @ khat8 @ A8.T
KG9 = A9 @ khat9 @ A9.T
KG10 = A10 @ khat10 @ A10.T
 
"""Total"""
KG = KG1 + KG2 + KG3 + KG4 + KG5 + KG6 + KG7 + KG8 + KG9 + KG10
 
 
"""Q Loads"""
low = -4880
med = -6560
high = -9280
vhigh = -12000
maxwind = -2349.001858441632 * 8
load = maxwind
 
 
"""Uniform Distrubuted Load"""
Ludl1 = L3
wdl1 = load
thetaudl1 = theta_rad3
Audl1 = A3
 
f_udl1 = np.array([[0], 
                  [(wdl1*Ludl1)/2], 
                  [(wdl1*Ludl1**2)/12], 
                  [0], 
                  [(wdl1*Ludl1)/2], 
                  [(-wdl1*Ludl1**2)/12]])
 
F_udl1 = upplambda(thetaudl1).T @ f_udl1
Q_udl1 = Audl1 @ F_udl1
 
"""Uniform Distrubuted Load"""
Ludl = L2
wdl = load
thetaudl = theta_rad2
Audl = A2
 
f_udl = np.array([[0], 
                  [(wdl*Ludl)/2], 
                  [(wdl*Ludl**2)/12], 
                  [0], 
                  [(wdl*Ludl)/2], 
                  [(-wdl*Ludl**2)/12]])
 
F_udl = upplambda(thetaudl).T @ f_udl
Q_udl = Audl @ F_udl
 
"""Q loads from structure"""
Qtot = Q_udl + Q_udl1
q = np.linalg.solve(KG,  Qtot)
 
"""Forces"""
F1 = khat1 @ A1.T @ q
F2 = khat2 @ A2.T @ q
F3 = khat3 @ A3.T @ q
F4 = khat4 @ A4.T @ q
F5 = khat5 @ A5.T @ q
F6 = khat6 @ A6.T @ q
F7 = khat7 @ A7.T @ q
F8 = khat8 @ A8.T @ q
F9 = khat9 @ A9.T @ q
F10 = khat10 @ A10.T @ q
 
"""Deflections"""
d1 = upplambda(theta_rad1) @ A1.T @ q
d2 = upplambda(theta_rad2) @ A2.T @ q
d3 = upplambda(theta_rad3) @ A3.T @ q
d4 = upplambda(theta_rad4) @ A4.T @ q
d5 = upplambda(theta_rad5) @ A5.T @ q
d6 = upplambda(theta_rad6) @ A6.T @ q
d7 = upplambda(theta_rad7) @ A7.T @ q
d8 = upplambda(theta_rad8) @ A8.T @ q
d9 = upplambda(theta_rad9) @ A9.T @ q
d10 = upplambda(theta_rad10) @ A10.T @ q
 
"""Internal forces"""
f1 = ke1 @ d1
f2 = ke2 @ d2
f3 = ke3 @ d3
f4 = ke4 @ d4
f5 = ke5 @ d5
f6 = ke6 @ d6
f7 = ke7 @ d7
f8 = ke8 @ d8
f9 = ke9 @ d9
f10 = ke10 @ d10
 
"""Tricky plot"""
def plot_deflected_shape(node1XG, node1YG, node2XG, node2YG, d_e, Disp_mag, N_points, L, alpha):
    """plot full deflected shape"""
    x_e = np.linspace(0, L, N_points)
    
    """Axial shape functions"""
    psi1 = 1 - x_e / L
    psi2 = x_e / L
    
    """Tansverse shape functions"""
    N1 = 1 - (3*x_e**2)/L**2 + (2*x_e**3)/L**3
    N2 = x_e**3/L**2 - (2*x_e**2)/L + x_e
    N3 = (3*x_e**2)/L**2 - (2*x_e**3)/L**3
    N4 =x_e**3/L**2 - x_e**2/L
    
    """Displacements"""
    u_x = psi1*d_e[0,0] + psi2*d_e[3,0]
    v_x = N1*d_e[1,0] + N2*d_e[2,0] + N3*d_e[4,0] + N4*d_e[5,0]
    
    """Delfections"""
    Deflections_XG = u_x*np.cos(alpha) - v_x*np.sin(alpha)
    Deflections_YG = u_x*np.sin(alpha) + v_x*np.cos(alpha)
    
    """Undeflected baseline of elements"""
    Undeflected_baseline_XG = np.linspace(node1XG, node2XG, N_points)
    Undeflected_baseline_YG = np.linspace(node1YG, node2YG, N_points)  
    
    """Deflected position of all points within element"""
    Deflected_XG = Undeflected_baseline_XG + Disp_mag*Deflections_XG
    Deflected_YG = Undeflected_baseline_YG + Disp_mag*Deflections_YG
    
    """Plot"""
    plt.plot(Undeflected_baseline_XG, Undeflected_baseline_YG, 'b.-')
    plt.plot(Deflected_XG,Deflected_YG, 'r.-')
    
Disp_mag = 10
N_points = 10
 
element_1 = plot_deflected_shape(0, 0, 0, 2.5, d1, Disp_mag, N_points, L1, theta_rad1)
element_2 = plot_deflected_shape(0, 2.5, 0, 5, d2, Disp_mag, N_points, L2, theta_rad2)
element_3 = plot_deflected_shape(0, 5, 0, 7.5, d3, Disp_mag, N_points, L3, theta_rad3)
element_4 = plot_deflected_shape(0, 2.5, 1, 0, d4, Disp_mag, N_points, L4, theta_rad4)
element_5 = plot_deflected_shape(0, 2.5, 1, 2.5, d5, Disp_mag, N_points, L5, theta_rad5)
element_6 = plot_deflected_shape(0, 2.5, 1, 5, d6, Disp_mag, N_points, L6, theta_rad6)
element_7 = plot_deflected_shape(0, 5, 1, 5, d7, Disp_mag, N_points, L7, theta_rad7)
element_8 = plot_deflected_shape(0, 7.5, 1, 5, d8, Disp_mag, N_points, L8, theta_rad8)
element_9 = plot_deflected_shape(1, 0, 1, 2.5, d9, Disp_mag, N_points, L9, theta_rad9)
element_10 = plot_deflected_shape(1, 2.5, 1, 5, d10, Disp_mag, N_points, L10, theta_rad10)
 
plt.grid(True)
plt.xlim(-2, 3)      # extend x-axis to 3
plt.ylim(0, 8)      # set y-axis limit to fit your structure
plt.show()
 
 
 
"""Deflection at the tip"""
print("Deflection in y and x direction at the tip")
print(f"x = {float(q[6,0]) * 1e3:.3f} mm, y = {float(q[7,0]) * 1e3:.3f} mm, rotation = {float(q[8,0]) * 1e3:.3f}mrad\n")
 
 
# --- Reactions at supports from element end-actions (no loads applied at base DOFs) ---
 
# Left support at (0,0): element 1, i-end (local DOFs 0..2)
R_left_x = float(F1[0,0])
R_left_y = float(F1[1,0])
R_left_M = float(F1[2,0])
 
# Right support at (1,0): element 4, j-end (local DOFs 3..5) + element 9, i-end (0..2)
R_right_x = (float(F4[3,0]) + float(F9[0,0]))
R_right_y = (float(F4[4,0]) + float(F9[1,0]))
R_right_M = (float(F4[5,0]) + float(F9[2,0]))
 
print("Left support")
print(f"x = {R_left_x:.0f} N, y = {R_left_y:.0f} N, M = {R_left_M:.0f} Nm\n")
print("Right support")
print(f"x = {R_right_x:.0f} N, y = {R_right_y:.0f} N, M = {R_right_M:.0f} Nm\n")
 
Fx_applied = float((Q_udl + Q_udl1).sum())  # sum of all global DOF x-entries; or sum only x-DOFs if you prefer
Fx_react   = R_left_x + R_right_x
print(f"ΣFx applied = {Fx_applied:.1f} N,  ΣRx = {Fx_react:.1f} N,  Σ = {Fx_applied + Fx_react:.1f} N\n")
 
 
 
print(f"The axial force of element 1: {float(abs(f1[0,0])):.3} N")
print(f"The axial force of element 2: {float(abs(f2[0,0])):.3} N")
print(f"The axial force of element 3: {float(abs(f3[0,0])):.3} N")
print(f"The axial force of element 4: {float(abs(f4[0,0])):.3} N")
print(f"The axial force of element 5: {float(abs(f5[0,0])):.3} N")
print(f"The axial force of element 6: {float(abs(f6[0,0])):.3} N")
print(f"The axial force of element 7: {float(abs(f7[0,0])):.3} N")
print(f"The axial force of element 8: {float(abs(f8[0,0])):.3} N")
print(f"The axial force of element 9: {float(abs(f9[0,0])):.3} N")
print(f"The axial force of element 10: {float(abs(f10[0,0])):.3} N\n")
 
 
"""Bending + Axial stress for element 1"""
c = 70e-3      # [m] distance to outer fibre
A = 4.256e-3   # [m^2] cross-sectional area
 
Me1 = max(abs(float(f1[2,0])), abs(float(f1[5,0])))
Be1 = (Me1 * c) / I1
Af1 = abs(float(f1[0,0])) / A
sigma_tot1 = (Be1 + Af1) / 1e6
 
print("The axial + bending stress for each element")
print(f"Element 1 combined stress: {sigma_tot1:.3f} MPa")
 
 
"""Bending + Axial stress for element 2"""
Me2 = max(abs(float(f2[2,0])), abs(float(f2[5,0])))
Be2 = (Me2 * c) / I2
Af2 = abs(float(f2[0,0])) / A
sigma_tot2 = (Be2 + Af2) / 1e6
print(f"Element 2 combined stress: {sigma_tot2:.3f} MPa")
 
 
"""Bending + Axial stress for element 3"""
Me3 = max(abs(float(f3[2,0])), abs(float(f3[5,0])))
Be3 = (Me3 * c) / I3
Af3 = abs(float(f3[0,0])) / A
sigma_tot3 = (Be3 + Af3) / 1e6
print(f"Element 3 combined stress: {sigma_tot3:.3f} MPa")
 
 
"""Bending + Axial stress for element 4"""
Me4 = max(abs(float(f4[2,0])), abs(float(f4[5,0])))
Be4 = (Me4 * c) / I4
Af4 = abs(float(f4[0,0])) / A
sigma_tot4 = (Be4 + Af4) / 1e6
print(f"Element 4 combined stress: {sigma_tot4:.3f} MPa")
 
 
"""Bending + Axial stress for element 5"""
Me5 = max(abs(float(f5[2,0])), abs(float(f5[5,0])))
Be5 = (Me5 * c) / I5
Af5 = abs(float(f5[0,0])) / A
sigma_tot5 = (Be5 + Af5) / 1e6
print(f"Element 5 combined stress: {sigma_tot5:.3f} MPa")
 
 
"""Bending + Axial stress for element 6"""
Me6 = max(abs(float(f6[2,0])), abs(float(f6[5,0])))
Be6 = (Me6 * c) / I6
Af6 = abs(float(f6[0,0])) / A
sigma_tot6 = (Be6 + Af6) / 1e6
print(f"Element 6 combined stress: {sigma_tot6:.3f} MPa")
 
 
"""Bending + Axial stress for element 7"""
Me7 = max(abs(float(f7[2,0])), abs(float(f7[5,0])))
Be7 = (Me7 * c) / I7
Af7 = abs(float(f7[0,0])) / A
sigma_tot7 = (Be7 + Af7) / 1e6
print(f"Element 7 combined stress: {sigma_tot7:.3f} MPa")
 
 
"""Bending + Axial stress for element 8"""
Me8 = max(abs(float(f8[2,0])), abs(float(f8[5,0])))
Be8 = (Me8 * c) / I8
Af8 = abs(float(f8[0,0])) / A
sigma_tot8 = (Be8 + Af8) / 1e6
print(f"Element 8 combined stress: {sigma_tot8:.3f} MPa")
 
 
"""Bending + Axial stress for element 9"""
Me9 = max(abs(float(f9[2,0])), abs(float(f9[5,0])))
Be9 = (Me9 * c) / I9
Af9 = abs(float(f9[0,0])) / A
sigma_tot9 = (Be9 + Af9) / 1e6
print(f"Element 9 combined stress: {sigma_tot9:.3f} MPa")
 
 
"""Bending + Axial stress for element 10"""
Me10 = max(abs(float(f10[2,0])), abs(float(f10[5,0])))
Be10 = (Me10 * c) / I10
Af10 = abs(float(f10[0,0])) / A
sigma_tot10 = (Be10 + Af10) / 1e6
print(f"Element 10 combined stress: {sigma_tot10:.3f} MPa\n")
 
 
"""Fos stuff"""
### The load has to be set to -8
FOS = 2.5
yeild = 350
max_allowed = yeild / FOS
p_max = max_allowed / sigma_tot1
wind_speed = np.sqrt(p_max / 0.6)
print("The maximum wind speed with a FOS = 2.5 is:")
print(f"{wind_speed:.2f}m/s")
print(p_max)
 
# Different stiffness matrix windy code
 
"""
Created on Thu Jul 31 14:16:22 2025
 
@author: camsmith
"""
import numpy as np 
import matplotlib.pyplot as plt
 
def local_frame(E, A, L, I, nu, As):
    # Shear modulus from Poisson's ratio
    G = 70e9
 
    
    # Timoshenko parameters
    stiff = (E * I) / (L**3)
    beta  = A * L**2 / I
    Phi   = (12.0 * E * I) / (G * As * L**2)
    inv   = 1.0 / (1.0 + Phi)
 
    ke = stiff * np.array([
        [  beta,          0.0,                   0.0,     -beta,          0.0,                    0.0 ],
        [   0.0,     12.0*inv,             6.0*L*inv,       0.0,    -12.0*inv,              6.0*L*inv ],
        [   0.0,     6.0*L*inv,  (4.0 + Phi)*L**2*inv,       0.0,    -6.0*L*inv,    (2.0 - Phi)*L**2*inv ],
        [ -beta,          0.0,                   0.0,      beta,          0.0,                    0.0 ],
        [   0.0,    -12.0*inv,            -6.0*L*inv,       0.0,     12.0*inv,             -6.0*L*inv ],
        [   0.0,     6.0*L*inv,  (2.0 - Phi)*L**2*inv,       0.0,    -6.0*L*inv,    (4.0 + Phi)*L**2*inv ],
    ])
    return ke
    
def upplambda(theta): #Uppercase lambda converts from local to global coord system
    c = np.cos(theta)
    s = np.sin(theta)
    lam = np.array([[c, s, 0], [-s, c, 0], [0, 0, 1]])
    upplambda = np.block([
        [lam,               np.zeros((3,3))],
        [np.zeros((3,3)),   lam]
    ])
    
    return upplambda
 
def global_frame(ke_local, upplambda):
    return upplambda.T @ ke_local @ upplambda #Returning the stiffness matrix in global coord
 
"""Global"""
As = 140e-3 * 12e-3
nu = 0.3
A = 4.496e-3
I = 13628778.666666666e-12
"""Element 1"""
E1, A1, L1, I1 = 200e9, A, 2.5, I #Material Props
theta_deg1 = 90 #Angle roatated
theta_rad1 = np.radians(theta_deg1)
ke1 = local_frame(E1, A1, L1, I1, nu, As)
khat1 = global_frame(ke1, upplambda(theta_rad1))
             #D#1  2  3  4  5  6
A1 = np.array([[0, 0, 0, 1, 0, 0], #q1
               [0, 0, 0, 0, 1, 0], #q2
               [0, 0, 0, 0, 0, 1], #q3
               [0, 0, 0, 0, 0, 0], #q4
               [0, 0, 0, 0, 0, 0], #q5
               [0, 0, 0, 0, 0, 0], #q6
               [0, 0, 0, 0, 0, 0], #q7
               [0, 0, 0, 0, 0, 0], #q8
               [0, 0, 0, 0, 0, 0], #q9
               [0, 0, 0, 0, 0, 0], #q10
               [0, 0, 0, 0, 0, 0], #q11
               [0, 0, 0, 0, 0, 0], #q12
               [0, 0, 0, 0, 0, 0], #q13
               [0, 0, 0, 0, 0, 0], #q14
               [0, 0, 0, 0, 0, 0]]) #q15
         
"""Element 2"""
E2, A2, L2, I2 = 200e9, A, 2.5, I #Material Props
theta_deg2 = 90 #Angle roatated
theta_rad2 = np.radians(theta_deg2)
ke2 = local_frame(E2, A2, L2, I2, nu, As)
khat2 = global_frame(ke2, upplambda(theta_rad2))
             #D#1  2  3  4  5  6
A2 = np.array([[1, 0, 0, 0, 0, 0], #q1
               [0, 1, 0, 0, 0, 0], #q2
               [0, 0, 1, 0, 0, 0], #q3
               [0, 0, 0, 1, 0, 0], #q4
               [0, 0, 0, 0, 1, 0], #q5
               [0, 0, 0, 0, 0, 1], #q6
               [0, 0, 0, 0, 0, 0], #q7
               [0, 0, 0, 0, 0, 0], #q8
               [0, 0, 0, 0, 0, 0], #q9
               [0, 0, 0, 0, 0, 0], #q10
               [0, 0, 0, 0, 0, 0], #q11
               [0, 0, 0, 0, 0, 0], #q12
               [0, 0, 0, 0, 0, 0], #q13
               [0, 0, 0, 0, 0, 0], #q14
               [0, 0, 0, 0, 0, 0]]) #q15
         
               
"""Element 3"""
E3, A3, L3, I3 = 200e9, A, 2.5, I #Material Props
theta_deg3 = 90 #Angle roatated
theta_rad3 = np.radians(theta_deg3)
ke3 = local_frame(E3, A3, L3, I3, nu, As)
khat3 = global_frame(ke3, upplambda(theta_rad3))
             #D#1  2  3  4  5  6
A3 = np.array([[0, 0, 0, 0, 0, 0], #q1
               [0, 0, 0, 0, 0, 0], #q2
               [0, 0, 0, 0, 0, 0], #q3
               [1, 0, 0, 0, 0, 0], #q4
               [0, 1, 0, 0, 0, 0], #q5
               [0, 0, 1, 0, 0, 0], #q6
               [0, 0, 0, 1, 0, 0], #q7
               [0, 0, 0, 0, 1, 0], #q8
               [0, 0, 0, 0, 0, 1], #q9
               [0, 0, 0, 0, 0, 0], #q10
               [0, 0, 0, 0, 0, 0], #q11
               [0, 0, 0, 0, 0, 0], #q12
               [0, 0, 0, 0, 0, 0], #q13
               [0, 0, 0, 0, 0, 0], #q14
               [0, 0, 0, 0, 0, 0]]) #q15
         
 
"""Element 4"""
E4, A4, L4, I4 = 200e9, A, 2.6925, I #Material Props
theta_deg4 = -68.2 #Angle roatated
theta_rad4 = np.radians(theta_deg4)
ke4 = local_frame(E4, A4, L4, I4, nu, As)
khat4 = global_frame(ke4, upplambda(theta_rad4))
             #D#1  2  3  4  5  6
A4 = np.array([[1, 0, 0, 0, 0, 0], #q1
               [0, 1, 0, 0, 0, 0], #q2
               [0, 0, 1, 0, 0, 0], #q3
               [0, 0, 0, 0, 0, 0], #q4
               [0, 0, 0, 0, 0, 0], #q5
               [0, 0, 0, 0, 0, 0], #q6
               [0, 0, 0, 0, 0, 0], #q7
               [0, 0, 0, 0, 0, 0], #q8
               [0, 0, 0, 0, 0, 0], #q9
               [0, 0, 0, 0, 0, 0], #q10
               [0, 0, 0, 0, 0, 0], #q11
               [0, 0, 0, 0, 0, 0], #q12
               [0, 0, 0, 0, 0, 0], #q13
               [0, 0, 0, 0, 0, 0], #q14
               [0, 0, 0, 0, 0, 0]]) #q15
         
 
"""Element 5"""
E5, A5, L5, I5 = 200e9, A, 1, I #Material Props
theta_deg5 = 0 #Angle roatated
theta_rad5 = np.radians(theta_deg5)
ke5 = local_frame(E5, A5, L5, I5, nu, As)
khat5 = global_frame(ke5, upplambda(theta_rad5))
             #D#1  2  3  4  5  6
A5 = np.array([[1, 0, 0, 0, 0, 0], #q1
               [0, 1, 0, 0, 0, 0], #q2
               [0, 0, 1, 0, 0, 0], #q3
               [0, 0, 0, 0, 0, 0], #q4
               [0, 0, 0, 0, 0, 0], #q5
               [0, 0, 0, 0, 0, 0], #q6
               [0, 0, 0, 0, 0, 0], #q7
               [0, 0, 0, 0, 0, 0], #q8
               [0, 0, 0, 0, 0, 0], #q9
               [0, 0, 0, 0, 0, 0], #q10
               [0, 0, 0, 0, 0, 0], #q11
               [0, 0, 0, 0, 0, 0], #q12
               [0, 0, 0, 1, 0, 0], #q13
               [0, 0, 0, 0, 1, 0], #q14
               [0, 0, 0, 0, 0, 1]]) #q15
         
 
"""Element 6"""
E6, A6, L6, I6 = 200e9, A, 2.6925, I #Material Props
theta_deg6 = 68.2 #Angle roatated
theta_rad6 = np.radians(theta_deg6)
ke6 = local_frame(E6, A6, L6, I6, nu, As)
khat6 = global_frame(ke6, upplambda(theta_rad6))
             #D#1  2  3  4  5  6
A6 = np.array([[1, 0, 0, 0, 0, 0], #q1
               [0, 1, 0, 0, 0, 0], #q2
               [0, 0, 1, 0, 0, 0], #q3
               [0, 0, 0, 0, 0, 0], #q4
               [0, 0, 0, 0, 0, 0], #q5
               [0, 0, 0, 0, 0, 0], #q6
               [0, 0, 0, 0, 0, 0], #q7
               [0, 0, 0, 0, 0, 0], #q8
               [0, 0, 0, 0, 0, 0], #q9
               [0, 0, 0, 1, 0, 0], #q10
               [0, 0, 0, 0, 1, 0], #q11
               [0, 0, 0, 0, 0, 1], #q12
               [0, 0, 0, 0, 0, 0], #q13
               [0, 0, 0, 0, 0, 0], #q14
               [0, 0, 0, 0, 0, 0]]) #q15
         
 
"""Element 7"""
E7, A7, L7, I7 = 200e9, A, 1, I #Material Props
theta_deg7 = 0 #Angle roatated
theta_rad7 = np.radians(theta_deg7)
ke7 = local_frame(E7, A7, L7, I7, nu, As)
khat7 = global_frame(ke7, upplambda(theta_rad7))
             #D#1  2  3  4  5  6
A7 = np.array([[0, 0, 0, 0, 0, 0], #q1
               [0, 0, 0, 0, 0, 0], #q2
               [0, 0, 0, 0, 0, 0], #q3
               [1, 0, 0, 0, 0, 0], #q4
               [0, 1, 0, 0, 0, 0], #q5
               [0, 0, 1, 0, 0, 0], #q6
               [0, 0, 0, 0, 0, 0], #q7
               [0, 0, 0, 0, 0, 0], #q8
               [0, 0, 0, 0, 0, 0], #q9
               [0, 0, 0, 1, 0, 0], #q10
               [0, 0, 0, 0, 1, 0], #q11
               [0, 0, 0, 0, 0, 1], #q12
               [0, 0, 0, 0, 0, 0], #q13
               [0, 0, 0, 0, 0, 0], #q14
               [0, 0, 0, 0, 0, 0]]) #q15
         
 
"""Element 8"""
E8, A8, L8, I8 = 200e9, A, 3, I #Material Props
theta_deg8 = -68.2 #Angle roatated
theta_rad8 = np.radians(theta_deg8)
ke8 = local_frame(E8, A8, L8, I8, nu, As)
khat8 = global_frame(ke8, upplambda(theta_rad8))
             #D#1  2  3  4  5  6
A8 = np.array([[0, 0, 0, 0, 0, 0], #q1
               [0, 0, 0, 0, 0, 0], #q2
               [0, 0, 0, 0, 0, 0], #q3
               [0, 0, 0, 0, 0, 0], #q4
               [0, 0, 0, 0, 0, 0], #q5
               [0, 0, 0, 0, 0, 0], #q6
               [1, 0, 0, 0, 0, 0], #q7
               [0, 1, 0, 0, 0, 0], #q8
               [0, 0, 1, 0, 0, 0], #q9
               [0, 0, 0, 1, 0, 0], #q10
               [0, 0, 0, 0, 1, 0], #q11
               [0, 0, 0, 0, 0, 1], #q12
               [0, 0, 0, 0, 0, 0], #q13
               [0, 0, 0, 0, 0, 0], #q14
               [0, 0, 0, 0, 0, 0]]) #q15
         
 
"""Element 9"""
E9, A9, L9, I9 = 200e9, A, 2.5, I #Material Props
theta_deg9 = 90 #Angle roatated
theta_rad9 = np.radians(theta_deg9)
ke9 = local_frame(E9, A9, L9, I9, nu, As)
khat9 = global_frame(ke9, upplambda(theta_rad9))
             #D#1  2  3  4  5  6
A9 = np.array([[0, 0, 0, 0, 0, 0], #q1
               [0, 0, 0, 0, 0, 0], #q2
               [0, 0, 0, 0, 0, 0], #q3
               [0, 0, 0, 0, 0, 0], #q4
               [0, 0, 0, 0, 0, 0], #q5
               [0, 0, 0, 0, 0, 0], #q6
               [0, 0, 0, 0, 0, 0], #q7
               [0, 0, 0, 0, 0, 0], #q8
               [0, 0, 0, 0, 0, 0], #q9
               [0, 0, 0, 0, 0, 0], #q10
               [0, 0, 0, 0, 0, 0], #q11
               [0, 0, 0, 0, 0, 0], #q12
               [0, 0, 0, 1, 0, 0], #q13
               [0, 0, 0, 0, 1, 0], #q14
               [0, 0, 0, 0, 0, 1]]) #q15
         
 
"""Element 10"""
E10, A10, L10, I10 = 200e9, A, 2.5, I #Material Props
theta_deg10 = 90 #Angle roatated
theta_rad10 = np.radians(theta_deg10)
ke10 = local_frame(E10, A10, L10, I10, nu, As)
khat10 = global_frame(ke10, upplambda(theta_rad10))
             #D#1  2  3  4  5  6
A10 = np.array([[0, 0, 0, 0, 0, 0], #q1
                [0, 0, 0, 0, 0, 0], #q2
                [0, 0, 0, 0, 0, 0], #q3
                [0, 0, 0, 0, 0, 0], #q4
                [0, 0, 0, 0, 0, 0], #q5
                [0, 0, 0, 0, 0, 0], #q6
                [0, 0, 0, 0, 0, 0], #q7
                [0, 0, 0, 0, 0, 0], #q8
                [0, 0, 0, 0, 0, 0], #q9
                [0, 0, 0, 1, 0, 0], #q10
                [0, 0, 0, 0, 1, 0], #q11
                [0, 0, 0, 0, 0, 1], #q12
                [1, 0, 0, 0, 0, 0], #q13
                [0, 1, 0, 0, 0, 0], #q14
                [0, 0, 1, 0, 0, 0]]) #q15
 
"""Indivdual KG Matrixs"""
KG1 = A1 @ khat1 @ A1.T
KG2 = A2 @ khat2 @ A2.T
KG3 = A3 @ khat3 @ A3.T
KG4 = A4 @ khat4 @ A4.T
KG5 = A5 @ khat5 @ A5.T
KG6 = A6 @ khat6 @ A6.T
KG7 = A7 @ khat7 @ A7.T
KG8 = A8 @ khat8 @ A8.T
KG9 = A9 @ khat9 @ A9.T
KG10 = A10 @ khat10 @ A10.T
 
"""Total"""
KG = KG1 + KG2 + KG3 + KG4 + KG5 + KG6 + KG7 + KG8 + KG9 + KG10
 
 
"""Q Loads"""
low = -4880
med = -6560
high = -9280
vhigh = -12000
maxwind = -2349.001858441632 * 8
load = maxwind
 
 
"""Uniform Distrubuted Load"""
Ludl1 = L3
wdl1 = load
thetaudl1 = theta_rad3
Audl1 = A3
 
f_udl1 = np.array([[0], 
                  [(wdl1*Ludl1)/2], 
                  [(wdl1*Ludl1**2)/12], 
                  [0], 
                  [(wdl1*Ludl1)/2], 
                  [(-wdl1*Ludl1**2)/12]])
 
F_udl1 = upplambda(thetaudl1).T @ f_udl1
Q_udl1 = Audl1 @ F_udl1
 
"""Uniform Distrubuted Load"""
Ludl = L2
wdl = load
thetaudl = theta_rad2
Audl = A2
 
f_udl = np.array([[0], 
                  [(wdl*Ludl)/2], 
                  [(wdl*Ludl**2)/12], 
                  [0], 
                  [(wdl*Ludl)/2], 
                  [(-wdl*Ludl**2)/12]])
 
F_udl = upplambda(thetaudl).T @ f_udl
Q_udl = Audl @ F_udl
 
"""Q loads from structure"""
Qtot = Q_udl + Q_udl1
q = np.linalg.solve(KG,  Qtot)
 
"""Forces"""
F1 = khat1 @ A1.T @ q
F2 = khat2 @ A2.T @ q
F3 = khat3 @ A3.T @ q
F4 = khat4 @ A4.T @ q
F5 = khat5 @ A5.T @ q
F6 = khat6 @ A6.T @ q
F7 = khat7 @ A7.T @ q
F8 = khat8 @ A8.T @ q
F9 = khat9 @ A9.T @ q
F10 = khat10 @ A10.T @ q
 
"""Deflections"""
d1 = upplambda(theta_rad1) @ A1.T @ q
d2 = upplambda(theta_rad2) @ A2.T @ q
d3 = upplambda(theta_rad3) @ A3.T @ q
d4 = upplambda(theta_rad4) @ A4.T @ q
d5 = upplambda(theta_rad5) @ A5.T @ q
d6 = upplambda(theta_rad6) @ A6.T @ q
d7 = upplambda(theta_rad7) @ A7.T @ q
d8 = upplambda(theta_rad8) @ A8.T @ q
d9 = upplambda(theta_rad9) @ A9.T @ q
d10 = upplambda(theta_rad10) @ A10.T @ q
 
"""Internal forces"""
f1 = ke1 @ d1
f2 = ke2 @ d2
f3 = ke3 @ d3
f4 = ke4 @ d4
f5 = ke5 @ d5
f6 = ke6 @ d6
f7 = ke7 @ d7
f8 = ke8 @ d8
f9 = ke9 @ d9
f10 = ke10 @ d10
 
"""Tricky plot"""
def plot_deflected_shape(node1XG, node1YG, node2XG, node2YG, d_e, Disp_mag, N_points, L, alpha):
    """plot full deflected shape"""
    x_e = np.linspace(0, L, N_points)
    
    """Axial shape functions"""
    psi1 = 1 - x_e / L
    psi2 = x_e / L
    
    """Tansverse shape functions"""
    N1 = 1 - (3*x_e**2)/L**2 + (2*x_e**3)/L**3
    N2 = x_e**3/L**2 - (2*x_e**2)/L + x_e
    N3 = (3*x_e**2)/L**2 - (2*x_e**3)/L**3
    N4 =x_e**3/L**2 - x_e**2/L
    
    """Displacements"""
    u_x = psi1*d_e[0,0] + psi2*d_e[3,0]
    v_x = N1*d_e[1,0] + N2*d_e[2,0] + N3*d_e[4,0] + N4*d_e[5,0]
    
    """Delfections"""
    Deflections_XG = u_x*np.cos(alpha) - v_x*np.sin(alpha)
    Deflections_YG = u_x*np.sin(alpha) + v_x*np.cos(alpha)
    
    """Undeflected baseline of elements"""
    Undeflected_baseline_XG = np.linspace(node1XG, node2XG, N_points)
    Undeflected_baseline_YG = np.linspace(node1YG, node2YG, N_points)  
    
    """Deflected position of all points within element"""
    Deflected_XG = Undeflected_baseline_XG + Disp_mag*Deflections_XG
    Deflected_YG = Undeflected_baseline_YG + Disp_mag*Deflections_YG
    
    """Plot"""
    plt.plot(Undeflected_baseline_XG, Undeflected_baseline_YG, 'b.-')
    plt.plot(Deflected_XG,Deflected_YG, 'r.-')
    
Disp_mag = 10
N_points = 10
 
element_1 = plot_deflected_shape(0, 0, 0, 2.5, d1, Disp_mag, N_points, L1, theta_rad1)
element_2 = plot_deflected_shape(0, 2.5, 0, 5, d2, Disp_mag, N_points, L2, theta_rad2)
element_3 = plot_deflected_shape(0, 5, 0, 7.5, d3, Disp_mag, N_points, L3, theta_rad3)
element_4 = plot_deflected_shape(0, 2.5, 1, 0, d4, Disp_mag, N_points, L4, theta_rad4)
element_5 = plot_deflected_shape(0, 2.5, 1, 2.5, d5, Disp_mag, N_points, L5, theta_rad5)
element_6 = plot_deflected_shape(0, 2.5, 1, 5, d6, Disp_mag, N_points, L6, theta_rad6)
element_7 = plot_deflected_shape(0, 5, 1, 5, d7, Disp_mag, N_points, L7, theta_rad7)
element_8 = plot_deflected_shape(0, 7.5, 1, 5, d8, Disp_mag, N_points, L8, theta_rad8)
element_9 = plot_deflected_shape(1, 0, 1, 2.5, d9, Disp_mag, N_points, L9, theta_rad9)
element_10 = plot_deflected_shape(1, 2.5, 1, 5, d10, Disp_mag, N_points, L10, theta_rad10)
 
plt.grid(True)
plt.xlim(-2, 3)      # extend x-axis to 3
plt.ylim(0, 8)      # set y-axis limit to fit your structure
plt.show()
 
 
 
"""Deflection at the tip"""
print("Deflection in y and x direction at the tip")
print(f"x = {float(q[6,0]) * 1e3:.3f} mm, y = {float(q[7,0]) * 1e3:.3f} mm, rotation = {float(q[8,0]) * 1e3:.3f}mrad\n")
 
 
# --- Reactions at supports from element end-actions (no loads applied at base DOFs) ---
 
# Left support at (0,0): element 1, i-end (local DOFs 0..2)
R_left_x = float(F1[0,0])
R_left_y = float(F1[1,0])
R_left_M = float(F1[2,0])
 
# Right support at (1,0): element 4, j-end (local DOFs 3..5) + element 9, i-end (0..2)
R_right_x = (float(F4[3,0]) + float(F9[0,0]))
R_right_y = (float(F4[4,0]) + float(F9[1,0]))
R_right_M = (float(F4[5,0]) + float(F9[2,0]))
 
print("Left support")
print(f"x = {R_left_x:.0f} N, y = {R_left_y:.0f} N, M = {R_left_M:.0f} Nm\n")
print("Right support")
print(f"x = {R_right_x:.0f} N, y = {R_right_y:.0f} N, M = {R_right_M:.0f} Nm\n")
 
Fx_applied = float((Q_udl + Q_udl1).sum())  # sum of all global DOF x-entries; or sum only x-DOFs if you prefer
Fx_react   = R_left_x + R_right_x
print(f"ΣFx applied = {Fx_applied:.1f} N,  ΣRx = {Fx_react:.1f} N,  Σ = {Fx_applied + Fx_react:.1f} N\n")
 
 
print(f"The axial force of element 1: {float(abs(f1[0,0])):.3} N")
print(f"The axial force of element 2: {float(abs(f2[0,0])):.3} N")
print(f"The axial force of element 3: {float(abs(f3[0,0])):.3} N")
print(f"The axial force of element 4: {float(abs(f4[0,0])):.3} N")
print(f"The axial force of element 5: {float(abs(f5[0,0])):.3} N")
print(f"The axial force of element 6: {float(abs(f6[0,0])):.3} N")
print(f"The axial force of element 7: {float(abs(f7[0,0])):.3} N")
print(f"The axial force of element 8: {float(abs(f8[0,0])):.3} N")
print(f"The axial force of element 9: {float(abs(f9[0,0])):.3} N")
print(f"The axial force of element 10: {float(abs(f10[0,0])):.3} N\n")
 
 
"""Bending + Axial stress for element 1"""
c = 70e-3      # [m] distance to outer fibre
A = 4.256e-3   # [m^2] cross-sectional area
 
Me1 = max(abs(float(f1[2,0])), abs(float(f1[5,0])))
Be1 = (Me1 * c) / I1
Af1 = abs(float(f1[0,0])) / A
sigma_tot1 = (Be1 + Af1) / 1e6
 
print("The axial + bending stress for each element")
print(f"Element 1 combined stress: {sigma_tot1:.3f} MPa")
 
 
"""Bending + Axial stress for element 2"""
Me2 = max(abs(float(f2[2,0])), abs(float(f2[5,0])))
Be2 = (Me2 * c) / I2
Af2 = abs(float(f2[0,0])) / A
sigma_tot2 = (Be2 + Af2) / 1e6
print(f"Element 2 combined stress: {sigma_tot2:.3f} MPa")
 
 
"""Bending + Axial stress for element 3"""
Me3 = max(abs(float(f3[2,0])), abs(float(f3[5,0])))
Be3 = (Me3 * c) / I3
Af3 = abs(float(f3[0,0])) / A
sigma_tot3 = (Be3 + Af3) / 1e6
print(f"Element 3 combined stress: {sigma_tot3:.3f} MPa")
 
 
"""Bending + Axial stress for element 4"""
Me4 = max(abs(float(f4[2,0])), abs(float(f4[5,0])))
Be4 = (Me4 * c) / I4
Af4 = abs(float(f4[0,0])) / A
sigma_tot4 = (Be4 + Af4) / 1e6
print(f"Element 4 combined stress: {sigma_tot4:.3f} MPa")
 
 
"""Bending + Axial stress for element 5"""
Me5 = max(abs(float(f5[2,0])), abs(float(f5[5,0])))
Be5 = (Me5 * c) / I5
Af5 = abs(float(f5[0,0])) / A
sigma_tot5 = (Be5 + Af5) / 1e6
print(f"Element 5 combined stress: {sigma_tot5:.3f} MPa")
 
 
"""Bending + Axial stress for element 6"""
Me6 = max(abs(float(f6[2,0])), abs(float(f6[5,0])))
Be6 = (Me6 * c) / I6
Af6 = abs(float(f6[0,0])) / A
sigma_tot6 = (Be6 + Af6) / 1e6
print(f"Element 6 combined stress: {sigma_tot6:.3f} MPa")
 
 
"""Bending + Axial stress for element 7"""
Me7 = max(abs(float(f7[2,0])), abs(float(f7[5,0])))
Be7 = (Me7 * c) / I7
Af7 = abs(float(f7[0,0])) / A
sigma_tot7 = (Be7 + Af7) / 1e6
print(f"Element 7 combined stress: {sigma_tot7:.3f} MPa")
 
 
"""Bending + Axial stress for element 8"""
Me8 = max(abs(float(f8[2,0])), abs(float(f8[5,0])))
Be8 = (Me8 * c) / I8
Af8 = abs(float(f8[0,0])) / A
sigma_tot8 = (Be8 + Af8) / 1e6
print(f"Element 8 combined stress: {sigma_tot8:.3f} MPa")
 
 
"""Bending + Axial stress for element 9"""
Me9 = max(abs(float(f9[2,0])), abs(float(f9[5,0])))
Be9 = (Me9 * c) / I9
Af9 = abs(float(f9[0,0])) / A
sigma_tot9 = (Be9 + Af9) / 1e6
print(f"Element 9 combined stress: {sigma_tot9:.3f} MPa")
 
 
"""Bending + Axial stress for element 10"""
Me10 = max(abs(float(f10[2,0])), abs(float(f10[5,0])))
Be10 = (Me10 * c) / I10
Af10 = abs(float(f10[0,0])) / A
sigma_tot10 = (Be10 + Af10) / 1e6
print(f"Element 10 combined stress: {sigma_tot10:.3f} MPa\n")
 
 
"""Fos stuff"""
### The load has to be set to -8
FOS = 2.5
yeild = 350
max_allowed = yeild / FOS
p_max = max_allowed / sigma_tot1
wind_speed = np.sqrt(p_max / 0.6)
print("The maximum wind speed with a FOS = 2.5 is:")
print(f"{wind_speed:.2f}m/s")
 
 
 
    
