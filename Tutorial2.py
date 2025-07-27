# Tutorial 2 focuses on the calculation of the:
# Bar element stiffness matrix for a 2D bar element
# Total stiffness matrix for a 2D bar element
# Displacement vector for a 2D bar element
# Element nodal forces in global coordinates
# Element nodal forces in local coordinates
# Element internal forces in local coordinates

import numpy as np
import matplotlib.pyplot as plt

# Element 1
diameter_1 = 0.1
youngs_modulus_1 = 2e11
length1 = 1
Q1_1 = 0
Q2_1 = 100000
alpha_1 = 0
assembly_matrix_1 = np.array([[0 , 0, 1, 0],
                           [0, 0, 0, 1]]) # Generated manually

# Element 2
diameter_2 = 0.1
youngs_modulus_2 = 2e11
length2 = 1.41
Q1_2 = 0
Q2_2 = 100000
alpha_2 = 45
assembly_matrix_2 = np.array([[0 , 0, 1, 0],
                            [0, 0, 0, 1]]) # Generated manually

cross_sectional_area_1 = np.pi * (diameter_1  / 2) ** 2
cross_sectional_area_2 = np.pi * (diameter_2  / 2) ** 2

def local_bar(youngs_modulus, cross_sectional_area, length):
    # Calculate the local stiffness matrix for a 2D bar element
    local_stiffness = (youngs_modulus * cross_sectional_area) / length * np.array([[1, -1], [-1, 1]])
    return local_stiffness

def global_bar(local_stiffness, angle):
    # Calculate the global stiffness matrix for a 2D bar element
    angle_rad = np.radians(angle)
    c = np.cos(angle_rad)
    s = np.sin(angle_rad)
    transformation_matrix = np.array([[c, s, 0, 0],
                                       [0, 0, c, s]])
    global_stiffness = np.transpose(transformation_matrix) @ local_stiffness @ transformation_matrix
    return global_stiffness, transformation_matrix

def local_stiffness_matrix(assembly_matrix, global_stiffness):
    # Calculate the stiffness matrix for a 2D bar element using the assembly matrix
    stiffness_matrix =  assembly_matrix @ global_stiffness @ np.transpose(assembly_matrix)
    return stiffness_matrix

def total_stiffness_matrix(stiffness_matrix_1, stiffness_matrix_2):
    # Calculate the total stiffness matrix for a 2D bar element
    total_stiffness = stiffness_matrix_1 + stiffness_matrix_2
    return total_stiffness  

def find_displacement_vector(Q1, Q2, total_stiffness):
    # Calculate the displacement vector for a 2D bar element
    force_vector = np.array([[Q1], [Q2]])
    displacement_vector = np.linalg.solve(total_stiffness, force_vector)
    return displacement_vector

def element_nodal_forces(displacement_vector, global_stiffness, transformation_matrix):
    # Calculate the element nodal forces in global coordinates
    local_displacement = transformation_matrix @ displacement_vector
    element_nodal_forces = global_stiffness @ local_displacement
    return element_nodal_forces

# To do:
# Finish up from element nodal forces function

# print(f"\n")
# print(f"K1 = 1 x 10^8 * {K1/1e8}\n")
# print(f"K1hat = 1 x 10^8 * {K1hat/1e8}\n")
# print(f"Lambda1 = {Lambda1}\n")

# print(f"K2 = 1 x 10^8 * {K2/1e8}\n")
# print(f"K2hat = 1 x 10^7 * {K2hat/1e7}\n")
# print(f"Lambda2 = {Lambda2}\n")

# print(f"KG1 = 1 x 10^8 * {KG1/1e8}\n")
# print(f"KG2 = 1 x 10^7 * {KG2/1e7}\n")
# print(f"KG = 1 x 10^8 * {KG/1e8}\n")