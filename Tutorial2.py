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

local_stiffness_1 = local_bar(youngs_modulus_1, cross_sectional_area_1, length1)
local_stiffness_2 = local_bar(youngs_modulus_2, cross_sectional_area_2, length2)

def global_bar(local_stiffness, angle):
    # Calculate the global stiffness matrix for a 2D bar element
    angle_rad = np.radians(angle)
    c = np.cos(angle_rad)
    s = np.sin(angle_rad)
    transformation_matrix = np.array([[c, s, 0, 0],
                                       [0, 0, c, s]])
    global_stiffness = np.transpose(transformation_matrix) @ local_stiffness @ transformation_matrix
    return global_stiffness, transformation_matrix

global_stiffness_1, transformation_matrix_1 = global_bar(local_stiffness_1, alpha_1)
global_stiffness_2, transformation_matrix_2 = global_bar(local_stiffness_2, alpha_2)

def local_stiffness_matrix(assembly_matrix, global_stiffness):
    # Calculate the stiffness matrix for a 2D bar element using the assembly matrix
    local_stiffness_matrix =  assembly_matrix @ global_stiffness @ np.transpose(assembly_matrix)
    return local_stiffness_matrix

local_stiffness_matrix_1 = local_stiffness_matrix(assembly_matrix_1, global_stiffness_1)
local_stiffness_matrix_2 = local_stiffness_matrix(assembly_matrix_2, global_stiffness_2)

def total_stiffness_matrix(stiffness_matrix_1, stiffness_matrix_2):
    # Calculate the total stiffness matrix for a 2D bar element
    total_stiffness = stiffness_matrix_1 + stiffness_matrix_2
    return total_stiffness  

total_stiffness_matrix = total_stiffness_matrix(local_stiffness_matrix_1, local_stiffness_matrix_2)

def find_displacement_vector(Q1, Q2, total_stiffness):
    # Calculate the displacement vector for a 2D bar element
    force_vector = np.array([[Q1], [Q2]])
    displacement_vector = np.linalg.solve(total_stiffness, force_vector)
    return displacement_vector

displacement_vector = find_displacement_vector(Q1_1, Q2_1, total_stiffness_matrix)

def element_nodal_forces(global_bar, displacement_vector, assembly_matrix):
    # Calculate the element nodal forces in global coordinates
    nodal_forces = global_bar @ np.transpose(assembly_matrix) @ displacement_vector
    return nodal_forces

element_nodal_forces_1 = element_nodal_forces(global_stiffness_1, displacement_vector, assembly_matrix_1)
element_nodal_forces_2 = element_nodal_forces(global_stiffness_2, displacement_vector, assembly_matrix_2)

def element_nodal_forces_local(displacement_vector, assembly_matrix, transformation_matrix):
    # Calculate the element nodal forces in local coordinates
    local_nodal_forces = transformation_matrix @ np.transpose(assembly_matrix) @ displacement_vector
    return local_nodal_forces

element_nodal_forces_local_1 = element_nodal_forces_local(displacement_vector, assembly_matrix_1, transformation_matrix_1)
element_nodal_forces_local_2 = element_nodal_forces_local(displacement_vector, assembly_matrix_2, transformation_matrix_2)

def element_strain(element_nodal_forces_local, length):
    # Calculate the element internal forces in local coordinates
    strain = (element_nodal_forces_local[1] - element_nodal_forces_local[0]) / length
    return strain

element_strain_1 = element_strain(element_nodal_forces_local_1, length1)
element_strain_2 = element_strain(element_nodal_forces_local_2, length2)

def element_nodal_force_local(element_nodal_forces_local, local_stiffness):
    # Calculate the element nodal forces in local coordinates
    local_nodal_force = local_stiffness @ element_nodal_forces_local
    return local_nodal_force

element_nodal_force_local_1 = element_nodal_force_local(element_nodal_forces_local_1, local_stiffness_1)
element_nodal_force_local_2 = element_nodal_force_local(element_nodal_forces_local_2, local_stiffness_2)

print(element_nodal_force_local_2)

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