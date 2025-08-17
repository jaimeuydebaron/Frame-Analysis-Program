# Tutorial 4 focuses on the calculation of the:
# Dealing with frames
# Element stiffness matrix in local coordinates for a frame element
# But this time it can handle UDL and ENL (Equivalent Nodal Loading)

import numpy as np
import matplotlib.pyplot as plt

# Element 1
length_1 = 10 # m
second_moment_of_area_1 = 5e-4 # m^4
area_1 = 1e-5 # m^2
youngs_modulus_1 = 200e9 # Pa
assembly_matrix_1 = np.array([[0, 0, 0, 1, 0, 0],
                              [0, 0, 0, 0, 0, 1],
                              [0, 0, 0, 0, 0, 0]]) # Generated manually

# Element 2 (same properties as Element 1)
length_2 = 10 # m
second_moment_of_area_2 = 5e-4 # m^4
area_2 = 1e-5 # m^2
youngs_modulus_2 = 200e9 # Pa
assembly_matrix_2 = np.array([[1, 0, 0, 0, 0, 0],
                              [0, 0, 1, 0, 0, 0],   
                              [0, 0, 0, 0, 0, 1]]) # Generated manually

global_forcing_vector = np.array([0, 140000, 0]) # Global forcing vector

def local_frame(youngs_modulus, area, second_moment_of_area, length):
    # Calculate the local stiffness matrix for a frame element
    constant = (youngs_modulus * second_moment_of_area) / length ** 3
    beta = area * length ** 2 / second_moment_of_area
    array = np.array([[beta, 0, 0, -beta, 0, 0],
                      [0, 12, 6 * length, 0, -12, 6 * length],
                      [0, 6 * length, 4 * length ** 2, 0, -6 * length, 2 * length ** 2],
                      [-beta, 0, 0, beta, 0, 0],
                      [0, -12, -6 * length, 0, 12, -6 * length],
                      [0, 6 * length, 2 * length ** 2, 0, -6 * length, 4 * length ** 2]])
    local_frame_stiffness = constant * array
    return local_frame_stiffness

local_frame_stiffness_1 = local_frame(youngs_modulus_1, area_1, second_moment_of_area_1, length_1)
local_frame_stiffness_2 = local_frame(youngs_modulus_2, area_2, second_moment_of_area_2, length_2)

def global_frame(local_stiffness, angle):
    # Calculate the global stiffness matrix for a frame element
    angle_rad = np.radians(angle)
    c = np.cos(angle_rad)
    s = np.sin(angle_rad)
    transformation_matrix = np.array([[c, s, 0, 0, 0, 0],
                                       [-s, c, 0, 0, 0, 0],
                                       [0, 0, 1, 0, 0, 0],
                                        [0, 0, 0, c, s, 0],
                                        [0, 0, 0, -s, c, 0],
                                        [0, 0, 0, 0, 0, 1]])
    global_stiffness = np.transpose(transformation_matrix) @ local_stiffness @ transformation_matrix

    # If a value is significantly smaller than 1e-10, set it to zero (CHECK THIS)
    global_stiffness[np.abs(global_stiffness) < 1e-9] = 0
    transformation_matrix[np.abs(transformation_matrix) < 1e-9] = 0

    # If a value is zero, set it to absolute zero (CHECK THIS)
    global_stiffness[np.isclose(global_stiffness, 0)] = 0
    transformation_matrix[np.isclose(transformation_matrix, 0)] = 0

    return global_stiffness, transformation_matrix

global_stiffness_1, transformation_matrix_1 = global_frame(local_frame_stiffness_1, 270)
global_stiffness_2, transformation_matrix_2 = global_frame(local_frame_stiffness_2, 0)

def local_stiffness_matrix(assembly_matrix, global_stiffness):
    # Calculate the stiffness matrix for a frame element using the assembly matrix
    local_stiffness_matrix = assembly_matrix @ global_stiffness @ np.transpose(assembly_matrix)
    return local_stiffness_matrix

local_stiffness_matrix_1 = local_stiffness_matrix(assembly_matrix_1, global_stiffness_1)
local_stiffness_matrix_2 = local_stiffness_matrix(assembly_matrix_2, global_stiffness_2)

def total_stiffness_matrix(stiffness_matrix_1, stiffness_matrix_2):
    # Calculate the total stiffness matrix for a frame element
    total_stiffness = stiffness_matrix_1 + stiffness_matrix_2
    return total_stiffness

total_stiffness_matrix = total_stiffness_matrix(local_stiffness_matrix_1, local_stiffness_matrix_2)

def find_displacement(stiffness_matrix, global_forcing_vector):
    # Calculate the displacement vector using the stiffness matrix and global forcing vector
    displacement = np.linalg.solve(stiffness_matrix, global_forcing_vector)
    return displacement

displacement_vector = find_displacement(total_stiffness_matrix, global_forcing_vector)

def structural_deflections(displacement_vector, assembly_matrix):
    # Calculate the element structural deflections in global coordinates
    structural_deflections = np.transpose(assembly_matrix) @ displacement_vector
    return structural_deflections

structural_deflections_1 = structural_deflections(displacement_vector, assembly_matrix_1)
structural_deflections_2 = structural_deflections(displacement_vector, assembly_matrix_2)

def element_structural_deflections_local(structural_deflection, transformation_matrix):
    # Calculate the element structural deflections in local coordinates
    local_structural_deflections = transformation_matrix @ structural_deflection
    return local_structural_deflections

element_structural_deflections_local_1 = element_structural_deflections_local(structural_deflections_1, transformation_matrix_1)
element_structural_deflections_local_2 = element_structural_deflections_local(structural_deflections_2, transformation_matrix_2) 

def element_nodal_forces_local(local_stiffness, element_structural_deflections_local):
    # Calculate the element nodal forces in global coordinates
    nodal_forces = local_stiffness @ element_structural_deflections_local
    return nodal_forces

element_nodal_forces_1 = element_nodal_forces_local(local_frame_stiffness_1, element_structural_deflections_local_1)
element_nodal_forces_2 = element_nodal_forces_local(local_frame_stiffness_2, element_structural_deflections_local_2)

def element_nodal_forces_global(global_stiffness, structural_deflections):
    # Calculate the element nodal forces in global coordinates
    nodal_forces = global_stiffness @ structural_deflections
    return nodal_forces

element_nodal_forces_global_1 = element_nodal_forces_global(global_stiffness_1, structural_deflections_1)
element_nodal_forces_global_2 = element_nodal_forces_global(global_stiffness_2, structural_deflections_2)

print(local_frame_stiffness_1)