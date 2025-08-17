# Tutorial 3 focuses on the calculation of the:
# Dealing with frames
# Element stiffness matrix in local coordinates for a frame element

import numpy as np
import matplotlib.pyplot as plt

# Element 1
length_1 = 3 # m
second_moment_of_area_1 = 1e-5 # m^4
area_1 = 5e-4 # m^2
youngs_modulus_1 = 200e9 # Pa
assembly_matrix_1 = np.array([[0, 0, 0, 1, 0, 0],
                              [0, 0, 0, 0, 1, 0],
                              [0, 0, 0, 0, 0, 1],
                              [0, 0, 0, 0, 0, 0],
                              [0, 0, 0, 0, 0, 0],
                              [0, 0, 0, 0, 0, 0]]) # Generated manually
angle_1 = 90

# Element 2 (same properties as Element 1)
length_2 = 4.5 # m
second_moment_of_area_2 = 1e-5 # m^4
area_2 = 5e-4 # m^2
youngs_modulus_2 = 200e9 # Pa
assembly_matrix_2 = np.array([[1, 0, 0, 0, 0, 0],
                              [0, 1, 0, 0, 0, 0],   
                              [0, 0, 1, 0, 0, 0],
                              [0, 0, 0, 1, 0, 0],
                              [0, 0, 0, 0, 1, 0],
                              [0, 0, 0, 0, 0, 1]]) # Generated manually
angle_2 = 0

# Element 3
length_3 = 3 # m
second_moment_of_area_3 = 1e-5 # m^4
area_3 = 5e-4 # m^2
youngs_modulus_3 = 200e9 # Pa
assembly_matrix_3 = np.array([[0, 0, 0, 0, 0, 0],
                              [0, 0, 0, 0, 0, 0],
                              [0, 0, 0, 0, 0, 0],
                              [1, 0, 0, 0, 0, 0],
                              [0, 1, 0, 0, 0, 0],   
                              [0, 0, 1, 0, 0, 0]]) # Generated manually
angle_3 = 270

# Building the global forcing vector


Q_nodal = np.array([10000, 0, 0, 10000, 0, 0]) # Nodal forces

# UDL
w_bar = -10000 # The uniform distributed load (UDL) in N/m
L_UDL = 4.5 # Length of the bar it is acting on

def build_UDL(w_bar, L, assembly_matrix, transformation_matrix):
    equivalent_forcing_element_coordinates = np.array([0,
                                                    w_bar * L / 2,
                                                    w_bar * L ** 2 / 12,
                                                    0,
                                                    w_bar * L / 2,
                                                    -w_bar * L ** 2 / 12])
    big_F = np.transpose(transformation_matrix) @ equivalent_forcing_element_coordinates
    Q = big_F @ assembly_matrix
    return Q

# PDL
P = -50000 # The point distributed load (PDL) in N
L = 4.5 # Length of the bar it is acting on

def build_PDL(P, L, assembly_matrix, transformation_matrix):
    equivalent_forcing_element_coordinates = np.array([0,
                                                    P / 2,
                                                    P * L / 8,
                                                    0,
                                                    P / 2,
                                                    -P * L / 8])
    big_F = np.transpose(transformation_matrix) @ equivalent_forcing_element_coordinates
    Q = big_F @ assembly_matrix
    return Q

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
local_frame_stiffness_3 = local_frame(youngs_modulus_3, area_3, second_moment_of_area_3, length_3)

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

global_stiffness_1, transformation_matrix_1 = global_frame(local_frame_stiffness_1, angle_1)
global_stiffness_2, transformation_matrix_2 = global_frame(local_frame_stiffness_2, angle_2)
global_stiffness_3, transformation_matrix_3 = global_frame(local_frame_stiffness_3, angle_3)

def local_stiffness_matrix(assembly_matrix, global_stiffness):
    # Calculate the stiffness matrix for a frame element using the assembly matrix
    local_stiffness_matrix = assembly_matrix @ global_stiffness @ np.transpose(assembly_matrix)
    return local_stiffness_matrix

local_stiffness_matrix_1 = local_stiffness_matrix(assembly_matrix_1, global_stiffness_1)
local_stiffness_matrix_2 = local_stiffness_matrix(assembly_matrix_2, global_stiffness_2)
local_stiffness_matrix_3 = local_stiffness_matrix(assembly_matrix_3, global_stiffness_3)

def total_stiffness_matrix(stiffness_matrix_1, stiffness_matrix_2, stiffness_matrix_3):
    # Calculate the total stiffness matrix for a frame element
    total_stiffness = stiffness_matrix_1 + stiffness_matrix_2 + stiffness_matrix_3
    return total_stiffness

total_stiffness_matrix = total_stiffness_matrix(local_stiffness_matrix_1, local_stiffness_matrix_2, local_stiffness_matrix_3)

def find_displacement(stiffness_matrix, global_forcing_vector):
    # Calculate the displacement vector using the stiffness matrix and global forcing vector
    displacement = np.linalg.solve(stiffness_matrix, global_forcing_vector)
    return displacement

global_forcing_vector_UDL = build_UDL(w_bar, L_UDL, assembly_matrix_2, transformation_matrix_2)
global_forcing_vector_PDL = build_PDL(P, L, assembly_matrix_2, transformation_matrix_2)
global_forcing_vector = global_forcing_vector_UDL + global_forcing_vector_PDL + Q_nodal

displacement_vector = find_displacement(total_stiffness_matrix, global_forcing_vector)

def structural_deflections(displacement_vector, assembly_matrix):
    # Calculate the element structural deflections in global coordinates
    structural_deflections = np.transpose(assembly_matrix) @ displacement_vector
    return structural_deflections

structural_deflections_1 = structural_deflections(displacement_vector, assembly_matrix_1)
structural_deflections_2 = structural_deflections(displacement_vector, assembly_matrix_2)
structural_deflections_3 = structural_deflections(displacement_vector, assembly_matrix_3)

def element_structural_deflections_local(structural_deflection, transformation_matrix):
    # Calculate the element structural deflections in local coordinates
    local_structural_deflections = transformation_matrix @ structural_deflection
    return local_structural_deflections

element_structural_deflections_local_1 = element_structural_deflections_local(structural_deflections_1, transformation_matrix_1)
element_structural_deflections_local_2 = element_structural_deflections_local(structural_deflections_2, transformation_matrix_2)
element_structural_deflections_local_3 = element_structural_deflections_local(structural_deflections_3, transformation_matrix_3)

def element_nodal_forces_local(local_stiffness, element_structural_deflections_local):
    # Calculate the element nodal forces in global coordinates
    nodal_forces = local_stiffness @ element_structural_deflections_local
    return nodal_forces

element_nodal_forces_1 = element_nodal_forces_local(local_frame_stiffness_1, element_structural_deflections_local_1)
element_nodal_forces_2 = element_nodal_forces_local(local_frame_stiffness_2, element_structural_deflections_local_2)
element_nodal_forces_3 = element_nodal_forces_local(local_frame_stiffness_3, element_structural_deflections_local_3)

def element_nodal_forces_global(global_stiffness, structural_deflections):
    # Calculate the element nodal forces in global coordinates
    nodal_forces = global_stiffness @ structural_deflections
    return nodal_forces

element_nodal_forces_global_1 = element_nodal_forces_global(global_stiffness_1, structural_deflections_1)
element_nodal_forces_global_2 = element_nodal_forces_global(global_stiffness_2, structural_deflections_2)
element_nodal_forces_global_3 = element_nodal_forces_global(global_stiffness_3, structural_deflections_3)

print(element_nodal_forces_global_1)