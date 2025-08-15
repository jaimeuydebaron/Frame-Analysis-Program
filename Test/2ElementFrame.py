import numpy as np
import matplotlib.pyplot as plt

# 2 Element Frame Analysis

# Function calculations
# Calculate the stiffness matrix for a frame element in LOCAL coordinates
def calc_local_coordinates_element_stiffness_matrix(youngs_modulus, area, second_moment_of_area, length):
    constant = (youngs_modulus * second_moment_of_area) / length ** 3
    beta = area * length ** 2 / second_moment_of_area
    array = np.array([[beta, 0, 0, -beta, 0, 0],
                      [0, 12, 6 * length, 0, -12, 6 * length],
                      [0, 6 * length, 4 * length ** 2, 0, -6 * length, 2 * length ** 2],
                      [-beta, 0, 0, beta, 0, 0],
                      [0, -12, -6 * length, 0, 12, -6 * length],
                      [0, 6 * length, 2 * length ** 2, 0, -6 * length, 4 * length ** 2]])
    local_stiffness = constant * array
    return local_stiffness

# Calculate the stiffness matrix for a frame element in GLOBAL coordinates
def calc_global_coordinates_element_stiffness_matrix(local_stiffness, angle):
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

# Calculate the global stifness matrix in the entire structure
def calc_global_stiffness_matrix(assembly_matrix, global_stiffness):
    return assembly_matrix @ global_stiffness @ np.transpose(assembly_matrix)

# Calculate the total stiffness matrix for a structure with THREE elements
def calc_total_global_stiffness_matrix(global_stiffness_matrix_1, global_stiffness_matrix_2):
    return global_stiffness_matrix_1 + global_stiffness_matrix_2

# Calculate the OVERALL structural deflections
def calc_structural_deflections(total_stiffness_global_coordinates, global_forcing_vector):
    return np.linalg.solve(total_stiffness_global_coordinates, global_forcing_vector)

# Calculate the element deflections in GLOBAL coordinates
def calc_global_element_structural_deflections(total_structural_deflections, assembly_matrix):
    return np.transpose(assembly_matrix) @ total_structural_deflections

# Calculate the element deflections in LOCAL coordinates
def calc_local_element_structural_deflections(global_structural_deflections, transformation_matrix):
    return transformation_matrix @ global_structural_deflections

# Calculate the element nodal forces in LOCAL coordinates
def calc_local_element_nodal_forces(local_stiffness, local_structural_deflections):
    return local_stiffness @ local_structural_deflections

# Calculate the element nodal forces in GLOBAL coordinates
def calc_global_element_nodal_forces(global_stiffness, global_structural_deflections):
    return global_stiffness @ global_structural_deflections

""" IF THERE IS A UDL: """

def build_UDL_global_forcing_vector(w_bar, L, assembly_matrix, transformation_matrix):
    equivalent_forcing_element_coordinates = np.array([0,
                                                    w_bar * L / 2,
                                                    w_bar * L ** 2 / 12,
                                                    0,
                                                    w_bar * L / 2,
                                                    -w_bar * L ** 2 / 12])
    big_F = np.transpose(transformation_matrix) @ equivalent_forcing_element_coordinates
    UDL_global_forcing_vector = big_F @ assembly_matrix
    return UDL_global_forcing_vector

""" IF THERE IS A PDL: """

def build_PDL_global_forcing_vector(P, L, assembly_matrix, transformation_matrix):
    equivalent_forcing_element_coordinates = np.array([0,
                                                    P / 2,
                                                    P * L / 8,
                                                    0,
                                                    P / 2,
                                                    -P * L / 8])
    big_F = np.transpose(transformation_matrix) @ equivalent_forcing_element_coordinates
    PDL_global_forcing_vector = big_F @ assembly_matrix
    return PDL_global_forcing_vector


# ___________________________________________________________________________________________________

# Element information

# Element 1
length_1 = 10 # m
second_moment_of_area_1 = 5e-4 # m^4
area_1 = 1e-5 # m^2
youngs_modulus_1 = 200e9 # Pa
assembly_matrix_1 = np.array([[0, 0, 0, 1, 0, 0],
                              [0, 0, 0, 0, 0, 1],
                              [0, 0, 0, 0, 0, 0]]) # Generated manually
angle_1 = 270
element_1_node_1_location = np.array([0, 0])  # Node 1 location in global coordinates
element_1_node_2_location = np.array([0, length_1])  # Node 2 location in global coordinates

# Element 2 (same properties as Element 1)
length_2 = 10 # m
second_moment_of_area_2 = 5e-4 # m^4
area_2 = 1e-5 # m^2
youngs_modulus_2 = 200e9 # Pa
assembly_matrix_2 = np.array([[1, 0, 0, 0, 0, 0],
                              [0, 0, 1, 0, 0, 0],   
                              [0, 0, 0, 0, 0, 1]]) # Generated manually
angle_2 = 0
element_2_node_1_location = np.array([0, 0])  # Node 1 location in global coordinates
element_2_node_2_location = np.array([length_2, 0])  # Node 2 location in global coordinates

""" UDL AND PDL INFORMATION """

# Global forcing vector NORMAL
normal_global_forcing_vector = np.array([0, 140000, 0])  # Nodal forces

# # Uniformly Distributed Load (UDL)
# w_bar = -10000  # The uniform distributed load (UDL) in N/m
# L_UDL = 4.5  # Length of the bar it is acting on

# global_forcing_vector_UDL = build_UDL_global_forcing_vector(w_bar, L_UDL, assembly_matrix_2, transformation_matrix_2)

# # Point Distributed Load (PDL)
# P = 20000  # The point distributed load (PDL) in N
# global_forcing_vector_PDL = build_PDL_global_forcing_vector(P, length_2, assembly_matrix_2, transformation_matrix_2)

# global_forcing_vector = normal_global_forcing_vector + global_forcing_vector_UDL + global_forcing_vector_PDL

global_forcing_vector = normal_global_forcing_vector  # Use only the normal forces for this example

# ___________________________________________________________________________________________________

# Perform calculations
local_stiffness_1 = calc_local_coordinates_element_stiffness_matrix(youngs_modulus_1, area_1, second_moment_of_area_1, length_1)
local_stiffness_2 = calc_local_coordinates_element_stiffness_matrix(youngs_modulus_2, area_2, second_moment_of_area_2, length_2)

global_stiffness_1, transformation_matrix_1 = calc_global_coordinates_element_stiffness_matrix(local_stiffness_1, angle_1)
global_stiffness_2, transformation_matrix_2 = calc_global_coordinates_element_stiffness_matrix(local_stiffness_2, angle_2)

global_stiffness_matrix_1 = calc_global_stiffness_matrix(assembly_matrix_1, global_stiffness_1)
global_stiffness_matrix_2 = calc_global_stiffness_matrix(assembly_matrix_2, global_stiffness_2)

total_global_stiffness_matrix = calc_total_global_stiffness_matrix(global_stiffness_matrix_1, global_stiffness_matrix_2)

total_structural_deflections = calc_structural_deflections(total_global_stiffness_matrix, global_forcing_vector)

global_structural_deflections_1 = calc_global_element_structural_deflections(total_structural_deflections, assembly_matrix_1)
global_structural_deflections_2 = calc_global_element_structural_deflections(total_structural_deflections, assembly_matrix_2)

local_structural_deflections_1 = calc_local_element_structural_deflections(global_structural_deflections_1, transformation_matrix_1)
local_structural_deflections_2 = calc_local_element_structural_deflections(global_structural_deflections_2, transformation_matrix_2)

element_nodal_forces_local_1 = calc_local_element_nodal_forces(local_stiffness_1, local_structural_deflections_1)
element_nodal_forces_local_2 = calc_local_element_nodal_forces(local_stiffness_2, local_structural_deflections_2)

element_nodal_forces_global_1 = calc_global_element_nodal_forces(global_stiffness_1, global_structural_deflections_1)
element_nodal_forces_global_2 = calc_global_element_nodal_forces(global_stiffness_2, global_structural_deflections_2)

# ___________________________________________________________________________________________________

print(element_nodal_forces_global_1)

# ___________________________________________________________________________________________________

# Plot the structure
plt.figure(figsize=(10, 6))
plt.plot([element_1_node_1_location[0], element_2_node_1_location[0]], 
         [element_1_node_1_location[1], element_2_node_1_location[1]], 'ro', label='Node 1')
plt.plot([element_1_node_2_location[0], element_2_node_2_location[0]], 
         [element_1_node_2_location[1], element_2_node_2_location[1]], 'bo', label='Node 2')
plt.plot([element_1_node_1_location[0], element_1_node_2_location[0]], 
         [element_1_node_1_location[1], element_1_node_2_location[1]], 'k-', label='Element 1')
plt.plot([element_2_node_1_location[0], element_2_node_2_location[0]], 
         [element_2_node_1_location[1], element_2_node_2_location[1]], 'g-', label='Element 2')
plt.title('2 Element Frame Structure')
plt.xlabel('X-axis (m)')
plt.ylabel('Y-axis (m)')
plt.legend()
plt.grid()
plt.axis('equal')

# Plot the deflected shape


plt.show()
# ___________________________________________________________________________________________________