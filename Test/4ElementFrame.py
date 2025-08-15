import numpy as np
import matplotlib.pyplot as plt

# 4 Element Frame Analysis

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
    return constant * array

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
    global_stiffness[np.abs(global_stiffness) < 1e-9] = 0
    transformation_matrix[np.abs(transformation_matrix) < 1e-9] = 0
    return global_stiffness, transformation_matrix

# Calculate the global stiffness matrix for an element
def calc_global_stiffness_matrix(assembly_matrix, global_stiffness):
    return assembly_matrix @ global_stiffness @ np.transpose(assembly_matrix)

# Calculate the total stiffness matrix for the structure
def calc_total_global_stiffness_matrix(global_stiffness_matrix_1, global_stiffness_matrix_2, global_stiffness_matrix_3, global_stiffness_matrix_4):
    return global_stiffness_matrix_1 + global_stiffness_matrix_2 + global_stiffness_matrix_3 + global_stiffness_matrix_4

# Calculate the overall structural deflections
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

# ___________________________________________________________________________________________________

# Element 1
length_1 = 10  # m
second_moment_of_area_1 = 5e-4  # m^4
area_1 = 1e-5  # m^2
youngs_modulus_1 = 200e9  # Pa
assembly_matrix_1 = np.array([[0, 0, 0, 1, 0, 0],
                              [0, 0, 0, 0, 0, 1],
                              [0, 0, 0, 0, 0, 0]])
angle_1 = 270
element_1_node_1_location = np.array([0, 0])
element_1_node_2_location = np.array([0, length_1])

# Element 2
length_2 = 10  # m
second_moment_of_area_2 = 5e-4
area_2 = 1e-5
youngs_modulus_2 = 200e9
assembly_matrix_2 = np.array([[1, 0, 0, 0, 0, 0],
                              [0, 0, 1, 0, 0, 0],
                              [0, 0, 0, 0, 0, 1]])
angle_2 = 0
element_2_node_1_location = np.array([0, 0])
element_2_node_2_location = np.array([length_2, 0])

# Element 3
length_3 = 10
second_moment_of_area_3 = 5e-4
area_3 = 1e-5
youngs_modulus_3 = 200e9
assembly_matrix_3 = np.array([[0, 0, 0, 0, 1, 0],
                              [0, 0, 0, 0, 0, 0],
                              [0, 0, 0, 1, 0, 0]])
angle_3 = 90
element_3_node_1_location = np.array([length_2, 0])
element_3_node_2_location = np.array([length_2, length_3])

# Element 4
length_4 = 10
second_moment_of_area_4 = 5e-4
area_4 = 1e-5
youngs_modulus_4 = 200e9
assembly_matrix_4 = np.array([[0, 0, 0, 0, 1, 0],
                              [0, 0, 0, 0, 0, 1],
                              [0, 0, 0, 0, 0, 0]])
angle_4 = 0
element_4_node_1_location = np.array([0, length_1])
element_4_node_2_location = np.array([length_4, length_1])

# ___________________________________________________________________________________________________

# Global forcing vector
global_forcing_vector = np.array([0, 140000, 0])

# ___________________________________________________________________________________________________

# Perform calculations
local_stiffness_1 = calc_local_coordinates_element_stiffness_matrix(youngs_modulus_1, area_1, second_moment_of_area_1, length_1)
local_stiffness_2 = calc_local_coordinates_element_stiffness_matrix(youngs_modulus_2, area_2, second_moment_of_area_2, length_2)
local_stiffness_3 = calc_local_coordinates_element_stiffness_matrix(youngs_modulus_3, area_3, second_moment_of_area_3, length_3)
local_stiffness_4 = calc_local_coordinates_element_stiffness_matrix(youngs_modulus_4, area_4, second_moment_of_area_4, length_4)

global_stiffness_1, transformation_matrix_1 = calc_global_coordinates_element_stiffness_matrix(local_stiffness_1, angle_1)
global_stiffness_2, transformation_matrix_2 = calc_global_coordinates_element_stiffness_matrix(local_stiffness_2, angle_2)
global_stiffness_3, transformation_matrix_3 = calc_global_coordinates_element_stiffness_matrix(local_stiffness_3, angle_3)
global_stiffness_4, transformation_matrix_4 = calc_global_coordinates_element_stiffness_matrix(local_stiffness_4, angle_4)

global_stiffness_matrix_1 = calc_global_stiffness_matrix(assembly_matrix_1, global_stiffness_1)
global_stiffness_matrix_2 = calc_global_stiffness_matrix(assembly_matrix_2, global_stiffness_2)
global_stiffness_matrix_3 = calc_global_stiffness_matrix(assembly_matrix_3, global_stiffness_3)
global_stiffness_matrix_4 = calc_global_stiffness_matrix(assembly_matrix_4, global_stiffness_4)

total_global_stiffness_matrix = calc_total_global_stiffness_matrix(global_stiffness_matrix_1, global_stiffness_matrix_2, global_stiffness_matrix_3, global_stiffness_matrix_4)

total_structural_deflections = calc_structural_deflections(total_global_stiffness_matrix, global_forcing_vector)

global_structural_deflections_1 = calc_global_element_structural_deflections(total_structural_deflections, assembly_matrix_1)
global_structural_deflections_2 = calc_global_element_structural_deflections(total_structural_deflections, assembly_matrix_2)
global_structural_deflections_3 = calc_global_element_structural_deflections(total_structural_deflections, assembly_matrix_3)
global_structural_deflections_4 = calc_global_element_structural_deflections(total_structural_deflections, assembly_matrix_4)

local_structural_deflections_1 = calc_local_element_structural_deflections(global_structural_deflections_1, transformation_matrix_1)
local_structural_deflections_2 = calc_local_element_structural_deflections(global_structural_deflections_2, transformation_matrix_2)
local_structural_deflections_3 = calc_local_element_structural_deflections(global_structural_deflections_3, transformation_matrix_3)
local_structural_deflections_4 = calc_local_element_structural_deflections(global_structural_deflections_4, transformation_matrix_4)

element_nodal_forces_local_1 = calc_local_element_nodal_forces(local_stiffness_1, local_structural_deflections_1)
element_nodal_forces_local_2 = calc_local_element_nodal_forces(local_stiffness_2, local_structural_deflections_2)
element_nodal_forces_local_3 = calc_local_element_nodal_forces(local_stiffness_3, local_structural_deflections_3)
element_nodal_forces_local_4 = calc_local_element_nodal_forces(local_stiffness_4, local_structural_deflections_4)

element_nodal_forces_global_1 = calc_global_element_nodal_forces(global_stiffness_1, global_structural_deflections_1)
element_nodal_forces_global_2 = calc_global_element_nodal_forces(global_stiffness_2, global_structural_deflections_2)
element_nodal_forces_global_3 = calc_global_element_nodal_forces(global_stiffness_3, global_structural_deflections_3)
element_nodal_forces_global_4 = calc_global_element_nodal_forces(global_stiffness_4, global_structural_deflections_4)

# ___________________________________________________________________________________________________

print(element_nodal_forces_global_1)
print(element_nodal_forces_global_2)
print(element_nodal_forces_global_3)
print(element_nodal_forces_global_4)

# ___________________________________________________________________________________________________

# Plot the structure
plt.figure(figsize=(10, 6))
plt.plot([element_1_node_1_location[0], element_1_node_2_location[0]], [element_1_node_1_location[1], element_1_node_2_location[1]], 'k-', label='Element 1')
plt.plot([element_2_node_1_location[0], element_2_node_2_location[0]], [element_2_node_1_location[1], element_2_node_2_location[1]], 'g-', label='Element 2')
plt.plot([element_3_node_1_location[0], element_3_node_2_location[0]], [element_3_node_1_location[1], element_3_node_2_location[1]], 'm-', label='Element 3')
plt.plot([element_4_node_1_location[0], element_4_node_2_location[0]], [element_4_node_1_location[1], element_4_node_2_location[1]], 'c-', label='Element 4')
plt.title('4 Element Frame Structure')
plt.xlabel('X-axis (m)')
plt.ylabel('Y-axis (m)')
plt.legend()
plt.grid()
plt.axis('equal')
plt.show()
