import numpy as np
import matplotlib.pyplot as plt

class FrameAnalysis:
    def __init__(self):
        self.np_precision = 4
        np.set_printoptions(precision=self.np_precision, suppress=True)
    
    # CORE CALCULATION FUNCTIONS
    def calc_local_coordinates_element_stiffness_matrix(self, youngs_modulus, area, second_moment_of_area, length):
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

    def calc_global_coordinates_element_stiffness_matrix(self, local_stiffness, angle):
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

        # Clean up small values
        global_stiffness[np.abs(global_stiffness) < 1e-9] = 0
        transformation_matrix[np.abs(transformation_matrix) < 1e-9] = 0
        global_stiffness[np.isclose(global_stiffness, 0)] = 0
        transformation_matrix[np.isclose(transformation_matrix, 0)] = 0

        return global_stiffness, transformation_matrix

    def calc_global_stiffness_matrix(self, assembly_matrix, global_stiffness):
        return assembly_matrix @ global_stiffness @ np.transpose(assembly_matrix)

    def calc_total_global_stiffness_matrix(self, global_stiffness_matrix_1, global_stiffness_matrix_2):
        return global_stiffness_matrix_1 + global_stiffness_matrix_2

    def calc_structural_deflections(self, total_stiffness_global_coordinates, global_forcing_vector):
        return np.linalg.solve(total_stiffness_global_coordinates, global_forcing_vector)

    def calc_global_element_structural_deflections(self, total_structural_deflections, assembly_matrix):
        return np.transpose(assembly_matrix) @ total_structural_deflections

    def calc_local_element_structural_deflections(self, global_structural_deflections, transformation_matrix):
        return transformation_matrix @ global_structural_deflections

    def calc_local_element_nodal_forces(self, local_stiffness, local_structural_deflections):
        return local_stiffness @ local_structural_deflections

    def calc_global_element_nodal_forces(self, global_stiffness, global_structural_deflections):
        return global_stiffness @ global_structural_deflections

    # LOAD FUNCTIONS
    def build_UDL_global_forcing_vector(self, w_bar, L, assembly_matrix, transformation_matrix):
        equivalent_forcing_element_coordinates = np.array([0,
                                                        w_bar * L / 2,
                                                        w_bar * L ** 2 / 12,
                                                        0,
                                                        w_bar * L / 2,
                                                        -w_bar * L ** 2 / 12])
        big_F = np.transpose(transformation_matrix) @ equivalent_forcing_element_coordinates
        UDL_global_forcing_vector = assembly_matrix @ big_F
        return UDL_global_forcing_vector

    def build_LVL_peakXL_global_forcing_vector(self, w_bar, L, assembly_matrix, transformation_matrix):
        equivalent_forcing_element_coordinates = np.array([0,
                                                        3 * w_bar * L / 20,
                                                        w_bar * L ** 2 / 30,
                                                        0,
                                                        7 * w_bar * L / 20,
                                                        - w_bar * L ** 2 / 20])
        big_F = np.transpose(transformation_matrix) @ equivalent_forcing_element_coordinates
        LVL_peakXL_global_forcing_vector = assembly_matrix @ big_F
        return LVL_peakXL_global_forcing_vector

    def build_LVL_peakXZERO_global_forcing_vector(self, w_bar, L, assembly_matrix, transformation_matrix):
        equivalent_forcing_element_coordinates = np.array([0,
                                                        7 * w_bar * L / 20,
                                                        w_bar * L ** 2 / 20,
                                                        0,
                                                        3 * w_bar * L / 20,
                                                        - w_bar * L ** 2 / 30])
        big_F = np.transpose(transformation_matrix) @ equivalent_forcing_element_coordinates
        LVL_peakXZERO_global_forcing_vector = assembly_matrix @ big_F
        return LVL_peakXZERO_global_forcing_vector

    def build_PDL_global_forcing_vector(self, P, L, assembly_matrix, transformation_matrix):
        equivalent_forcing_element_coordinates = np.array([0,
                                                        P / 2,
                                                        P * L / 8,
                                                        0,
                                                        P / 2,
                                                        -P * L / 8])
        big_F = np.transpose(transformation_matrix) @ equivalent_forcing_element_coordinates
        PDL_global_forcing_vector = assembly_matrix @ big_F
        return PDL_global_forcing_vector

    # ENHANCED PRINTING FUNCTIONS
    def find_common_factor(self, matrix):
        """Find common factor to simplify matrix display"""
        non_zero_values = matrix[matrix != 0]
        if len(non_zero_values) == 0:
            return 1
        
        # Check for common powers of 10
        abs_values = np.abs(non_zero_values)
        min_abs = np.min(abs_values)
        
        # Find the largest power of 10 that divides into most values
        for power in range(-12, 13):
            factor = 10**power
            if min_abs >= factor:
                # Check if this factor simplifies most values significantly
                simplified = non_zero_values / factor
                if np.all(np.abs(simplified) >= 0.1) and np.all(np.abs(simplified) < 10000):
                    return factor
        return 1

    def print_matrix_vertical(self, name, matrix, units=""):
        factor = self.find_common_factor(matrix)
        factor_text = f" (×{factor:.0e})" if factor != 1 else ""
        
        print(f"\n{name}{factor_text} {units}:")
        print("-" * (len(name) + len(factor_text) + len(units) + 1))
        
        simplified_matrix = matrix / factor
        
        if matrix.ndim == 1:
            # Print vectors vertically with indices
            for i, val in enumerate(simplified_matrix):
                print(f"  [{i:2d}]: {val:12.4f}")
        else:
            # Print matrices in standard bracket format
            print(np.array_str(simplified_matrix, precision=4, suppress_small=True))
        print()

    def print_matrix_standard(self, name, matrix, units=""):
        factor = self.find_common_factor(matrix)
        factor_text = f" (×{factor:.0e})" if factor != 1 else ""
        
        print(f"\n{name}{factor_text} {units}:")
        print("-" * (len(name) + len(factor_text) + len(units) + 1))
        
        simplified_matrix = matrix / factor
        
        if matrix.ndim == 1:
            # Print vectors vertically with indices
            for i, val in enumerate(simplified_matrix):
                print(f"  [{i:2d}]: {val:12.4f}")
        else:
            # Print matrices in standard bracket format
            print(np.array_str(simplified_matrix, precision=4, suppress_small=True))
        print()

    # PLOTTING WITH TOTAL STRUCTURAL DEFLECTIONS
    def plot_deflected_shape_with_total_deflections(self, ax, node1XG, node1YG, node2XG, node2YG, 
                                                  total_structural_deflections, assembly_matrix, 
                                                  disp_mag, N_points, angle, element_name, show_values=True):
        # Extract element deflections from total structural deflections using assembly matrix
        element_deflections_global = np.transpose(assembly_matrix) @ total_structural_deflections
        
        # Transform to local coordinates
        angle_rad = np.radians(angle)
        c = np.cos(angle_rad)
        s = np.sin(angle_rad)
        transformation_matrix = np.array([[c, s, 0, 0, 0, 0],
                                         [-s, c, 0, 0, 0, 0],
                                         [0, 0, 1, 0, 0, 0],
                                         [0, 0, 0, c, s, 0],
                                         [0, 0, 0, -s, c, 0],
                                         [0, 0, 0, 0, 0, 1]])
        
        element_deflections_local = transformation_matrix @ element_deflections_global
        
        # Now plot the deflected shape using these deflections
        L = np.sqrt((node2XG - node1XG) ** 2 + (node2YG - node1YG) ** 2)
        x_e = np.linspace(0, L, N_points)

        # Shape functions
        axial_shape_func_1 = (1 - x_e / L)
        axial_shape_func_2 = x_e / L
        transverse_shape_func_1 = 1 - 3 * x_e ** 2 / L ** 2 + 2 * x_e ** 3 / L ** 3
        transverse_shape_func_2 = x_e ** 3 / L ** 2 - 2 * x_e ** 2 / L + x_e
        transverse_shape_func_3 = 3 * x_e ** 2 / L ** 2 - 2 * x_e ** 3 / L ** 3
        transverse_shape_func_4 = x_e ** 3 / L ** 2 - x_e ** 2 / L

        # Displacements in local coordinates using total structural deflections
        d_e = element_deflections_local
        axial_displacement = axial_shape_func_1 * d_e[0] + axial_shape_func_2 * d_e[3]
        transverse_displacement = (transverse_shape_func_1 * d_e[1] +
                                    transverse_shape_func_2 * d_e[2] +
                                    transverse_shape_func_3 * d_e[4] +
                                    transverse_shape_func_4 * d_e[5])

        # Rotate back to global coordinates
        deflections_XG = axial_displacement * np.cos(angle_rad) - transverse_displacement * np.sin(angle_rad)
        deflections_YG = axial_displacement * np.sin(angle_rad) + transverse_displacement * np.cos(angle_rad)

        # Undeformed coordinates
        undeflected_XG = np.linspace(node1XG, node2XG, N_points)
        undeflected_YG = np.linspace(node1YG, node2YG, N_points)

        # Deflected coordinates
        deflected_XG = undeflected_XG + disp_mag * deflections_XG
        deflected_YG = undeflected_YG + disp_mag * deflections_YG

        # Plot lines
        ax.plot(undeflected_XG, undeflected_YG, 'k.-', linewidth=2, 
               label='Undeformed Shape' if element_name == "Element 1" else "")
        ax.plot(deflected_XG, deflected_YG, 'r.-', linewidth=2, 
               label='Deflected Shape' if element_name == "Element 1" else "")

        if show_values:
            # Add deflection labels at key points
            key_points = [0, N_points//4, N_points//2, 3*N_points//4, N_points-1]
            for i in key_points:
                if i < len(deflected_XG):
                    # Calculate actual deflection magnitude
                    actual_deflection_x = deflections_XG[i]
                    actual_deflection_y = deflections_YG[i]
                    deflection_magnitude = np.sqrt(actual_deflection_x**2 + actual_deflection_y**2)
                    
                    if deflection_magnitude > 1e-6:  # Only label significant deflections
                        ax.annotate(f'{deflection_magnitude:.2e}m', 
                                  (deflected_XG[i], deflected_YG[i]),
                                  xytext=(5, 5), textcoords='offset points',
                                  fontsize=8, color='blue',
                                  bbox=dict(boxstyle='round,pad=0.3', facecolor='yellow', alpha=0.7))

        # Add element label
        mid_idx = N_points // 2
        ax.text(undeflected_XG[mid_idx], undeflected_YG[mid_idx] + 0.1, element_name, 
               fontsize=10, ha='center', weight='bold')

    def run_analysis(self, config):
        """Main analysis function with configuration dictionary"""
        
        print("\n" + "="*100)
        print("                        INTERACTIVE FRAME ANALYSIS RESULTS")
        print("="*100)
        
        # Extract configuration
        elements = config['elements']
        loads = config['loads']
        normal_forces = config['normal_forces']
        display_magnitude = config['display_magnitude']
        show_vertical_matrices = config['show_vertical_matrices']
        show_deflection_labels = config['show_deflection_labels']
        
        # NEW: Extract DOF definitions from config
        dof_definitions = config.get('dof_definitions', [])
        
        # Calculate element stiffness matrices
        results = {}
        total_global_forcing_vector = np.array(normal_forces, dtype=float)
        
        for i, elem in enumerate(elements, 1):
            elem_key = f'element_{i}'
            results[elem_key] = {}
            
            # Local stiffness
            local_stiffness = self.calc_local_coordinates_element_stiffness_matrix(
                elem['youngs_modulus'], elem['area'], elem['second_moment_of_area'], elem['length'])
            results[elem_key]['local_stiffness'] = local_stiffness
            
            # Global stiffness and transformation
            global_stiffness, transformation_matrix = self.calc_global_coordinates_element_stiffness_matrix(
                local_stiffness, elem['angle'])
            results[elem_key]['global_stiffness'] = global_stiffness
            results[elem_key]['transformation_matrix'] = transformation_matrix
            
            # Global stiffness matrix contribution
            global_stiffness_matrix = self.calc_global_stiffness_matrix(elem['assembly_matrix'], global_stiffness)
            results[elem_key]['global_stiffness_matrix'] = global_stiffness_matrix
            
            if i == 1:
                total_global_stiffness_matrix = global_stiffness_matrix.copy()
            else:
                total_global_stiffness_matrix += global_stiffness_matrix
        
        # Add distributed loads
        for load in loads:
            elem_idx = load['element'] - 1
            elem = elements[elem_idx]
            elem_key = f'element_{load["element"]}'
            
            if load['type'] == 'UDL':
                load_vector = self.build_UDL_global_forcing_vector(
                    load['magnitude'], elem['length'], 
                    elem['assembly_matrix'], results[elem_key]['transformation_matrix'])
            elif load['type'] == 'LVL_peakXL':
                load_vector = self.build_LVL_peakXL_global_forcing_vector(
                    load['magnitude'], elem['length'],
                    elem['assembly_matrix'], results[elem_key]['transformation_matrix'])
            elif load['type'] == 'LVL_peakXZERO':
                load_vector = self.build_LVL_peakXZERO_global_forcing_vector(
                    load['magnitude'], elem['length'],
                    elem['assembly_matrix'], results[elem_key]['transformation_matrix'])
            elif load['type'] == 'PDL':
                load_vector = self.build_PDL_global_forcing_vector(
                    load['magnitude'], elem['length'],
                    elem['assembly_matrix'], results[elem_key]['transformation_matrix'])
            
            total_global_forcing_vector += load_vector
        
        # Solve for deflections
        total_structural_deflections = self.calc_structural_deflections(
            total_global_stiffness_matrix, total_global_forcing_vector)
        
        # Calculate element deflections and forces
        for i, elem in enumerate(elements, 1):
            elem_key = f'element_{i}'
            
            # Global deflections
            global_deflections = self.calc_global_element_structural_deflections(
                total_structural_deflections, elem['assembly_matrix'])
            results[elem_key]['global_deflections'] = global_deflections
            
            # Local deflections
            local_deflections = self.calc_local_element_structural_deflections(
                global_deflections, results[elem_key]['transformation_matrix'])
            results[elem_key]['local_deflections'] = local_deflections
            
            # Forces
            local_forces = self.calc_local_element_nodal_forces(
                results[elem_key]['local_stiffness'], local_deflections)
            global_forces = self.calc_global_element_nodal_forces(
                results[elem_key]['global_stiffness'], global_deflections)
            results[elem_key]['local_forces'] = local_forces
            results[elem_key]['global_forces'] = global_forces
        
        # Print results
        print_func = self.print_matrix_vertical if show_vertical_matrices else self.print_matrix_standard
        
        for i, elem in enumerate(elements, 1):
            elem_key = f'element_{i}'
            print(f"\n{'-'*60}")
            print(f"ELEMENT {i}")
            print(f"{'-'*60}")
            
            print_func(f"Local stiffness matrix [k]_local_{i}", results[elem_key]['local_stiffness'], "[N/m, N⋅m/rad]")
            print_func(f"Global stiffness matrix [k]_global_{i}", results[elem_key]['global_stiffness'], "[N/m, N⋅m/rad]")
            print_func(f"Transformation matrix [T]_{i}", results[elem_key]['transformation_matrix'], "[-]")
            print_func(f"Global element deflections {{d}}_global_{i}", results[elem_key]['global_deflections'], "[m, rad]")
            print_func(f"Local element deflections {{d}}_local_{i}", results[elem_key]['local_deflections'], "[m, rad]")
            print_func(f"Element nodal forces LOCAL {{F}}_local_{i}", results[elem_key]['local_forces'], "[N, N⋅m]")
            print_func(f"Element nodal forces GLOBAL {{F}}_global_{i}", results[elem_key]['global_forces'], "[N, N⋅m]")
        
        # Structure level results
        print(f"\n{'-'*60}")
        print("STRUCTURE LEVEL")
        print(f"{'-'*60}")
        
        # Print individual element contributions to global stiffness matrix
        for i, elem in enumerate(elements, 1):
            elem_key = f'element_{i}'
            print_func(f"Global stiffness matrix (Element {i} contribution) [K]_{i}", 
                      results[elem_key]['global_stiffness_matrix'], "[N/m, N⋅m/rad]")
        
        print_func("TOTAL global stiffness matrix [K]_total", total_global_stiffness_matrix, "[N/m, N⋅m/rad]")
        print_func("Total structural deflections {d}_total", total_structural_deflections, "[m, rad]")
        print_func("Global forcing vector {F}", total_global_forcing_vector, "[N, N⋅m]")
        
        print("\n" + "="*100)
        print("                                END OF RESULTS")
        print("="*100 + "\n")
        
        # Plot results
        fig, ax = plt.subplots(figsize=(14, 10))
        
        # Plot each element's deflected shape using total structural deflections
        for i, elem in enumerate(elements, 1):
            elem_key = f'element_{i}'
            node1_pos = elem['node_positions'][0]
            node2_pos = elem['node_positions'][1]
            
            self.plot_deflected_shape_with_total_deflections(
                ax, node1_pos[0], node1_pos[1], node2_pos[0], node2_pos[1],
                total_structural_deflections, elem['assembly_matrix'], 
                display_magnitude, 30, elem['angle'], f"Element {i}", False)  # Don't show individual labels
        
        # Create text box with total structural deflections using correct units
        if show_deflection_labels:
            # Use DOF definitions if provided, otherwise create fallback labels
            if dof_definitions:
                # Ensure we have enough DOF definitions
                while len(dof_definitions) < len(total_structural_deflections):
                    dof_definitions.append({'name': f'DOF_{len(dof_definitions)}', 'type': 'unknown'})
                
                # Create deflection text with correct units
                deflection_text = "TOTAL STRUCTURAL DEFLECTIONS:\n" + "="*40 + "\n"
                for i, (deflection, dof_def) in enumerate(zip(total_structural_deflections, dof_definitions)):
                    # Determine unit based on DOF type
                    if dof_def['type'] == 'translation':
                        unit = "[m]"
                    elif dof_def['type'] == 'rotation':
                        unit = "[rad]"
                    else:
                        unit = "[?]"
                    
                    deflection_text += f"DOF {i:2d} ({dof_def['name']:12s}): {deflection:10.4e} {unit}\n"
            else:
                # Fallback: try to guess based on magnitude (this is less reliable)
                deflection_text = "TOTAL STRUCTURAL DEFLECTIONS:\n" + "="*40 + "\n"
                deflection_text += "Warning: No DOF definitions provided. Units may be incorrect.\n"
                for i, deflection in enumerate(total_structural_deflections):
                    # Simple heuristic: smaller values might be rotations, larger might be translations
                    # This is not reliable and should not be used in practice
                    if abs(deflection) < 1e-3:
                        unit = "[rad?]"
                        dof_name = f"DOF_{i}_rot?"
                    else:
                        unit = "[m?]"
                        dof_name = f"DOF_{i}_trans?"
                    
                    deflection_text += f"DOF {i:2d} ({dof_name:12s}): {deflection:10.4e} {unit}\n"
            
            # Add text box to plot
            ax.text(0.02, 0.98, deflection_text, transform=ax.transAxes, fontsize=10,
                   verticalalignment='top', bbox=dict(boxstyle='round,pad=0.5', 
                   facecolor='lightblue', alpha=0.8), fontfamily='monospace')
            
            # Print to console as well
            print("\nTOTAL STRUCTURAL DEFLECTIONS:")
            print("="*50)
            if dof_definitions:
                for i, (deflection, dof_def) in enumerate(zip(total_structural_deflections, dof_definitions)):
                    if dof_def['type'] == 'translation':
                        unit = "[m]"
                    elif dof_def['type'] == 'rotation':
                        unit = "[rad]"
                    else:
                        unit = "[?]"
                    
                    print(f"  DOF {i:2d} ({dof_def['name']:12s}): {deflection:12.6e} {unit}")
            else:
                print("  Warning: No DOF definitions provided. Units may be incorrect.")
                for i, deflection in enumerate(total_structural_deflections):
                    if abs(deflection) < 1e-3:
                        unit = "[rad?]"
                        dof_name = f"DOF_{i}_rot?"
                    else:
                        unit = "[m?]"
                        dof_name = f"DOF_{i}_trans?"
                    
                    print(f"  DOF {i:2d} ({dof_name:12s}): {deflection:12.6e} {unit}")
            print()
        
        ax.set_aspect('equal')
        ax.grid(True, alpha=0.3)
        ax.legend()
        ax.set_xlabel('X (m)')
        ax.set_ylabel('Y (m)')
        ax.set_title(f'Frame Deflection Analysis with Total Structural Deflections\n(Magnification Factor: {display_magnitude}x)')
        plt.tight_layout()
        plt.show()
        
        return results


# USAGE EXAMPLE - CONFIGURE YOUR ANALYSIS HERE
def main():
    analyzer = FrameAnalysis()
    
    # CONFIGURATION - EDIT THIS SECTION FOR YOUR PROBLEM
    config = {
        # Element properties (you still define length, angle, assembly matrix)
        'elements': [
            {  # Element 1
                'length': 4,                    # m
                'second_moment_of_area': 5e-6,  # m^4
                'area': 2e-4,                   # m^2
                'youngs_modulus': 200e9,        # Pa
                'angle': 0,                     # degrees
                'node_positions': [np.array([0, 0]), np.array([4, 0])],  # Global coordinates
                'assembly_matrix': np.array([[0, 0, 0, 1, 0, 0],
                                             [0, 0, 0, 0, 0, 1],
                                             [0, 0, 0, 0, 0, 0],
                                             [0, 0, 0, 0, 0, 0],
                                             [0, 0, 0, 0, 0, 0],])
            },
            {  # Element 2
                'length': 5,                    # m
                'second_moment_of_area': 5e-6,  # m^4
                'area': 2e-4,                   # m^2
                'youngs_modulus': 200e9,        # Pa
                'angle': 36.87,                     # degrees
                'node_positions': [np.array([0, 0]), np.array([4, 3])],  # Global coordinates
                'assembly_matrix': np.array([[0, 0, 0, 0, 0, 0],
                                             [0, 0, 0, 0, 0, 0],   
                                             [0, 0, 0, 1, 0, 0],
                                             [0, 0, 0, 0, 1, 0],
                                             [0, 0, 0, 0, 0, 1]])
            },
            {  # Element 3
                'length': 3,                    # m
                'second_moment_of_area': 5e-6,  # m^4
                'area': 2e-4,                   # m^2
                'youngs_modulus': 200e9,        # Pa
                'angle': 90,                     # degrees
                'node_positions': [np.array([4, 0]), np.array([4, 3])],  # Global coordinates
                'assembly_matrix': np.array([[1, 0, 0, 0, 0, 0],
                                             [0, 0, 1, 0, 0, 0],   
                                             [0, 0, 0, 1, 0, 0],
                                             [0, 0, 0, 0, 1, 0],
                                             [0, 0, 0, 0, 0, 1]])
            }
        ],
        
        # NEW: DOF definitions - specify what each DOF represents
        'dof_definitions': [
            {'name': 'Node2_ux', 'type': 'translation'},     # DOF 0: Node 2 x-displacement
            {'name': 'Node1_rz', 'type': 'rotation'},        # DOF 1: Node 1 z-rotation  
            {'name': 'Node2_uy', 'type': 'translation'},     # DOF 2: Node 2 y-displacement
            {'name': 'Node2_uz', 'type': 'translation'},        # DOF 3: Node 2 z-rotation
            {'name': 'Node3_rz', 'type': 'rotation'}         # DOF 4: Node 3 z-rotation
        ],
        
        # Distributed loads (UDL, LVL_peakXL, LVL_peakXZERO, PDL)
        'loads': [
            {'type': 'UDL', 'element': 1, 'magnitude': -40000}  # N/m on element 1
            # Add more loads as needed:
            # {'type': 'PDL', 'element': 2, 'magnitude': 20000},  # N on element 2
            # {'type': 'LVL_peakXL', 'element': 1, 'magnitude': 10000}  # N/m peak at x=L
        ],
        
        # Normal nodal forces [F1x, F1y, M1z, F2x, F2y] - adjust size based on your DOFs
        'normal_forces': [0, 0, 25000, -10000, 0],  # N, N⋅m
        
        # Display options
        'display_magnitude': 100,        # Magnification factor for deflection plot
        'show_vertical_matrices': True,  # True for vertical matrix display, False for standard
        'show_deflection_labels': True   # True to show deflection values on plot
    }
    
    # Run the analysis
    results = analyzer.run_analysis(config)
    
    return results

if __name__ == "__main__":
    results = main()