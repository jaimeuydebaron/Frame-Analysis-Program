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

    def extract_node_coordinates_and_dofs(self, elements, dof_definitions, total_structural_deflections):
        """
        Extract unique node coordinates and their corresponding DOFs from elements.
        Returns a dictionary mapping node coordinates to their DOF information.
        """
        nodes = {}
        
        for elem in elements:
            node1_pos = tuple(elem['node_positions'][0])  # Convert to tuple for dict key
            node2_pos = tuple(elem['node_positions'][1])
            
            # Initialize nodes if not seen before
            if node1_pos not in nodes:
                nodes[node1_pos] = {'dofs': []}
            if node2_pos not in nodes:
                nodes[node2_pos] = {'dofs': []}
        
        # Map DOFs to nodes based on naming convention
        for dof_idx, dof_def in enumerate(dof_definitions):
            if dof_idx < len(total_structural_deflections):
                dof_name = dof_def['name']
                dof_type = dof_def['type']
                dof_value = total_structural_deflections[dof_idx]
                
                # Try to extract node information from DOF name
                # Expected format: "Node1_ux", "Node2_ry", etc.
                if 'Node' in dof_name:
                    try:
                        # Extract node number from DOF name
                        node_part = dof_name.split('_')[0]  # Gets "Node1", "Node2", etc.
                        node_num = int(node_part.replace('Node', ''))
                        
                        # Map node number to actual coordinates
                        node_coords = list(nodes.keys())
                        if node_num <= len(node_coords):
                            target_node = node_coords[node_num - 1]  # Node1 -> index 0, etc.
                            
                            unit = "[m]" if dof_type == 'translation' else "[rad]" if dof_type == 'rotation' else "[?]"
                            
                            nodes[target_node]['dofs'].append({
                                'name': dof_name,
                                'value': dof_value,
                                'unit': unit,
                                'dof_idx': dof_idx
                            })
                    except (ValueError, IndexError):
                        # If parsing fails, try to assign to first available node
                        if nodes:
                            first_node = list(nodes.keys())[0]
                            unit = "[m]" if dof_type == 'translation' else "[rad]" if dof_type == 'rotation' else "[?]"
                            nodes[first_node]['dofs'].append({
                                'name': dof_name,
                                'value': dof_value,
                                'unit': unit,
                                'dof_idx': dof_idx
                            })
        
        return nodes

    def add_node_dof_textboxes(self, ax, nodes):
        """
        Add text boxes next to each node showing their DOF values.
        """
        for node_coords, node_data in nodes.items():
            if node_data['dofs']:  # Only add textbox if node has DOFs
                # Create text content
                text_lines = [f"Node at ({node_coords[0]:.1f}, {node_coords[1]:.1f}):"]
                text_lines.append("─" * 25)
                
                for dof in node_data['dofs']:
                    text_lines.append(f"{dof['name']}: {dof['value']:.3e} {dof['unit']}")
                
                text_content = "\n".join(text_lines)
                
                # Position text box offset from node
                offset_x = 0.3  # Horizontal offset
                offset_y = 0.2  # Vertical offset
                
                # Adjust offset based on node position to avoid overlaps
                if node_coords[0] > 2.5:  # Nodes on the right side
                    offset_x = -1.2  # Move text to the left
                if node_coords[1] > 0:  # Nodes above origin
                    offset_y = -0.5  # Move text below
                
                # Add text box
                ax.text(node_coords[0] + offset_x, node_coords[1] + offset_y, 
                       text_content, fontsize=9, 
                       bbox=dict(boxstyle='round,pad=0.5', facecolor='lightgreen', alpha=0.8),
                       verticalalignment='top', fontfamily='monospace')
                
                # Add a small circle to mark the node
                ax.plot(node_coords[0], node_coords[1], 'ko', markersize=8, zorder=10)
                ax.plot(node_coords[0], node_coords[1], 'wo', markersize=6, zorder=11)

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
        show_matrices = config['show_matrices']
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
        print_func = self.print_matrix_vertical if show_matrices else self.print_matrix_standard
        
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
                display_magnitude, 30, elem['angle'], f"Element {i}", False)
        
        # NEW: Extract node information and add DOF textboxes
        if show_deflection_labels and dof_definitions:
            nodes = self.extract_node_coordinates_and_dofs(elements, dof_definitions, total_structural_deflections)
            self.add_node_dof_textboxes(ax, nodes)
            
            # Print to console as well
            print("\nNODE DOF INFORMATION:")
            print("="*50)
            for node_coords, node_data in nodes.items():
                print(f"\nNode at ({node_coords[0]:.1f}, {node_coords[1]:.1f}):")
                if node_data['dofs']:
                    for dof in node_data['dofs']:
                        print(f"  {dof['name']:15s}: {dof['value']:12.6e} {dof['unit']}")
                else:
                    print("  No DOFs assigned")
            print()
        
        ax.set_aspect('equal')
        ax.grid(True, alpha=0.3)
        ax.legend()
        ax.set_xlabel('X (m)')
        ax.set_ylabel('Y (m)')
        ax.set_title(f'Frame Deflection Analysis with Node DOF Information\n(Magnification Factor: {display_magnitude}x)')
        plt.tight_layout()
        plt.show()
        
        return results