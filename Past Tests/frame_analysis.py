import numpy as np
import matplotlib.pyplot as plt

class FrameAnalysis:
    def __init__(self):
        self.np_precision = 4
        np.set_printoptions(precision=self.np_precision, suppress=True)

    def _process_nodes(self, nodes):
        constraints = {}
        for node_id, node in nodes.items():
            n = (node_id - 1) * 3
            if node['type'] == 'fixed':
                constraints.update({n: 0.0, n+1: 0.0, n+2: 0.0})
            elif node['type'] == 'roller_x':
                constraints.update({n+1: 0.0})
            elif node['type'] == 'roller_y':
                constraints.update({n: 0.0})
            elif node['type'] == 'pinned':
                constraints.update({n: 0.0, n+1: 0.0})
            elif node['type'] == 'free':
                pass
            else:
                raise ValueError(f"Unknown node type: {node['type']}")
        return constraints
   
    
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
    def build_UDL_global_forcing_vector(self, w_bar, L, assembly_matrix, transformation_matrix, direction='local'):
        """
        Build global forcing vector for a UDL.
        
        Args:
            w_bar: Load intensity [N/m]
            L: Element length [m]
            assembly_matrix: Element assembly matrix
            transformation_matrix: Element transformation matrix
            direction: 'local' (perpendicular to element) or 'vertical' (vertically downward)
        """
        if direction == 'local':
            # Original implementation - load acts perpendicular to element
            equivalent_forcing_element_coordinates = np.array([0,
                                                            w_bar * L / 2,
                                                            w_bar * L ** 2 / 12,
                                                            0,
                                                            w_bar * L / 2,
                                                            -w_bar * L ** 2 / 12])
        elif direction == 'vertical':
            # Extract angle from transformation matrix
            angle_rad = np.arctan2(transformation_matrix[1, 0], transformation_matrix[0, 0])
            
            # Decompose vertical load into local components
            w_transverse = w_bar * np.cos(angle_rad)  # Perpendicular to element
            w_axial = w_bar * np.sin(angle_rad)       # Along element axis
            
            # Combine transverse and axial components
            equivalent_forcing_element_coordinates = np.array([w_axial * L / 2,           # Axial at node 1
                                                            w_transverse * L / 2,       # Transverse at node 1
                                                            w_transverse * L ** 2 / 12, # Moment at node 1
                                                            w_axial * L / 2,           # Axial at node 2
                                                            w_transverse * L / 2,       # Transverse at node 2
                                                            -w_transverse * L ** 2 / 12]) # Moment at node 2
        else:
            raise ValueError(f"Unknown direction: {direction}. Use 'local' or 'vertical'")
        
        big_F = np.transpose(transformation_matrix) @ equivalent_forcing_element_coordinates
        UDL_global_forcing_vector = assembly_matrix @ big_F
        return UDL_global_forcing_vector

    def build_LVL_peakXZERO_global_forcing_vector(self, w_bar, L, assembly_matrix, transformation_matrix, direction='local'):
        """
        Build global forcing vector for a linearly varying load (peak at x=0, zero at x=L).
        """
        # Original implementation
        equivalent_forcing_element_coordinates = np.array([0,
                                                        7 * w_bar * L / 20,
                                                        w_bar * L ** 2 / 20,
                                                        0,
                                                        3 * w_bar * L / 20,
                                                        -w_bar * L ** 2 / 30])
        
        big_F = np.transpose(transformation_matrix) @ equivalent_forcing_element_coordinates
        LVL_peakXZERO_global_forcing_vector = assembly_matrix @ big_F
        return LVL_peakXZERO_global_forcing_vector

    def build_LVL_peakXL_global_forcing_vector(self, w_bar, L, assembly_matrix, transformation_matrix, direction='local'):
        """
        Build global forcing vector for a linearly varying load (zero at x=0, peak at x=L).
        """
        # Original implementation
        equivalent_forcing_element_coordinates = np.array([0,
                                                        3 * w_bar * L / 20,
                                                        w_bar * L ** 2 / 30,
                                                        0,
                                                        7 * w_bar * L / 20,
                                                        -w_bar * L ** 2 / 20])
        
        big_F = np.transpose(transformation_matrix) @ equivalent_forcing_element_coordinates
        LVL_peakXL_global_forcing_vector = assembly_matrix @ big_F
        return LVL_peakXL_global_forcing_vector

    def build_PDL_global_forcing_vector(self, P, L, a, assembly_matrix, transformation_matrix, direction='local'):
        """
        Build global forcing vector for a point load at mid-span.
        """
        if direction == 'local':
            # Original implementation
            equivalent_forcing_element_coordinates = np.array([0,
                                                            P / 2,
                                                            P * L / 8,
                                                            0,
                                                            P / 2,
                                                            -P * L / 8])
        elif direction == 'vertical':
            # Extract angle from transformation matrix
            angle_rad = np.arctan2(transformation_matrix[1, 0], transformation_matrix[0, 0])
            
            # Decompose vertical load into local components
            P_transverse = P * np.cos(angle_rad)
            P_axial = P * np.sin(angle_rad)
            
            # Combine components
            equivalent_forcing_element_coordinates = np.array([P_axial * (1 - a / L),           # Half axial at node 1
                                                            P_transverse / 2,       # Half transverse at node 1
                                                            P_transverse * L / 8,   # Moment at node 1
                                                            P_axial * a/ L,           # Half axial at node 2
                                                            P_transverse / 2,       # Half transverse at node 2
                                                            -P_transverse * L / 8]) # Moment at node 2
        else:
            raise ValueError(f"Unknown direction: {direction}. Use 'local' or 'vertical'")
        
        big_F = np.transpose(transformation_matrix) @ equivalent_forcing_element_coordinates
        PDL_global_forcing_vector = assembly_matrix @ big_F
        return PDL_global_forcing_vector

    def build_DAL_global_forcing_vector(self, P, L, assembly_matrix, transformation_matrix):
        equivalent_forcing_element_coordinates = np.array([P * L / 2,
                                                        0,
                                                        0,
                                                        P * L / 2,
                                                        0,
                                                        0])
        big_F = np.transpose(transformation_matrix) @ equivalent_forcing_element_coordinates
        DAL_global_forcing_vector = assembly_matrix @ big_F
        return DAL_global_forcing_vector
    
    def build_CAL_global_forcing_vector(self, P, L, a, assembly_matrix, transformation_matrix):
        equivalent_forcing_element_coordinates = np.array([P * (1 - a / L),
                                                        0,
                                                        0,
                                                        P * a / L,
                                                        0,
                                                        0])
        big_F = np.transpose(transformation_matrix) @ equivalent_forcing_element_coordinates
        CAL_global_forcing_vector = assembly_matrix @ big_F
        return CAL_global_forcing_vector

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

    def extract_node_coordinates_and_dofs(self, nodes_config, total_structural_deflections):
        """
        Extract node coordinates and their corresponding DOFs using the actual node configuration.
        This handles reduced DOF systems where constrained DOFs are eliminated.
        """
        nodes = {}
        
        # Use the actual node configuration to get coordinates
        for node_id, node_info in nodes_config.items():
            node_coords = tuple(node_info['coords'])
            nodes[node_coords] = {'dofs': [], 'node_id': node_id}
        
        # Create a mapping of unconstrained DOFs to global DOF indices
        # This accounts for the reduced system where constrained DOFs are eliminated
        dof_names = ['ux', 'uy', 'rz']
        dof_types = ['translation', 'translation', 'rotation']
        
        # Create list of all DOFs with their constraint status
        all_dofs = []
        for node_id, node_info in nodes_config.items():
            node_type = node_info['type']
            for local_dof in range(3):
                dof_name = f"Node{node_id}_{dof_names[local_dof]}"
                
                # Determine if this DOF is constrained
                is_constrained = False
                if node_type == 'fixed':
                    is_constrained = True
                elif node_type == 'roller_x' and local_dof == 1:  # uy constrained
                    is_constrained = True
                elif node_type == 'roller_y' and local_dof == 0:  # ux constrained
                    is_constrained = True
                elif node_type == 'pinned' and local_dof in [0, 1]:  # ux, uy constrained
                    is_constrained = True
                
                all_dofs.append({
                    'node_id': node_id,
                    'local_dof': local_dof,
                    'name': dof_name,
                    'type': dof_types[local_dof],
                    'constrained': is_constrained
                })
        
        # Map unconstrained DOFs to solution indices
        unconstrained_dof_idx = 0
        for dof_info in all_dofs:
            node_id = dof_info['node_id']
            node_coords = None
            for nid, ninfo in nodes_config.items():
                if nid == node_id:
                    node_coords = tuple(ninfo['coords'])
                    break
            
            if dof_info['constrained']:
                # Constrained DOF - value is zero
                dof_value = 0.0
            else:
                # Unconstrained DOF - get value from solution
                if unconstrained_dof_idx < len(total_structural_deflections):
                    dof_value = total_structural_deflections[unconstrained_dof_idx]
                    unconstrained_dof_idx += 1
                else:
                    dof_value = 0.0
            
            # Add to node's DOF list
            unit = "[m]" if dof_info['type'] == 'translation' else "[rad]"
            nodes[node_coords]['dofs'].append({
                'name': dof_info['name'],
                'value': dof_value,
                'unit': unit,
                'constrained': dof_info['constrained']
            })
        
        return nodes

    def add_node_dof_textboxes(self, ax, nodes, nodes_config):
        """
        Add text boxes next to each node showing their DOF values.
        Skip fully fixed nodes, but show nodes with partial constraints or free nodes.
        """
        for node_coords, node_data in nodes.items():
            node_id = node_data['node_id']
            node_type = nodes_config[node_id]['type']
            
            if node_data['dofs']:  # Only consider nodes with DOFs
                # Skip fully fixed nodes (they have no interesting deflections to show)
                if node_type == 'fixed':
                    continue
                
                # For other node types, show textbox if there are any deflections (even small ones)
                # Use a much smaller tolerance for "non-zero" check
                non_zero_dofs = [dof for dof in node_data['dofs'] if not np.isclose(dof['value'], 0.0, atol=1e-15)]
                
                # Show all DOFs for non-fixed nodes, regardless of magnitude
                display_dofs = node_data['dofs']  # Show all DOFs, not just non-zero ones
                
                if not display_dofs:
                    continue

                # Create text content
                text_lines = [f"Node {node_id} ({node_type}):"]
                text_lines.append("─" * 25)
                
                for dof in display_dofs:
                    text_lines.append(f"{dof['name']}: {dof['value']:.3e} {dof['unit']}")
                
                text_content = "\n".join(text_lines)
                
                # Position text box offset from node
                offset_x = 0.3
                offset_y = 0.2
                if node_coords[0] > 2.5:
                    offset_x = -1.2
                if node_coords[1] > 0:
                    offset_y = -0.5
                
                ax.text(node_coords[0] + offset_x, node_coords[1] + offset_y, 
                    text_content, fontsize=9, 
                    bbox=dict(boxstyle='round,pad=0.5', facecolor='lightgreen', alpha=0.8),
                    verticalalignment='top', fontfamily='monospace')
                
                # Add a small circle to mark the node
                ax.plot(node_coords[0], node_coords[1], 'ko', markersize=8, zorder=10)
                ax.plot(node_coords[0], node_coords[1], 'wo', markersize=6, zorder=11)

    # PLOTTING WITH TOTAL STRUCTURAL DEFLECTIONS
    def plot_deflected_shape_with_total_deflections(self, ax, start_node_coords, end_node_coords, 
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
        
        # Extract node coordinates
        node1XG, node1YG = start_node_coords[0], start_node_coords[1]
        node2XG, node2YG = end_node_coords[0], end_node_coords[1]
        
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
        
        constraints = self._process_nodes(config['nodes'])
        
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
                # Expect 'a' to be provided in the load dictionary
                a = load.get('a', 0)
                load_vector = self.build_PDL_global_forcing_vector(
                    load['magnitude'], elem['length'], a,
                    elem['assembly_matrix'], results[elem_key]['transformation_matrix'])
            elif load['type'] == 'DAL':
                load_vector = self.build_DAL_global_forcing_vector(
                    load['magnitude'], elem['length'],
                    elem['assembly_matrix'], results[elem_key]['transformation_matrix'])
            elif load['type'] == 'CAL':
                # Expect 'a' to be provided in the load dictionary
                a = load.get('a', 0)
                load_vector = self.build_CAL_global_forcing_vector(
                    load['magnitude'], elem['length'], a,
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
            
            # Get node coordinates from the nodes configuration
            start_node_id = elem['start_node']
            end_node_id = elem['end_node']
            start_node_coords = config['nodes'][start_node_id]['coords']
            end_node_coords = config['nodes'][end_node_id]['coords']
            
            self.plot_deflected_shape_with_total_deflections(
                ax, start_node_coords, end_node_coords,
                total_structural_deflections, elem['assembly_matrix'], 
                display_magnitude, 30, elem['angle'], f"Element {i}", False)
            
        # Extract node deflections using simplified approach
        nodes_with_dofs = self.extract_node_coordinates_and_dofs(
            config['nodes'], 
            total_structural_deflections
        )

        # Add DOF textboxes to plot
        self.add_node_dof_textboxes(ax, nodes_with_dofs, config['nodes'])

        # Plot symbols for node types
        for node_id, node in config['nodes'].items():
            x, y = node['coords']
            n_type = node['type']
            
            # Choose marker style based on type
            if n_type == 'fixed':
                marker = 's'  # square
                color = 'red'
            elif n_type == 'pinned':
                marker = 'o'  # circle
                color = 'blue'
            elif n_type == 'roller_x' or n_type == 'roller_y':
                marker = 'D'  # diamond
                color = 'green'
            elif n_type == 'free':
                marker = '^'  # triangle
                color = 'orange'
            else:
                marker = 'x'
                color = 'black'
            
            ax.plot(x, y, marker=marker, color=color, markersize=12, label=f"Node {node_id} ({n_type})")

        # Remove duplicate legend entries
        handles, labels = ax.get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        ax.legend(by_label.values(), by_label.keys())
        
        ax.set_aspect('equal')
        ax.grid(True, alpha=0.3)
        ax.legend()
        ax.set_xlabel('X (m)')
        ax.set_ylabel('Y (m)')
        ax.set_title(f'Frame Deflection Analysis with Node DOF Information\n(Magnification Factor: {display_magnitude}x)')
        plt.tight_layout()
        plt.show()
        
        return results