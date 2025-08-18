import numpy as np
import matplotlib.pyplot as plt

class FrameAnalysis:
    def __init__(self):
        self.np_precision = 4
        np.set_printoptions(precision=self.np_precision, suppress=True)
    
    def parse_node_input(self, nodes_input):
        """
        Parse compact node input format into full node dictionaries.
        Input format: [x, y, 'support_type', 'label'] or [x, y, 'support_type']
        """
        nodes = []
        for i, node_data in enumerate(nodes_input):
            if len(node_data) == 3:
                x, y, support_type = node_data
                label = f'Node {i + 1}'
            elif len(node_data) == 4:
                x, y, support_type, label = node_data
            else:
                raise ValueError(f"Node {i + 1}: Expected format [x, y, 'support_type'] or [x, y, 'support_type', 'label']")
            
            nodes.append({
                'position': np.array([x, y]),
                'type': support_type,
                'label': label
            })
        
        return nodes
    def get_expected_dofs_for_node_type(self, node_type):
        """Return expected DOFs for each node type"""
        dof_types = {
            'fixed': [],  # No DOFs - all restrained
            'pin': ['translation'],  # Only rotational DOF (no translations)
            'free': ['translation', 'rotation']  # All DOFs available
        }
        return dof_types.get(node_type, ['translation', 'rotation'])
    
    def validate_node_dof_assignments(self, nodes, dof_definitions):
        """Validate that DOF assignments match node support types"""
        print("\nNODE TYPE VALIDATION:")
        print("="*60)
        
        validation_results = []
        
        for node_idx, node in enumerate(nodes):
            node_type = node['type']
            node_label = node.get('label', f'Node {node_idx + 1}')
            expected_dof_types = self.get_expected_dofs_for_node_type(node_type)
            
            # Find DOFs assigned to this node
            assigned_dofs = [dof for dof in dof_definitions if dof.get('node_index') == node_idx]
            assigned_dof_types = [dof['type'] for dof in assigned_dofs]
            
            print(f"\n{node_label} ({node_type.upper()} support):")
            print(f"  Position: ({node['position'][0]:.1f}, {node['position'][1]:.1f})")
            print(f"  Expected DOF types: {expected_dof_types if expected_dof_types else 'None (fully restrained)'}")
            print(f"  Assigned DOF types: {assigned_dof_types if assigned_dof_types else 'None'}")
            
            # Check for inconsistencies
            issues = []
            
            if node_type == 'fixed' and assigned_dofs:
                issues.append("FIXED node should have no DOFs assigned")
            
            elif node_type == 'pin':
                # Pin should only have rotation DOF, no translations
                if 'translation' in assigned_dof_types:
                    issues.append("PIN node should not have translation DOFs")
                if 'rotation' not in assigned_dof_types:
                    issues.append("PIN node should have rotation DOF")
            
            elif node_type == 'free':
                # Free node should have all DOF types
                if 'translation' not in assigned_dof_types:
                    issues.append("FREE node should have translation DOFs")
                if 'rotation' not in assigned_dof_types:
                    issues.append("FREE node should have rotation DOF")
            
            if issues:
                print(f"  ⚠️  VALIDATION ISSUES:")
                for issue in issues:
                    print(f"     - {issue}")
                validation_results.append({
                    'node': node_label,
                    'status': 'ISSUES',
                    'issues': issues
                })
            else:
                print(f"  ✅ Validation passed")
                validation_results.append({
                    'node': node_label,
                    'status': 'PASSED',
                    'issues': []
                })
        
        # Summary
        passed = sum(1 for result in validation_results if result['status'] == 'PASSED')
        total = len(validation_results)
        print(f"\nValidation Summary: {passed}/{total} nodes passed validation")
        
        if passed != total:
            print("\n⚠️  Please review the issues above to ensure your model is correct!")
        
        return validation_results
    
    def plot_node_supports(self, ax, nodes):
        """Add visual indicators for different support types"""
        support_symbols = {
            'fixed': {'marker': 's', 'color': 'red', 'size': 120, 'label': 'Fixed Support'},
            'pin': {'marker': '^', 'color': 'blue', 'size': 100, 'label': 'Pin Support'},
            'free': {'marker': 'o', 'color': 'green', 'size': 80, 'label': 'Free Node'}
        }
        
        # Track which support types we've plotted for legend
        plotted_types = set()
        
        for node in nodes:
            pos = node['position']
            node_type = node['type']
            symbol = support_symbols[node_type]
            
            # Plot support symbol
            ax.scatter(pos[0], pos[1], 
                      marker=symbol['marker'], 
                      c=symbol['color'], 
                      s=symbol['size'],
                      alpha=0.7, 
                      edgecolors='black',
                      linewidth=2,
                      zorder=15,
                      label=symbol['label'] if node_type not in plotted_types else "")
            
            plotted_types.add(node_type)
            
            # Add node label
            ax.text(pos[0], pos[1] - 0.15, node.get('label', ''), 
                   fontsize=10, ha='center', weight='bold')
    
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

    def extract_node_coordinates_and_dofs(self, nodes, dof_definitions, total_structural_deflections):
        """
        Extract node coordinates and their corresponding DOFs from nodes and dof_definitions.
        Returns a dictionary mapping node coordinates to their DOF information.
        """
        node_dofs = {}
        
        for node_idx, node in enumerate(nodes):
            node_pos = tuple(node['position'])  # Convert to tuple for dict key
            node_type = node['type']
            node_label = node.get('label', f'Node {node_idx + 1}')
            
            node_dofs[node_pos] = {
                'dofs': [],
                'type': node_type,
                'label': node_label,
                'index': node_idx
            }
        
        # Map DOFs to nodes based on node_index in dof_definitions
        for dof_idx, dof_def in enumerate(dof_definitions):
            if dof_idx < len(total_structural_deflections):
                dof_name = dof_def['name']
                dof_type = dof_def['type']
                dof_value = total_structural_deflections[dof_idx]
                node_index = dof_def.get('node_index', 0)
                
                # Find the node position for this node index
                if node_index < len(nodes):
                    node_pos = tuple(nodes[node_index]['position'])
                    unit = "[m]" if dof_type == 'translation' else "[rad]" if dof_type == 'rotation' else "[?]"
                    
                    node_dofs[node_pos]['dofs'].append({
                        'name': dof_name,
                        'value': dof_value,
                        'unit': unit,
                        'dof_idx': dof_idx
                    })
        
        return node_dofs

    def add_node_dof_textboxes(self, ax, node_dofs):
        """
        Add text boxes next to each node showing their DOF values and support type.
        """
        for node_coords, node_data in node_dofs.items():
            # Create text content
            text_lines = [f"{node_data['label']} ({node_data['type'].upper()})"]
            text_lines.append(f"Position: ({node_coords[0]:.1f}, {node_coords[1]:.1f})")
            text_lines.append("─" * 30)
            
            if node_data['dofs']:
                for dof in node_data['dofs']:
                    text_lines.append(f"{dof['name']}: {dof['value']:.3e} {dof['unit']}")
            else:
                text_lines.append("No DOFs (Fully Restrained)")
            
            text_content = "\n".join(text_lines)
            
            # Position text box offset from node
            offset_x = 0.4  # Horizontal offset
            offset_y = 0.3  # Vertical offset
            
            # Adjust offset based on node position to avoid overlaps
            if node_coords[0] > 1.0:  # Nodes on the right side
                offset_x = -1.5  # Move text to the left
            if node_coords[1] > 1.0:  # Nodes above middle
                offset_y = -0.7  # Move text below
            
            # Choose text box color based on node type
            colors = {
                'fixed': 'lightcoral',
                'pin': 'lightblue', 
                'free': 'lightgreen'
            }
            box_color = colors.get(node_data['type'], 'lightgray')
            
            # Add text box
            ax.text(node_coords[0] + offset_x, node_coords[1] + offset_y, 
                   text_content, fontsize=9, 
                   bbox=dict(boxstyle='round,pad=0.5', facecolor=box_color, alpha=0.9),
                   verticalalignment='top', fontfamily='monospace')

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
        nodes_input = config['nodes']
        nodes = self.parse_node_input(nodes_input)  # Parse compact format
        elements = config['elements']
        loads = config['loads']
        normal_forces = config['normal_forces']
        display_magnitude = config['display_magnitude']
        show_matrices = config['show_matrices']
        show_deflection_labels = config['show_deflection_labels']
        validate_node_types = config.get('validate_node_types', False)
        
        # Extract DOF definitions from config
        dof_definitions = config.get('dof_definitions', [])
        
        # Validate node type assignments if requested
        if validate_node_types and dof_definitions:
            self.validate_node_dof_assignments(nodes, dof_definitions)
        
        # Update elements with node positions from nodes array
        for elem in elements:
            node1_idx, node2_idx = elem['node_indices']
            elem['node_positions'] = [nodes[node1_idx]['position'], nodes[node2_idx]['position']]
        
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
        
        # Plot support symbols first (so they appear behind elements)
        self.plot_node_supports(ax, nodes)
        
        # Plot each element's deflected shape using total structural deflections
        for i, elem in enumerate(elements, 1):
            elem_key = f'element_{i}'
            node1_pos = elem['node_positions'][0]
            node2_pos = elem['node_positions'][1]
            
            self.plot_deflected_shape_with_total_deflections(
                ax, node1_pos[0], node1_pos[1], node2_pos[0], node2_pos[1],
                total_structural_deflections, elem['assembly_matrix'], 
                display_magnitude, 30, elem['angle'], f"Element {i}", False)
        
        # Extract node information and add DOF textboxes
        if show_deflection_labels and dof_definitions:
            node_dofs = self.extract_node_coordinates_and_dofs(nodes, dof_definitions, total_structural_deflections)
            self.add_node_dof_textboxes(ax, node_dofs)
            
            # Print to console as well
            print("\nNODE DOF INFORMATION:")
            print("="*50)
            for node_coords, node_data in node_dofs.items():
                print(f"\n{node_data['label']} at ({node_coords[0]:.1f}, {node_coords[1]:.1f}) - {node_data['type'].upper()} support:")
                if node_data['dofs']:
                    for dof in node_data['dofs']:
                        print(f"  {dof['name']:15s}: {dof['value']:12.6e} {dof['unit']}")
                else:
                    print("  No DOFs assigned (Fully restrained)")
            print()
        
        ax.set_aspect('equal')
        ax.grid(True, alpha=0.3)
        ax.legend()
        ax.set_xlabel('X (m)')
        ax.set_ylabel('Y (m)')
        ax.set_title(f'Frame Analysis with Node Support Types and DOF Information\n(Magnification Factor: {display_magnitude}x)')
        plt.tight_layout()
        plt.show()
        
        return results