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
    def build_UDL_global_forcing_vector(self, w_bar, L, assembly_matrix, transformation_matrix, direction):
        """
        Build global forcing vector for a UDL.
        
        Args:
            w_bar: Load intensity [N/m]
            L: Element length [m]
            assembly_matrix: Element assembly matrix
            transformation_matrix: Element transformation matrix
            direction: 'local' (perpendicular to element) or 'vertical' (vertically downward)
        
        Returns:
            tuple: (UDL_global_forcing_vector, equivalent_forcing_element_coordinates)
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
            angle_rad = - 1 * np.arctan2(transformation_matrix[1, 0], transformation_matrix[0, 0])

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
        return UDL_global_forcing_vector, equivalent_forcing_element_coordinates, big_F

    def build_LVL_peakXZERO_global_forcing_vector(self, w_bar, L, assembly_matrix, transformation_matrix):
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
        return LVL_peakXZERO_global_forcing_vector, equivalent_forcing_element_coordinates, big_F

    def build_LVL_peakXL_global_forcing_vector(self, w_bar, L, assembly_matrix, transformation_matrix):
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
        return LVL_peakXL_global_forcing_vector, equivalent_forcing_element_coordinates, big_F

    def build_PDL_global_forcing_vector(self, P, L, a, assembly_matrix, transformation_matrix, direction):
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
        return PDL_global_forcing_vector, equivalent_forcing_element_coordinates, big_F

    def build_DAL_global_forcing_vector(self, P, L, assembly_matrix, transformation_matrix):
        equivalent_forcing_element_coordinates = np.array([P * L / 2,
                                                        0,
                                                        0,
                                                        P * L / 2,
                                                        0,
                                                        0])
        big_F = np.transpose(transformation_matrix) @ equivalent_forcing_element_coordinates
        DAL_global_forcing_vector = assembly_matrix @ big_F
        return DAL_global_forcing_vector, equivalent_forcing_element_coordinates, big_F
    
    def build_CAL_global_forcing_vector(self, P, L, a, assembly_matrix, transformation_matrix):
        equivalent_forcing_element_coordinates = np.array([P * (1 - a / L),
                                                        0,
                                                        0,
                                                        P * a / L,
                                                        0,
                                                        0])
        big_F = np.transpose(transformation_matrix) @ equivalent_forcing_element_coordinates
        CAL_global_forcing_vector = assembly_matrix @ big_F
        return CAL_global_forcing_vector, equivalent_forcing_element_coordinates, big_F

    # NEW REACTION CALCULATION FUNCTIONS
    def calculate_reactions(self, config, results):
        """
        Calculate reaction forces at all constrained nodes.
        
        Reaction = Element Forces - Load Equivalent Forces
        """
        reactions = {}
        nodes = config['nodes']
        elements = config['elements']
        
        # Initialize reaction storage for all nodes
        for node_id, node_info in nodes.items():
            if self.is_constrained_node(node_info['type']):
                reactions[node_id] = {
                    'type': node_info['type'],
                    'coords': node_info['coords'],
                    'forces': np.zeros(3),  # [Rx, Ry, Mz]
                    'constrained_dofs': self.get_constrained_dofs(node_info['type'])
                }
        
        # For each constrained node, calculate reactions
        for node_id in reactions.keys():
            element_forces = np.zeros(3)  # Forces from connected elements
            load_forces = np.zeros(3)     # Forces from loads on connected elements
            
            # Find all elements connected to this node
            for elem_idx, elem in enumerate(elements, 1):
                elem_key = f'element_{elem_idx}'
                
                # Check if this element connects to the current node
                if elem['start_node'] == node_id:
                    # Node is at start of element (indices 0-2)
                    element_forces += results[elem_key]['global_forces'][:3]
                    
                    # Add load contributions from this element
                    for load_info in results[elem_key]['loads']:
                        load_forces += load_info['equiv_global'][:3]
                        
                elif elem['end_node'] == node_id:
                    # Node is at end of element (indices 3-5) 
                    element_forces += results[elem_key]['global_forces'][3:6]
                    
                    # Add load contributions from this element
                    for load_info in results[elem_key]['loads']:
                        load_forces += load_info['equiv_global'][3:6]
            
            # Calculate reaction: Element forces - Load equivalent forces
            reactions[node_id]['forces'] = element_forces - load_forces
            
            # Zero out reactions for unconstrained DOFs
            constrained_dofs = reactions[node_id]['constrained_dofs']
            for i in range(3):
                if not constrained_dofs[i]:
                    reactions[node_id]['forces'][i] = 0.0
        
        return reactions
    
    def is_constrained_node(self, node_type):
        """Check if a node has any constraints"""
        return node_type in ['fixed', 'pinned', 'roller_x', 'roller_y', 
                           'symmetrical_roller_x', 'symmetrical_roller_y']
    
    def get_constrained_dofs(self, node_type):
        """Get which DOFs are constrained for each node type"""
        constraints = {
            'fixed': [True, True, True],              # ux, uy, rz all constrained
            'pinned': [True, True, False],            # ux, uy constrained, rz free
            'roller_x': [False, True, False],         # uy constrained
            'roller_y': [True, False, False],         # ux constrained  
            'symmetrical_roller_x': [False, True, True],  # uy, rz constrained
            'symmetrical_roller_y': [True, False, True], # ux, rz constrained
            'free': [False, False, False]             # All free
        }
        return constraints.get(node_type, [False, False, False])

    def print_reactions(self, reactions):
        """Print reaction forces in a formatted table"""
        if not reactions:
            print("\nNo reaction forces to display (no constrained nodes)")
            return
        
        print(f"\n{'='*80}")
        print(f"{'REACTION FORCES AT CONSTRAINED NODES':^80}")
        print(f"{'='*80}")
        
        print(f"{'Node':<6} {'Type':<20} {'Rx [N]':<15} {'Ry [N]':<15} {'Mz [N⋅m]':<15}")
        print(f"{'-'*80}")
        
        for node_id, reaction_data in reactions.items():
            node_type = reaction_data['type']
            forces = reaction_data['forces']
            constrained = reaction_data['constrained_dofs']
            
            # Format values, show only constrained DOFs
            rx_str = f"{forces[0]:12.3f}" if constrained[0] else "---"
            ry_str = f"{forces[1]:12.3f}" if constrained[1] else "---" 
            mz_str = f"{forces[2]:12.3f}" if constrained[2] else "---"
            
            print(f"{node_id:<6} {node_type:<20} {rx_str:<15} {ry_str:<15} {mz_str:<15}")
        
        print(f"{'-'*80}\n")

    def add_reaction_arrows_to_plot(self, ax, reactions, scale_factor=1.0):
        """Add reaction force arrows to the plot"""
        if not reactions:
            return
        
        arrow_length = 0.5  # Fixed arrow length
        
        for node_id, reaction_data in reactions.items():
            coords = reaction_data['coords']
            forces = reaction_data['forces']
            constrained = reaction_data['constrained_dofs']
            
            x, y = coords[0], coords[1]
            
            # Draw force arrows for constrained DOFs only
            if constrained[0] and abs(forces[0]) > 1e-6:  # Rx
                # Use fixed arrow length, direction based on force sign
                dx = arrow_length if forces[0] > 0 else -arrow_length
                ax.arrow(x, y, dx, 0, head_width=0.1, head_length=0.05, 
                        fc='red', ec='red', linewidth=2, alpha=0.8)
                ax.text(x + dx/2, y + 0.15, f'Rx={forces[0]:.1f}N', 
                    ha='center', va='bottom', fontsize=8, color='red', weight='bold')
            
            if constrained[1] and abs(forces[1]) > 1e-6:  # Ry  
                # Use fixed arrow length, direction based on force sign
                dy = arrow_length if forces[1] > 0 else -arrow_length
                ax.arrow(x, y, 0, dy, head_width=0.1, head_length=0.05,
                        fc='blue', ec='blue', linewidth=2, alpha=0.8)
                ax.text(x - 0.15, y + dy/2, f'Ry={forces[1]:.1f}N',
                    ha='right', va='center', fontsize=8, color='blue', weight='bold',
                    rotation=90)
            
            if constrained[2] and abs(forces[2]) > 1e-6:  # Mz
                # Draw moment as a curved arrow 
                circle = plt.Circle((x, y), 0.2, fill=False, color='green', linewidth=2, alpha=0.8)
                ax.add_patch(circle)
                
                # Add moment magnitude text
                ax.text(x + 0.25, y, f'Mz={forces[2]:.1f}N⋅m',
                    ha='left', va='center', fontsize=8, color='green', weight='bold')

    # ENHANCED PRINTING FUNCTIONS
    def find_common_factor(self, matrix):
        """Find common factor to simplify matrix display using engineering increments (10^3)"""
        # Handle edge cases
        if matrix.size == 0:
            return 1
        
        non_zero_values = matrix[matrix != 0]
        if len(non_zero_values) == 0:
            return 1
        
        # Handle very small values that might be numerical noise
        abs_values = np.abs(non_zero_values)
        significant_values = abs_values[abs_values > 1e-15]
        if len(significant_values) == 0:
            return 1
        
        # Use median for more robust factor selection
        median_val = np.median(significant_values)
        
        # Engineering powers of 10^3: ..., 1e-6, 1e-3, 1e0, 1e3, 1e6, 1e9, ...
        # This gives us: micro, milli, base, kilo, mega, giga units
        engineering_powers = range(-6, 10)  # Covers 1e-6 to 1e9
        
        best_factor = 1
        best_score = float('inf')
        
        for power in engineering_powers:
            factor = 10**(power * 3)  # 10^(-6), 10^(-3), 10^0, 10^3, 10^6, etc.
            
            simplified = significant_values / factor
            
            # Calculate score: prefer factors that put most values in range [0.1, 999]
            # This gives nice readable numbers for engineering
            in_range = np.sum((simplified >= 0.1) & (simplified <= 999))
            total_values = len(simplified)
            
            # Bonus for values in ideal engineering range [1, 99]
            in_ideal_range = np.sum((simplified >= 1) & (simplified <= 99))
            
            # Score: lower is better
            score = (total_values - in_range) * 10 + (total_values - in_ideal_range)
            
            # Heavy penalty for factors that create numbers outside engineering comfort zone
            if np.any(simplified > 9999) or np.any(simplified < 0.001):
                score += 1000
            
            if score < best_score:
                best_score = score
                best_factor = factor
        
        return best_factor

    def get_engineering_units(self, base_units, factor):
        """Convert base units to engineering units based on the factor"""
        # Define unit mappings for different factors
        unit_mappings = {
            # Force units
            "[N]": {
                1e-3: "[mN]",
                1: "[N]", 
                1e3: "[kN]",
                1e6: "[MN]"
            },
            "[N/m]": {
                1e-3: "[mN/m]",
                1: "[N/m]",
                1e3: "[kN/m]", 
                1e6: "[MN/m]"
            },
            "[N⋅m]": {
                1e-3: "[mN⋅m]",
                1: "[N⋅m]",
                1e3: "[kN⋅m]",
                1e6: "[MN⋅m]"
            },
            "[N⋅m/rad]": {
                1e-3: "[mN⋅m/rad]",
                1: "[N⋅m/rad]",
                1e3: "[kN⋅m/rad]",
                1e6: "[MN⋅m/rad]"
            },
            "[N/m, N⋅m/rad]": {
                1e-3: "[mN/m, mN⋅m/rad]",
                1: "[N/m, N⋅m/rad]",
                1e3: "[kN/m, kN⋅m/rad]",
                1e6: "[MN/m, MN⋅m/rad]"
            },
            "[N, N⋅m]": {
                1e-3: "[mN, mN⋅m]",
                1: "[N, N⋅m]",
                1e3: "[kN, kN⋅m]",
                1e6: "[MN, MN⋅m]"
            },
            # Length units
            "[m]": {
                1e-6: "[μm]",
                1e-3: "[mm]",
                1: "[m]",
                1e3: "[km]"
            },
            "[m, rad]": {
                1e-6: "[μm, μrad]",
                1e-3: "[mm, mrad]", 
                1: "[m, rad]",
                1e3: "[km, krad]"
            },
            # Dimensionless
            "[-]": {
                1: "[-]"  # Transformation matrices stay dimensionless
            }
        }
        
        # Return the engineering unit if mapping exists, otherwise return original
        if base_units in unit_mappings and factor in unit_mappings[base_units]:
            return unit_mappings[base_units][factor]
        else:
            return base_units

    def print_matrix_vertical(self, name, matrix, units=""):
        factor = self.find_common_factor(matrix)
        engineering_units = self.get_engineering_units(units, factor)
        factor_text = f" (×{factor:.0e})" if factor != 1 else ""
        
        print(f"\n{name}{factor_text} {engineering_units}:")
        print("-" * (len(name) + len(factor_text) + len(engineering_units) + 1))
        
        simplified_matrix = matrix / factor
        
        if matrix.ndim == 1:
            # Print vectors vertically with indices
            for i, val in enumerate(simplified_matrix):
                print(f"  [{i + 1:2d}]: {val:12.4f}")
        else:
            # Print matrices in standard bracket format
            print(np.array_str(simplified_matrix, precision=4, suppress_small=True))
        print()

    def print_matrix_standard(self, name, matrix, units=""):
        factor = self.find_common_factor(matrix)
        engineering_units = self.get_engineering_units(units, factor)
        factor_text = f" (×{factor:.0e})" if factor != 1 else ""
        
        print(f"\n{name}{factor_text} {engineering_units}:")
        print("-" * (len(name) + len(factor_text) + len(engineering_units) + 1))
        
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
                
                # Determine if this DOF is constrained based on support type
                is_constrained = False
                if node_type == 'fixed':
                    is_constrained = True
                elif node_type == 'roller_x' and local_dof == 1:  # uy constrained
                    is_constrained = True
                elif node_type == 'roller_y' and local_dof == 0:  # ux constrained
                    is_constrained = True
                elif node_type == 'pinned' and local_dof in [0, 1]:  # ux, uy constrained
                    is_constrained = True
                elif node_type == 'symmetrical_roller_x' and local_dof in [1, 2]:  # NEW: uy, rz constrained, ux free
                    is_constrained = True
                elif node_type == 'symmetrical_roller_y' and local_dof in [0, 2]:  # NEW: ux, rz constrained, uy free
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

    def get_load_type_label(self, load_type):
        """Get the proper label for each load type"""
        load_labels = {
            'UDL': 'UDL',
            'LVL_peakXL': 'LVL (peak at x=L)',
            'LVL_peakXZERO': 'LVL (peak at x=0)',
            'PDL': 'PDL',
            'DAL': 'DAL',
            'CAL': 'CAL'
        }
        return load_labels.get(load_type, load_type)

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
        
        # Calculate element stiffness matrices
        results = {}
        total_global_forcing_vector = np.array(normal_forces, dtype=float)
        
        # Store load information for each element
        element_loads = {i: [] for i in range(1, len(elements) + 1)}
        
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
        
        # Add distributed loads and store forcing vectors
        all_load_vectors = []  # Store all load vectors for structure level printing
        
        for load in loads:
            elem_idx = load['element'] - 1
            elem = elements[elem_idx]
            elem_key = f'element_{load["element"]}'
            
            load_info = {
                'type': load['type'],
                'magnitude': load['magnitude'],
                'element': load['element']
            }
            
            if load['type'] == 'UDL':
                direction = load.get('direction', 'local')
                load_vector, equiv_local, equiv_global = self.build_UDL_global_forcing_vector(
                    load['magnitude'], elem['length'], 
                    elem['assembly_matrix'], results[elem_key]['transformation_matrix'], direction)
                load_info['description'] = f"Uniformly Distributed Load: {load['magnitude']} N/m (direction: {direction})"
                load_info['direction'] = direction
                
            elif load['type'] == 'LVL_peakXL':
                load_vector, equiv_local, equiv_global = self.build_LVL_peakXL_global_forcing_vector(
                    load['magnitude'], elem['length'],
                    elem['assembly_matrix'], results[elem_key]['transformation_matrix'])
                load_info['description'] = f"Linearly Varying Load (peak at x=L): {load['magnitude']} N/m"
                
            elif load['type'] == 'LVL_peakXZERO':
                load_vector, equiv_local, equiv_global = self.build_LVL_peakXZERO_global_forcing_vector(
                    load['magnitude'], elem['length'],
                    elem['assembly_matrix'], results[elem_key]['transformation_matrix'])
                load_info['description'] = f"Linearly Varying Load (peak at x=0): {load['magnitude']} N/m"
                
            elif load['type'] == 'PDL':
                direction = load.get('direction', 'local')
                a = load.get('a', elem['length']/2)  # Default to mid-span
                load_vector, equiv_local, equiv_global = self.build_PDL_global_forcing_vector(
                    load['magnitude'], elem['length'], a,
                    elem['assembly_matrix'], results[elem_key]['transformation_matrix'], direction)
                load_info['description'] = f"Point Load: {load['magnitude']} N at x={a} m"
                load_info['a'] = a
                
            elif load['type'] == 'DAL':
                load_vector, equiv_local, equiv_global = self.build_DAL_global_forcing_vector(
                    load['magnitude'], elem['length'],
                    elem['assembly_matrix'], results[elem_key]['transformation_matrix'])
                load_info['description'] = f"Distributed Axial Load: {load['magnitude']} N/m"
                
            elif load['type'] == 'CAL':
                a = load.get('a', elem['length']/2)  # Default to mid-span
                load_vector, equiv_local, equiv_global = self.build_CAL_global_forcing_vector(
                    load['magnitude'], elem['length'], a,
                    elem['assembly_matrix'], results[elem_key]['transformation_matrix'])
                load_info['description'] = f"Concentrated Axial Load: {load['magnitude']} N at x={a} m"
                load_info['a'] = a

            # Store the forcing vectors in the load info
            load_info['equiv_local'] = equiv_local
            load_info['equiv_global'] = equiv_global
            load_info['global_vector'] = load_vector
            
            # Add to element loads and all loads
            element_loads[load['element']].append(load_info)
            all_load_vectors.append(load_info)
            
            total_global_forcing_vector += load_vector
        
        # Store load information in results for printing
        for i in range(1, len(elements) + 1):
            elem_key = f'element_{i}'
            results[elem_key]['loads'] = element_loads[i]
        
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
        
        # NEW: Calculate reaction forces
        reactions = self.calculate_reactions(config, results)
        
        # Print results
        print_func = self.print_matrix_vertical if show_matrices else self.print_matrix_standard
        
        for i, elem in enumerate(elements, 1):
            elem_key = f'element_{i}'
            print(f"\n{'-'*60}")
            print(f"ELEMENT {i}")
            print(f"{'-'*60}")
            
            print_func(f"LOCAL stiffness matrix K_Local_{i}", results[elem_key]['local_stiffness'], "[N/m, N⋅m/rad]")
            print_func(f"GLOBAL (coordinates) stiffness matrix K_Global_{i}", results[elem_key]['global_stiffness'], "[N/m, N⋅m/rad]")
            print_func(f"Transformation matrix T_{i}", results[elem_key]['transformation_matrix'], "[-]")
            print_func(f"GLOBAL element deflections d_Global_{i}", results[elem_key]['global_deflections'], "[m, rad]")
            print_func(f"LOCAL element deflections d_Local_{i}", results[elem_key]['local_deflections'], "[m, rad]")
            print_func(f"Element nodal forces LOCAL F_Local_{i}", results[elem_key]['local_forces'], "[N, N⋅m]")
            print_func(f"Element nodal forces GLOBAL F_Global_{i}", results[elem_key]['global_forces'], "[N, N⋅m]")

            # Print load equivalent forcing vectors at the very end of each element section
            if results[elem_key]['loads']:
                print(f"\n{'EXTRA FORCING VECTORS FOR ELEMENT ' + str(i)}")
                print(f"{'-'*60}")
                for j, load_info in enumerate(results[elem_key]['loads'], 1):
                    load_label = self.get_load_type_label(load_info['type'])
                    print(f"\n{load_label} Load {j}: {load_info['description']}")
                    print_func(f"Equivalent LOCAL forces f_Local_equiv_{load_label}_{i}_{j}", 
                              load_info['equiv_local'], "[N, N⋅m]")
                    print_func(f"Equivalent GLOBAL forces F_Global_equiv_{load_label}_{i}_{j}", 
                              load_info['equiv_global'], "[N, N⋅m]")
                    print_func(f"GLOBAL forcing vector Q_Global_vector_{load_label}_{i}_{j}", 
                              load_info['global_vector'], "[N, N⋅m]")
        
        # Structure level results
        print(f"\n{'-'*60}")
        print("STRUCTURE LEVEL")
        print(f"{'-'*60}")
        
        # Print individual element contributions to global stiffness matrix
        for i, elem in enumerate(elements, 1):
            elem_key = f'element_{i}'
            print_func(f"Element {i} stiffness matrix contribution) [K]_{i}", 
                      results[elem_key]['global_stiffness_matrix'], "[N/m, N⋅m/rad]")
        
        print_func("TOTAL global stiffness matrix [K]_total", total_global_stiffness_matrix, "[N/m, N⋅m/rad]")
        
        # Print all load contributions at structure level
        if all_load_vectors:
            print(f"\n{'STRUCTURE LEVEL LOAD CONTRIBUTIONS':^60}")
            print(f"{'-'*60}")
            for j, load_info in enumerate(all_load_vectors, 1):
                load_label = self.get_load_type_label(load_info['type'])
                print(f"\n{load_label} Structure Load {j} (Element {load_info['element']}): {load_info['description']}")
                print_func(f"GLOBAL forcing vector contribution Q_{load_label}_{j}", 
                          load_info['global_vector'], "[N, N⋅m]")
        
        print_func("Total structural deflections d_total", total_structural_deflections, "[m, rad]")
        print_func("TOTAL global forcing vector Q_total", total_global_forcing_vector, "[N, N⋅m]")
        
        # NEW: Print reaction forces
        self.print_reactions(reactions)
        
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

        # NEW: Add reaction arrows to plot
        self.add_reaction_arrows_to_plot(ax, reactions, scale_factor=5000)

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
            elif n_type == 'symmetrical_roller_x':  # NEW: Symmetrical roller in X direction
                marker = '<'  # left triangle
                color = 'purple'
            elif n_type == 'symmetrical_roller_y':  # NEW: Symmetrical roller in Y direction
                marker = 'v'  # down triangle
                color = 'purple'
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
        ax.set_title(f'Frame Deflection Analysis with Node DOF Information and Reactions\n(Magnification Factor: {display_magnitude}x)')
        plt.tight_layout()
        plt.show()
        
        return results