import numpy as np
import matplotlib.pyplot as plt

class BarAnalysis:
    def __init__(self):
        self.np_precision = 4
        np.set_printoptions(precision=self.np_precision, suppress=True)
   
    # CORE CALCULATION FUNCTIONS FOR 2D BAR ELEMENTS
    def calc_local_coordinates_element_stiffness_matrix(self, youngs_modulus, area, length):
        """Calculate local stiffness matrix for 2D bar element (axial only)"""
        local_stiffness = (youngs_modulus * area) / length * np.array([[1, -1], [-1, 1]])
        return local_stiffness

    def calc_global_coordinates_element_stiffness_matrix(self, local_stiffness, angle):
        """Transform 2D bar stiffness matrix from local to global coordinates"""
        angle_rad = np.radians(angle)
        c = np.cos(angle_rad)
        s = np.sin(angle_rad)
        
        # Transformation matrix for 2D bar (2x4: local to global DOFs)
        transformation_matrix = np.array([[c, s, 0, 0],
                                         [0, 0, c, s]])
        
        # Global stiffness: T^T * K_local * T (4x4 matrix)
        global_stiffness = np.transpose(transformation_matrix) @ local_stiffness @ transformation_matrix

        # Clean up small values
        global_stiffness[np.abs(global_stiffness) < 1e-9] = 0
        transformation_matrix[np.abs(transformation_matrix) < 1e-9] = 0
        global_stiffness[np.isclose(global_stiffness, 0)] = 0
        transformation_matrix[np.isclose(transformation_matrix, 0)] = 0

        return global_stiffness, transformation_matrix

    def calc_global_stiffness_matrix(self, assembly_matrix, global_stiffness):
        """Assemble element stiffness matrix into structural stiffness matrix"""
        return assembly_matrix @ global_stiffness @ np.transpose(assembly_matrix)

    def calc_total_global_stiffness_matrix(self, global_stiffness_matrix_1, global_stiffness_matrix_2):
        """Sum multiple element stiffness matrices"""
        return global_stiffness_matrix_1 + global_stiffness_matrix_2

    def calc_structural_deflections(self, total_stiffness_global_coordinates, global_forcing_vector):
        """Solve for structural deflections"""
        return np.linalg.solve(total_stiffness_global_coordinates, global_forcing_vector)

    def calc_global_element_structural_deflections(self, total_structural_deflections, assembly_matrix):
        """Extract element deflections from structural deflections"""
        return np.transpose(assembly_matrix) @ total_structural_deflections

    def calc_local_element_structural_deflections(self, global_structural_deflections, transformation_matrix):
        """Transform global element deflections to local coordinates"""
        return transformation_matrix @ global_structural_deflections

    def calc_local_element_nodal_forces(self, local_stiffness, local_structural_deflections):
        """Calculate local element forces"""
        return local_stiffness @ local_structural_deflections

    def calc_global_element_nodal_forces(self, global_stiffness, global_structural_deflections):
        """Calculate global element forces"""
        return global_stiffness @ global_structural_deflections

    def calc_element_normal_stress(self, local_forces, area):
        """Calculate normal stress in 2D bar element using sigma = F/A
        
        Args:
            local_forces: Array of local nodal forces [F1, F2] where F1 is axial force at node 1
            area: Cross-sectional area of the element
            
        Returns:
            normal_stress: Normal stress in the element (positive = tension, negative = compression)
        """
        # For a bar element, the axial force is the same throughout the element
        # Use the force at the first node (with appropriate sign convention)
        axial_force = local_forces[0]  # Force at node 1
        
        # Normal stress = Force / Area
        normal_stress = axial_force / area
        
        return normal_stress

    # LOAD FUNCTIONS FOR 2D BAR ELEMENTS
    def build_DAL_global_forcing_vector(self, w_bar, L, assembly_matrix, transformation_matrix):
        """Build global forcing vector for distributed axial load"""
        # Equivalent nodal forces in local coordinates (axial only)
        equivalent_forcing_element_coordinates = np.array([w_bar * L / 2, w_bar * L / 2])
        
        # Transform to global coordinates
        big_F = np.transpose(transformation_matrix) @ equivalent_forcing_element_coordinates
        
        # Assemble into structural force vector
        DAL_global_forcing_vector = assembly_matrix @ big_F
        return DAL_global_forcing_vector, equivalent_forcing_element_coordinates, big_F
    
    def build_CAL_global_forcing_vector(self, P, L, a, assembly_matrix, transformation_matrix):
        """Build global forcing vector for concentrated axial load at distance 'a' from start"""
        # Equivalent nodal forces in local coordinates
        equivalent_forcing_element_coordinates = np.array([P * (1 - a / L), P * a / L])
        
        # Transform to global coordinates  
        big_F = np.transpose(transformation_matrix) @ equivalent_forcing_element_coordinates
        
        # Assemble into structural force vector
        CAL_global_forcing_vector = assembly_matrix @ big_F
        return CAL_global_forcing_vector, equivalent_forcing_element_coordinates, big_F

    def build_UDL_normal_global_forcing_vector(self, w_bar, L, assembly_matrix, transformation_matrix):
        """Build global forcing vector for UDL normal to bar axis (creates axial forces due to projection)"""
        # Extract angle from transformation matrix
        angle_rad = np.arctan2(transformation_matrix[0, 1], transformation_matrix[0, 0])
        
        # Project normal load onto bar axis
        w_axial = w_bar * np.sin(angle_rad)
        
        # Equivalent nodal forces in local coordinates
        equivalent_forcing_element_coordinates = np.array([w_axial * L / 2, w_axial * L / 2])
        
        # Transform to global coordinates
        big_F = np.transpose(transformation_matrix) @ equivalent_forcing_element_coordinates
        
        # Assemble into structural force vector
        UDL_normal_global_forcing_vector = assembly_matrix @ big_F
        return UDL_normal_global_forcing_vector, equivalent_forcing_element_coordinates, big_F

    # REACTION CALCULATION FUNCTIONS
    def calculate_reactions(self, config, results):
        """Calculate reaction forces at all constrained nodes for 2D bars"""
        reactions = {}
        nodes = config['nodes']
        elements = config['elements']
        
        # Initialize reaction storage for all nodes
        for node_id, node_info in nodes.items():
            if self.is_constrained_node(node_info['type']):
                reactions[node_id] = {
                    'type': node_info['type'],
                    'coords': node_info['coords'],
                    'forces': np.zeros(2),  # [Rx, Ry] for 2D bars
                    'constrained_dofs': self.get_constrained_dofs(node_info['type'])
                }
        
        # For each constrained node, calculate reactions
        for node_id in reactions.keys():
            element_forces = np.zeros(2)  # Forces from connected elements
            load_forces = np.zeros(2)     # Forces from loads on connected elements
            
            # Find all elements connected to this node
            for elem_idx, elem in enumerate(elements, 1):
                elem_key = f'element_{elem_idx}'
                
                # Check if this element connects to the current node
                if elem['start_node'] == node_id:
                    # Node is at start of element (indices 0-1)
                    element_forces += results[elem_key]['global_forces'][:2]
                    
                    # Add load contributions from this element
                    for load_info in results[elem_key]['loads']:
                        load_forces += load_info['equiv_global'][:2]
                        
                elif elem['end_node'] == node_id:
                    # Node is at end of element (indices 2-3)
                    element_forces += results[elem_key]['global_forces'][2:4]
                    
                    # Add load contributions from this element
                    for load_info in results[elem_key]['loads']:
                        load_forces += load_info['equiv_global'][2:4]
            
            # Calculate reaction: Element forces - Load equivalent forces
            reactions[node_id]['forces'] = element_forces - load_forces
            
            # Zero out reactions for unconstrained DOFs
            constrained_dofs = reactions[node_id]['constrained_dofs']
            for i in range(2):
                if not constrained_dofs[i]:
                    reactions[node_id]['forces'][i] = 0.0
        
        return reactions
    
    def is_constrained_node(self, node_type):
        """Check if a node has any constraints"""
        return node_type in ['fixed', 'pinned', 'roller_x', 'roller_y']
    
    def get_constrained_dofs(self, node_type):
        """Get which DOFs are constrained for each node type (2D bars only have translations)"""
        constraints = {
            'fixed': [True, True],      # ux, uy both constrained
            'pinned': [True, True],     # ux, uy both constrained (same as fixed for bars)
            'roller_x': [False, True],  # uy constrained
            'roller_y': [True, False],  # ux constrained  
            'free': [False, False]      # All free
        }
        return constraints.get(node_type, [False, False])

    def print_reactions(self, reactions):
        """Print reaction forces in a formatted table for 2D bars"""
        if not reactions:
            print("\nNo reaction forces to display (no constrained nodes)")
            return
        
        print(f"\n{'='*60}")
        print(f"{'REACTION FORCES AT CONSTRAINED NODES':^60}")
        print(f"{'='*60}")
        
        print(f"{'Node':<6} {'Type':<15} {'Rx [N]':<15} {'Ry [N]':<15}")
        print(f"{'-'*60}")
        
        for node_id, reaction_data in reactions.items():
            node_type = reaction_data['type']
            forces = reaction_data['forces']
            constrained = reaction_data['constrained_dofs']
            
            # Format values, show only constrained DOFs
            rx_str = f"{forces[0]:12.3f}" if constrained[0] else "---"
            ry_str = f"{forces[1]:12.3f}" if constrained[1] else "---"
            
            print(f"{node_id:<6} {node_type:<15} {rx_str:<15} {ry_str:<15}")
        
        print(f"{'-'*60}\n")

    def print_element_stresses(self, results, elements):
        """Print normal stresses for all elements in a formatted table"""
        print(f"\n{'='*80}")
        print(f"{'ELEMENT NORMAL STRESSES':^80}")
        print(f"{'='*80}")
        
        print(f"{'Element':<8} {'Area [m²]':<12} {'Axial Force [N]':<16} {'Normal Stress [Pa]':<18} {'Status':<12}")
        print(f"{'-'*80}")
        
        for i, elem in enumerate(elements, 1):
            elem_key = f'element_{i}'
            area = elem['area']
            axial_force = results[elem_key]['local_forces'][0]
            normal_stress = results[elem_key]['normal_stress']
            
            # Determine stress status
            if abs(normal_stress) < 1e-6:
                status = "No Load"
            elif normal_stress > 0:
                status = "Tension"
            else:
                status = "Compression"
            
            print(f"{i:<8} {area:11.3e} {axial_force:15.3e} {normal_stress:17.3e} {status:<12}")
        
        print(f"{'-'*80}\n")

    def add_reaction_arrows_to_plot(self, ax, reactions, scale_factor=1.0):
        """Add reaction force arrows to the plot for 2D bars"""
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
                dx = arrow_length if forces[0] > 0 else -arrow_length
                ax.arrow(x, y, dx, 0, head_width=0.1, head_length=0.05, 
                        fc='red', ec='red', linewidth=2, alpha=0.8)
                ax.text(x + dx/2, y + 0.15, f'Rx={forces[0]:.1f}N', 
                    ha='center', va='bottom', fontsize=8, color='red', weight='bold')
            
            if constrained[1] and abs(forces[1]) > 1e-6:  # Ry  
                dy = arrow_length if forces[1] > 0 else -arrow_length
                ax.arrow(x, y, 0, dy, head_width=0.1, head_length=0.05,
                        fc='blue', ec='blue', linewidth=2, alpha=0.8)
                ax.text(x - 0.15, y + dy/2, f'Ry={forces[1]:.1f}N',
                    ha='right', va='center', fontsize=8, color='blue', weight='bold',
                    rotation=90)

    # MATRIX PRINTING FUNCTIONS (simplified for bars)
    def find_common_factor(self, matrix):
        """Find common factor to simplify matrix display"""
        if matrix.size == 0:
            return 1
        
        non_zero_values = matrix[matrix != 0]
        if len(non_zero_values) == 0:
            return 1
        
        abs_values = np.abs(non_zero_values)
        significant_values = abs_values[abs_values > 1e-15]
        if len(significant_values) == 0:
            return 1
        
        median_val = np.median(significant_values)
        
        # Engineering powers of 10^3
        engineering_powers = range(-6, 10)
        
        best_factor = 1
        best_score = float('inf')
        
        for power in engineering_powers:
            factor = 10**(power * 3)
            
            simplified = significant_values / factor
            
            in_range = np.sum((simplified >= 0.1) & (simplified <= 999))
            total_values = len(simplified)
            in_ideal_range = np.sum((simplified >= 1) & (simplified <= 99))
            
            score = (total_values - in_range) * 10 + (total_values - in_ideal_range)
            
            if np.any(simplified > 9999) or np.any(simplified < 0.001):
                score += 1000
            
            if score < best_score:
                best_score = score
                best_factor = factor
        
        return best_factor

    def get_engineering_units(self, base_units, factor):
        """Convert base units to engineering units based on the factor"""
        unit_mappings = {
            "[N/m]": {1e-3: "[mN/m]", 1: "[N/m]", 1e3: "[kN/m]", 1e6: "[MN/m]"},
            "[N]": {1e-3: "[mN]", 1: "[N]", 1e3: "[kN]", 1e6: "[MN]"},
            "[m]": {1e-6: "[μm]", 1e-3: "[mm]", 1: "[m]", 1e3: "[km]"},
            "[Pa]": {1e-3: "[mPa]", 1: "[Pa]", 1e3: "[kPa]", 1e6: "[MPa]", 1e9: "[GPa]"},
            "[-]": {1: "[-]"}
        }
        
        if base_units in unit_mappings and factor in unit_mappings[base_units]:
            return unit_mappings[base_units][factor]
        else:
            return base_units

    def print_matrix_standard(self, name, matrix, units=""):
        factor = self.find_common_factor(matrix)
        engineering_units = self.get_engineering_units(units, factor)
        factor_text = f" (×{factor:.0e})" if factor != 1 else ""
        
        print(f"\n{name}{factor_text} {engineering_units}:")
        print("-" * (len(name) + len(factor_text) + len(engineering_units) + 1))
        
        simplified_matrix = matrix / factor
        
        if matrix.ndim == 1:
            for i, val in enumerate(simplified_matrix):
                print(f"  [{i:2d}]: {val:12.4f}")
        else:
            print(np.array_str(simplified_matrix, precision=4, suppress_small=True))
        print()

    def extract_node_coordinates_and_dofs(self, nodes_config, total_structural_deflections):
        """Extract node coordinates and DOFs for 2D bars (only translations)"""
        nodes = {}
        
        for node_id, node_info in nodes_config.items():
            node_coords = tuple(node_info['coords'])
            nodes[node_coords] = {'dofs': [], 'node_id': node_id}
        
        dof_names = ['ux', 'uy']
        dof_types = ['translation', 'translation']
        
        all_dofs = []
        for node_id, node_info in nodes_config.items():
            node_type = node_info['type']
            for local_dof in range(2):  # Only 2 DOFs per node for bars
                dof_name = f"Node{node_id}_{dof_names[local_dof]}"
                
                is_constrained = False
                if node_type in ['fixed', 'pinned']:
                    is_constrained = True
                elif node_type == 'roller_x' and local_dof == 1:  # uy constrained
                    is_constrained = True
                elif node_type == 'roller_y' and local_dof == 0:  # ux constrained
                    is_constrained = True
                
                all_dofs.append({
                    'node_id': node_id,
                    'local_dof': local_dof,
                    'name': dof_name,
                    'type': dof_types[local_dof],
                    'constrained': is_constrained
                })
        
        unconstrained_dof_idx = 0
        for dof_info in all_dofs:
            node_id = dof_info['node_id']
            node_coords = None
            for nid, ninfo in nodes_config.items():
                if nid == node_id:
                    node_coords = tuple(ninfo['coords'])
                    break
            
            if dof_info['constrained']:
                dof_value = 0.0
            else:
                if unconstrained_dof_idx < len(total_structural_deflections):
                    dof_value = total_structural_deflections[unconstrained_dof_idx]
                    unconstrained_dof_idx += 1
                else:
                    dof_value = 0.0
            
            unit = "[m]"  # Only translations for bars
            nodes[node_coords]['dofs'].append({
                'name': dof_info['name'],
                'value': dof_value,
                'unit': unit,
                'constrained': dof_info['constrained']
            })
        
        return nodes

    def add_node_dof_textboxes(self, ax, nodes, nodes_config):
        """Add text boxes showing DOF values at each node for 2D bars"""
        for node_coords, node_data in nodes.items():
            node_id = node_data['node_id']
            node_type = nodes_config[node_id]['type']
            
            if node_data['dofs']:
                if node_type in ['fixed', 'pinned']:
                    continue  # Skip fully constrained nodes
                
                display_dofs = node_data['dofs']
                
                if not display_dofs:
                    continue

                text_lines = [f"Node {node_id} ({node_type}):"]
                text_lines.append("─" * 20)
                
                for dof in display_dofs:
                    text_lines.append(f"{dof['name']}: {dof['value']:.3e} {dof['unit']}")
                
                text_content = "\n".join(text_lines)
                
                offset_x = 0.3
                offset_y = 0.2
                if node_coords[0] > 2.5:
                    offset_x = -1.0
                if node_coords[1] > 0:
                    offset_y = -0.3
                
                ax.text(node_coords[0] + offset_x, node_coords[1] + offset_y, 
                    text_content, fontsize=9, 
                    bbox=dict(boxstyle='round,pad=0.5', facecolor='lightblue', alpha=0.8),
                    verticalalignment='top', fontfamily='monospace')
                
                ax.plot(node_coords[0], node_coords[1], 'ko', markersize=8, zorder=10)
                ax.plot(node_coords[0], node_coords[1], 'wo', markersize=6, zorder=11)

    def plot_deflected_shape_with_total_deflections(self, ax, start_node_coords, end_node_coords, 
                                                  total_structural_deflections, assembly_matrix, 
                                                  disp_mag, element_name):
        """Plot deflected shape for 2D bar elements (straight lines only)"""
        # Extract element deflections
        element_deflections_global = np.transpose(assembly_matrix) @ total_structural_deflections
        
        # Node coordinates
        node1XG, node1YG = start_node_coords[0], start_node_coords[1]
        node2XG, node2YG = end_node_coords[0], end_node_coords[1]
        
        # Undeformed coordinates
        undeflected_XG = np.array([node1XG, node2XG])
        undeflected_YG = np.array([node1YG, node2YG])
        
        # Deflected coordinates (linear for bars)
        deflected_XG = undeflected_XG + disp_mag * np.array([element_deflections_global[0], element_deflections_global[2]])
        deflected_YG = undeflected_YG + disp_mag * np.array([element_deflections_global[1], element_deflections_global[3]])
        
        # Plot lines
        ax.plot(undeflected_XG, undeflected_YG, 'k-', linewidth=3, 
               label='Undeformed Shape' if element_name == "Element 1" else "")
        ax.plot(deflected_XG, deflected_YG, 'r-', linewidth=3, 
               label='Deflected Shape' if element_name == "Element 1" else "")
        
        # Add element label
        mid_x = (undeflected_XG[0] + undeflected_XG[1]) / 2
        mid_y = (undeflected_YG[0] + undeflected_YG[1]) / 2
        ax.text(mid_x, mid_y + 0.1, element_name, 
               fontsize=10, ha='center', weight='bold')

    def get_load_type_label(self, load_type):
        """Get the proper label for each load type"""
        load_labels = {
            'DAL': 'DAL (Distributed Axial Load)',
            'CAL': 'CAL (Concentrated Axial Load)',
            'UDL_normal': 'UDL Normal to Axis'
        }
        return load_labels.get(load_type, load_type)

    def run_analysis(self, config):
        """Main analysis function for 2D bar structures"""
        
        print("\n" + "="*100)
        print("                        2D BAR ANALYSIS RESULTS")
        print("="*100)
        
        # Extract configuration
        elements = config['elements']
        loads = config['loads']
        normal_forces = config['normal_forces']
        display_magnitude = config['display_magnitude']
        show_matrices = config.get('show_matrices', False)
        
        # Calculate element stiffness matrices
        results = {}
        total_global_forcing_vector = np.array(normal_forces, dtype=float)
        
        # Store load information for each element
        element_loads = {i: [] for i in range(1, len(elements) + 1)}
        
        for i, elem in enumerate(elements, 1):
            elem_key = f'element_{i}'
            results[elem_key] = {}
            
            # Local stiffness (2x2 for bars)
            local_stiffness = self.calc_local_coordinates_element_stiffness_matrix(
                elem['youngs_modulus'], elem['area'], elem['length'])
            results[elem_key]['local_stiffness'] = local_stiffness
            
            # Global stiffness and transformation (4x4 and 2x4)
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
        
        # Add loads
        all_load_vectors = []
        
        for load in loads:
            elem_idx = load['element'] - 1
            elem = elements[elem_idx]
            elem_key = f'element_{load["element"]}'
            
            load_info = {
                'type': load['type'],
                'magnitude': load['magnitude'],
                'element': load['element']
            }
            
            if load['type'] == 'DAL':
                load_vector, equiv_local, equiv_global = self.build_DAL_global_forcing_vector(
                    load['magnitude'], elem['length'], 
                    elem['assembly_matrix'], results[elem_key]['transformation_matrix'])
                load_info['description'] = f"Distributed Axial Load: {load['magnitude']} N/m"
                
            elif load['type'] == 'CAL':
                a = load.get('a', elem['length']/2)
                load_vector, equiv_local, equiv_global = self.build_CAL_global_forcing_vector(
                    load['magnitude'], elem['length'], a,
                    elem['assembly_matrix'], results[elem_key]['transformation_matrix'])
                load_info['description'] = f"Concentrated Axial Load: {load['magnitude']} N at x={a} m"
                load_info['a'] = a
                
            elif load['type'] == 'UDL_normal':
                load_vector, equiv_local, equiv_global = self.build_UDL_normal_global_forcing_vector(
                    load['magnitude'], elem['length'],
                    elem['assembly_matrix'], results[elem_key]['transformation_matrix'])
                load_info['description'] = f"UDL Normal to Axis: {load['magnitude']} N/m"

            load_info['equiv_local'] = equiv_local
            load_info['equiv_global'] = equiv_global
            load_info['global_vector'] = load_vector
            
            element_loads[load['element']].append(load_info)
            all_load_vectors.append(load_info)
            
            total_global_forcing_vector += load_vector
        
        # Store load information in results
        for i in range(1, len(elements) + 1):
            elem_key = f'element_{i}'
            results[elem_key]['loads'] = element_loads[i]
        
        # Solve for deflections
        total_structural_deflections = self.calc_structural_deflections(
            total_global_stiffness_matrix, total_global_forcing_vector)
        
        # Calculate element deflections and forces
        for i, elem in enumerate(elements, 1):
            elem_key = f'element_{i}'
            
            global_deflections = self.calc_global_element_structural_deflections(
                total_structural_deflections, elem['assembly_matrix'])
            results[elem_key]['global_deflections'] = global_deflections
            
            local_deflections = self.calc_local_element_structural_deflections(
                global_deflections, results[elem_key]['transformation_matrix'])
            results[elem_key]['local_deflections'] = local_deflections
            
            local_forces = self.calc_local_element_nodal_forces(
                results[elem_key]['local_stiffness'], local_deflections)
            global_forces = self.calc_global_element_nodal_forces(
                results[elem_key]['global_stiffness'], global_deflections)
            results[elem_key]['local_forces'] = local_forces
            results[elem_key]['global_forces'] = global_forces
            
            # Calculate normal stress for this element
            normal_stress = self.calc_element_normal_stress(local_forces, elem['area'])
            results[elem_key]['normal_stress'] = normal_stress
        
        # Calculate reaction forces
        reactions = self.calculate_reactions(config, results)
        
        # Print results
        for i, elem in enumerate(elements, 1):
            elem_key = f'element_{i}'
            print(f"\n{'-'*60}")
            print(f"ELEMENT {i}")
            print(f"{'-'*60}")
            
            self.print_matrix_standard(f"LOCAL stiffness matrix K_Local_{i}", results[elem_key]['local_stiffness'], "[N/m]")
            self.print_matrix_standard(f"GLOBAL stiffness matrix K_Global_{i}", results[elem_key]['global_stiffness'], "[N/m]")
            self.print_matrix_standard(f"Transformation matrix T_{i}", results[elem_key]['transformation_matrix'], "[-]")
            self.print_matrix_standard(f"GLOBAL element deflections d_Global_{i}", results[elem_key]['global_deflections'], "[m]")
            self.print_matrix_standard(f"LOCAL element deflections d_Local_{i}", results[elem_key]['local_deflections'], "[m]")
            self.print_matrix_standard(f"LOCAL element forces F_Local_{i}", results[elem_key]['local_forces'], "[N]")
            self.print_matrix_standard(f"GLOBAL element forces F_Global_{i}", results[elem_key]['global_forces'], "[N]")
            self.print_matrix_standard(f"Normal stress σ_{i}", np.array([results[elem_key]['normal_stress']]), "[Pa]")

            # Print load equivalent forcing vectors
            if results[elem_key]['loads']:
                print(f"\n{'LOAD FORCING VECTORS FOR ELEMENT ' + str(i)}")
                print(f"{'-'*60}")
                for j, load_info in enumerate(results[elem_key]['loads'], 1):
                    load_label = self.get_load_type_label(load_info['type'])
                    print(f"\n{load_label} {j}: {load_info['description']}")
                    self.print_matrix_standard(f"Equivalent LOCAL forces f_Local_{i}_{j}", 
                              load_info['equiv_local'], "[N]")
                    self.print_matrix_standard(f"Equivalent GLOBAL forces F_Global_{i}_{j}", 
                              load_info['equiv_global'], "[N]")
                    self.print_matrix_standard(f"GLOBAL forcing vector Q_Global_{i}_{j}", 
                              load_info['global_vector'], "[N]")
        
        # Structure level results
        print(f"\n{'-'*60}")
        print("STRUCTURE LEVEL")
        print(f"{'-'*60}")
        
        for i, elem in enumerate(elements, 1):
            elem_key = f'element_{i}'
            self.print_matrix_standard(f"Element {i} stiffness matrix contribution [K]_{i}", 
                      results[elem_key]['global_stiffness_matrix'], "[N/m]")
        
        self.print_matrix_standard("TOTAL global stiffness matrix [K]_total", total_global_stiffness_matrix, "[N/m]")
        
        # Print all load contributions at structure level
        if all_load_vectors:
            print(f"\n{'STRUCTURE LEVEL LOAD CONTRIBUTIONS':^60}")
            print(f"{'-'*60}")
            for j, load_info in enumerate(all_load_vectors, 1):
                load_label = self.get_load_type_label(load_info['type'])
                print(f"\n{load_label} Structure Load {j} (Element {load_info['element']}): {load_info['description']}")
                self.print_matrix_standard(f"GLOBAL forcing vector contribution Q_{load_label}_{j}", 
                          load_info['global_vector'], "[N]")
        
        self.print_matrix_standard("Total structural deflections d_total", total_structural_deflections, "[m]")
        self.print_matrix_standard("TOTAL global forcing vector Q_total", total_global_forcing_vector, "[N]")
        
        # Print reaction forces
        self.print_reactions(reactions)
        
        # Print element stresses
        self.print_element_stresses(results, elements)
        
        print("\n" + "="*100)
        print("                                END OF RESULTS")
        print("="*100 + "\n")
        
        # Plot results
        fig, ax = plt.subplots(figsize=(12, 8))
        
        # Plot each element's deflected shape
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
                display_magnitude, f"Element {i}")
            
        # Extract node deflections
        nodes_with_dofs = self.extract_node_coordinates_and_dofs(
            config['nodes'], 
            total_structural_deflections
        )

        # Add DOF textboxes to plot
        self.add_node_dof_textboxes(ax, nodes_with_dofs, config['nodes'])

        # Add reaction arrows to plot
        self.add_reaction_arrows_to_plot(ax, reactions)

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
                color = 'cyan'
            else:
                marker = 'x'
                color = 'black'
            
            ax.plot(x, y, marker=marker, color=color, markersize=12, label=f"Node {node_id} ({n_type})")

        # Remove duplicate legend entries
        handles, labels = ax.get_legend_handles_labels()
        by_label = dict(zip(labels, handles))
        ax.legend(by_label.values(), by_label.keys())
        
        ax.grid(True, alpha=0.3)
        ax.set_xlabel('X (m)')
        ax.set_ylabel('Y (m)')
        ax.set_title(f'2D Bar Analysis - Deflected Shape with Node DOFs and Reactions\n(Magnification Factor: {display_magnitude}x)')
        plt.tight_layout()
        plt.show()
        
        return results