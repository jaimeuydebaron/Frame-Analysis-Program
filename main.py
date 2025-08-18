import numpy as np
from frame_analysis import FrameAnalysis

def main():
    """
    Example usage of the Frame Analysis tool.
    Configure your problem here and run the analysis.
    """
    analyzer = FrameAnalysis()
    
    # CONFIGURATION - EDIT THIS SECTION FOR YOUR PROBLEM
    config = {
        # Node definitions: [x, y, 'support_type', 'label'] or [x, y, 'support_type']
        'nodes': [
            [0, 0, 'fixed', 'Node 1'],      # Node 1: Fixed support at origin
            [0, 2.5, 'free', 'Node 2'],     # Node 2: Free node
            [2, 1, 'free', 'Node 3'],       # Node 3: Free node
            [2, 0, 'fixed', 'Node 4']         # Node 4: Pin support
        ],
        
        # Element properties
        'elements': [
            {  # Element 1
                'node_indices': [0, 1],  # Connect Node 1 to Node 2 (0-indexed)
                'length': 2.5,                    # m
                'second_moment_of_area': 7.5e-6,  # m^4
                'area': 9e-4,                   # m^2
                'youngs_modulus': 200e9,        # Pa
                'angle': 90,                     # degrees
                'assembly_matrix': np.array([[0, 0, 0, 1, 0, 0],
                                             [0, 0, 0, 0, 1, 0],
                                             [0, 0, 0, 0, 0, 1],
                                             [0, 0, 0, 0, 0, 0],
                                             [0, 0, 0, 0, 0, 0],
                                             [0, 0, 0, 0, 0, 0]])
            },
            {  # Element 2
                'node_indices': [1, 2],  # Connect Node 2 to Node 3 (0-indexed)
                'length': 2.5,                    # m
                'second_moment_of_area': 7.5e-6,  # m^4
                'area': 9e-4,                   # m^2
                'youngs_modulus': 200e9,        # Pa
                'angle': 360 - 36.87,                     # degrees
                'assembly_matrix': np.array([[1, 0, 0, 0, 0, 0],
                                             [0, 1, 0, 0, 0, 0],   
                                             [0, 0, 1, 0, 0, 0],
                                             [0, 0, 0, 1, 0, 0],
                                             [0, 0, 0, 0, 1, 0],
                                             [0, 0, 0, 0, 0, 1]])
            },
            {  # Element 3
                'node_indices': [2, 3],  # Connect Node 3 to Node 4 (0-indexed)
                'length': 1,                    # m
                'second_moment_of_area': 7.5e-6,  # m^4
                'area': 9e-4,                   # m^2
                'youngs_modulus': 200e9,        # Pa
                'angle': 270,                     # degrees
                'assembly_matrix': np.array([[0, 0, 0, 0, 0, 0],
                                             [0, 0, 0, 0, 0, 0],   
                                             [0, 0, 0, 0, 0, 0],
                                             [1, 0, 0, 0, 0, 0],
                                             [0, 1, 0, 0, 0, 0],
                                             [0, 0, 1, 0, 0, 0]])
            }
        ],
        
        # DOF definitions - specify what each DOF represents and its type
        'dof_definitions': [
            {'name': 'Node2_ux', 'type': 'translation', 'node_index': 1},     # DOF 0: Node 2 x-displacement
            {'name': 'Node2_uy', 'type': 'translation', 'node_index': 1},     # DOF 1: Node 2 y-displacement
            {'name': 'Node2_rz', 'type': 'rotation', 'node_index': 1},        # DOF 2: Node 2 z-rotation  
            {'name': 'Node3_ux', 'type': 'translation', 'node_index': 2},     # DOF 3: Node 3 x-displacement
            {'name': 'Node3_uy', 'type': 'translation', 'node_index': 2},     # DOF 4: Node 3 y-displacement
            {'name': 'Node3_rz', 'type': 'rotation', 'node_index': 2}         # DOF 5: Node 3 z-rotation
        ],
        
        # Distributed loads (UDL, LVL_peakXL, LVL_peakXZERO, PDL)
        'loads': [
            {'type': 'LVL_peakXZERO', 'element': 1, 'magnitude': -73575},  # N/m on element 1
            # Add more loads as needed:
            # {'type': 'PDL', 'element': 2, 'magnitude': 20000},  # N on element 2
            # {'type': 'LVL_peakXL', 'element': 1, 'magnitude': 10000}  # N/m peak at x=L
        ],
        
        # Normal nodal forces - adjust size based on your DOFs
        'normal_forces': [0, 0, 0, 0, 0, 0],  # [N, N⋅m] corresponding to DOFs
        
        # Display options
        'display_magnitude': 100,        # Magnification factor for deflection plot
        'show_matrices': True,  # True for vertical matrix display, False for standard
        'show_deflection_labels': True,   # True to show deflection values on plot
        'validate_node_types': True      # True to validate DOF assignments against node types
    }
    
    # Run the analysis
    results = analyzer.run_analysis(config)
    
    # Optional: Access specific results
    # print("Element 1 local deflections:", results['element_1']['local_deflections'])
    # print("Element 2 global forces:", results['element_2']['global_forces'])
    
    return results

if __name__ == "__main__":
    results = main()