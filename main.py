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
        # Element properties
        'elements': [
            {  # Element 1
                'length': 3,                    # m
                'second_moment_of_area': 1e-5,  # m^4
                'area': 5e-4,                   # m^2
                'youngs_modulus': 200e9,        # Pa
                'angle': 90,                     # degrees
                'node_positions': [np.array([0, 0]), np.array([0, 3])],  # Global coordinates
                'assembly_matrix': np.array([[0, 0, 0, 1, 0, 0],
                                           [0, 0, 0, 0, 0, 1],
                                           [0, 0, 0, 0, 0, 0],
                                           [0, 0, 0, 0, 0, 0]])
            },
            {  # Element 2
                 'length': 4.5,                    # m
                'second_moment_of_area': 1e-5,  # m^4
                'area': 5e-4,                   # m^2
                'youngs_modulus': 200e9,        # Pa
                'angle': 0,                     # degrees
                'node_positions': [np.array([0, 3]), np.array([4.5, 3])],  # Global coordinates
                'assembly_matrix': np.array([[1, 0, 0, 0, 0, 0],
                                           [0, 0, 1, 0, 0, 0],   
                                           [0, 0, 0, 1, 0, 0],
                                           [0, 0, 0, 0, 0, 1]])
            },
            {  # Element 3
                'length': 3,                    # m
                'second_moment_of_area': 1e-5,  # m^4
                'area': 5e-4,                   # m^2
                'youngs_modulus': 200e9,        # Pa
                'angle': 270,                     # degrees
                'node_positions': [np.array([4.5, 3]), np.array([4.5, 0])],  # Global coordinates
                'assembly_matrix': np.array([[0, 0, 0, 0, 0, 0],
                                           [0, 0, 0, 0, 0, 0],   
                                           [1, 0, 0, 0, 0, 0],
                                           [0, 0, 1, 0, 0, 0]])
            }
        ],
        
        # DOF definitions - specify what each DOF represents and its type
        'dof_definitions': [
            {'name': 'Node1_rx', 'type': 'translation'},     # DOF 0: Node 2 x-displacement
            {'name': 'Node1_rz', 'type': 'rotation'},        # DOF 1: Node 1 z-rotation  
            {'name': 'Node2_rx', 'type': 'translation'},     # DOF 2: Node 2 y-displacement
            # {'name': 'Node2_ry', 'type': 'translation'},        # DOF 3: Node 2 z-rotation
            {'name': 'Node2_rz', 'type': 'rotation'}         # DOF 4: Node 3 z-rotation
            
        ],
        
        # Distributed loads (UDL, LVL_peakXL, LVL_peakXZERO, PDL)
        'loads': [
            # {'type': 'UDL', 'element': 1, 'magnitude': -40000},  # N/m on element 1
            # Add more loads as needed:
            # {'type': 'PDL', 'element': 2, 'magnitude': 20000},  # N on element 2
            # {'type': 'LVL_peakXL', 'element': 1, 'magnitude': 10000}  # N/m peak at x=L
        ],
        
        # Normal nodal forces - adjust size based on your DOFs
        'normal_forces': [10000, 0, 10000, 0],  # [N, N⋅m] corresponding to DOFs
        
        # Display options
        'display_magnitude': 10,        # Magnification factor for deflection plot
        'show_matrices': True,  # True for vertical matrix display, False for standard
        'show_deflection_labels': True   # True to show deflection values on plot
    }
    
    # Run the analysis
    results = analyzer.run_analysis(config)
    
    # Optional: Access specific results
    # print("Element 1 local deflections:", results['element_1']['local_deflections'])
    # print("Element 2 global forces:", results['element_2']['global_forces'])
    
    return results

if __name__ == "__main__":
    results = main()