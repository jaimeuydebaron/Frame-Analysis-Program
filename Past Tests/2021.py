import numpy as np
from frame_analysis import FrameAnalysis

def main():
    """
    Example usage of the Frame Analysis tool.
    Configure your problem here and run the analysis.
    """
    analyzer = FrameAnalysis()

    E = 200e9
    I = 1e-6
    A = 1e-4
    
    config = {
        # Define nodes with coordinates and boundary conditions
        'nodes': {

            1: {'coords': np.array([0, 0]),   'type': 'fixed'},
            2: {'coords': np.array([3, 0]),   'type': 'roller_x'},
            3: {'coords': np.array([5, 0]),   'type': 'free'},
        },

        # Element properties - reference nodes by ID
        'elements': [
            {  # Element 1: connects node 1 to node 2
                'start_node': 1,
                'end_node': 2,
                'length': 3,                    # m
                'angle': 0,                     # degrees

                'assembly_matrix': np.array([[0, 0, 0, 1, 0, 0],
                                             [0, 0, 0, 0, 0, 1],
                                             [0, 0, 0, 0, 0, 0],
                                             [0, 0, 0, 0, 0, 0],
                                             [0, 0, 0, 0, 0, 0]]),

                'second_moment_of_area': I,  # m^4
                'area': A,                   # m^2
                'youngs_modulus': E,        # Pa
            },
            {  # Element 2: connects node 1 to node 3
                'start_node': 2,
                'end_node': 3,
                'length': 2,                    # m
                'angle': 0,                 # degrees
                'assembly_matrix': np.array([[1, 0, 0, 0, 0, 0],
                                             [0, 0, 1, 0, 0, 0],   
                                             [0, 0, 0, 1, 0, 0],
                                             [0, 0, 0, 0, 1, 0],
                                             [0, 0, 0, 0, 0, 1]]),

                'second_moment_of_area': I,  # m^4
                'area': A,                   # m^2
                'youngs_modulus': E,        # Pa

            },
            # {  # Element 3: connects node 2 to node 3
            #     'start_node': 2,
            #     'end_node': 3,
            #     'length': 3,                    # m
            #     'angle': 90,                    # degrees
            #     'assembly_matrix': np.array([[1, 0, 0, 0, 0, 0],
            #                                  [0, 0, 1, 0, 0, 0],   
            #                                  [0, 0, 0, 1, 0, 0],
            #                                  [0, 0, 0, 0, 1, 0],
            #                                  [0, 0, 0, 0, 0, 1]]),

                # 'second_moment_of_area': I,  # m^4
                # 'area': A,                   # m^2
                # 'youngs_modulus': E,        # Pa
            # }
        ],

        # Distributed loads (UDL, LVL_peakXL, LVL_peakXZERO, PDL)
        'loads': [
            {'type': 'UDL', 'element': 1, 'magnitude': -40000, 'direction': 'vertical'},  # N/m on element 1
            # {'type': 'LVL_peakXZERO', 'element': 1, 'magnitude': 10000},  # N/m peak at x=0 on element 1
            # {'type': 'LVL_peakXL', 'element': 1, 'magnitude': 10000}  # N/m peak at x=L
            # {'type': 'PDL', 'element': 2, 'magnitude': 20000, 'direction': 'vertical', 'a': 1.5},  # N on element 2
            # {'type': 'DAL', 'element': 2, 'magnitude': 15000},  # N on element 2
            # {'type': 'CAL', 'element': 2, 'magnitude': 10000, 'a': 1.5},  # N on element 2
        ],
        
        # Normal nodal forces - adjust size based on your DOFs
        'normal_forces': [0, 0, 25000, -10000, 0],  # [N, N⋅m] corresponding to DOFs
        
        # Display options
        'display_magnitude': 10,        # Magnification factor for deflection plot
        'show_matrices': True,           # True for vertical matrix display, False for standard
        'show_deflection_labels': True   # True to show deflection values on plot
    }

    results = analyzer.run_analysis(config)

if __name__ == "__main__":
    main()