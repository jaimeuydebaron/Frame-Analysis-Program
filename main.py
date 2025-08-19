import numpy as np
from frame_analysis import FrameAnalysis

def main():
    analyzer = FrameAnalysis()

    # Material and section properties - replace with your values
    E = 200e9  # Young's modulus [Pa] - typical steel: 200 GPa
    I = 6e-6   # Second moment of area [m^4]
    A = 2e-4   # Cross-sectional area [m^2]
    
    config = {
        # NODES: Define all nodes in your structure
        # Each node needs: coordinates [x, y] in meters and boundary condition type
        'nodes': {
            # Node ID: {'coords': np.array([x_coordinate, y_coordinate]), 'type': 'boundary_condition_type'}
            1: {'coords': np.array([0, 0]),   'type': 'pinned'},                    
            2: {'coords': np.array([3, 4]),   'type': 'free'},                      
            3: {'coords': np.array([6, 4]),   'type': 'symmetrical_roller_y'},     
            4: {'coords': np.array([6, 0]),   'type': 'symmetrical_roller_y'},    
        },

        # ELEMENTS: Define all frame elements connecting your nodes
        'elements': [
            {  # Element connecting specified start_node to end_node
                'start_node': 2, 'end_node': 3,          
                'length': 3,            # Element length [m]
                'angle': 0,             # Element angle [degrees]
                'assembly_matrix': np.array([[0, 0, 0, 0, 0, 0],  
                                             [1, 0, 0, 0, 0, 0],  
                                             [0, 1, 0, 0, 0, 0],
                                             [0, 0, 1, 0, 0, 0], 
                                             [0, 0, 0, 0, 1, 0], 
                                             [0, 0, 0, 0, 0, 0]]),

                'second_moment_of_area': I, 'area': A, 'youngs_modulus': E,
            },
            {  # Second element
                'start_node': 1, 'end_node': 2,         
                'length': 5,            # Element length [m]
                'angle': 53.13,         # Element angle [degrees]
                'assembly_matrix': np.array([[0, 0, 1, 0, 0, 0],
                                             [0, 0, 0, 1, 0, 0],   
                                             [0, 0, 0, 0, 1, 0],
                                             [0, 0, 0, 0, 0, 1],
                                             [0, 0, 0, 0, 0, 0],
                                             [0, 0, 0, 0, 0, 0]]),

                'second_moment_of_area': I, 'area': A, 'youngs_modulus': E,
            },
            {  # Third element
                'start_node': 2, 'end_node': 4,       
                'length': 5,            # Element length [m]
                'angle': 306.87,        # Element angle [degrees]
                'assembly_matrix': np.array([[0, 0, 0, 0, 0, 0],
                                             [1, 0, 0, 0, 0, 0],   
                                             [0, 1, 0, 0, 0, 0],
                                             [0, 0, 1, 0, 0, 0],
                                             [0, 0, 0, 0, 0, 0],
                                             [0, 0, 0, 0, 1, 0]]),

                'second_moment_of_area': I, 'area': A, 'youngs_modulus': E,
            },
            {  # Fourth element
                'start_node': 1, 'end_node': 4,   
                'length': 6,            # Element length [m]
                'angle': 0,             # Element angle [degrees]
                'assembly_matrix': np.array([[0, 0, 1, 0, 0, 0],
                                             [0, 0, 0, 0, 0, 0],   
                                             [0, 0, 0, 0, 0, 0],
                                             [0, 0, 0, 0, 0, 0],
                                             [0, 0, 0, 0, 0, 0],
                                             [0, 0, 0, 0, 1, 0]]),

                'second_moment_of_area': I, 'area': A, 'youngs_modulus': E,         
            }
        ],

        # LOADS: Define distributed and point loads on elements
        'loads': [
            # Distributed loads - applied to specific elements
            {'type': 'UDL', 'element': 4, 'magnitude': -40000, 'direction': 'local'},  # Uniformly distributed load
            
            # Other load types (uncomment and modify as needed):
            # {'type': 'UDL', 'element': 2, 'magnitude': -5000, 'direction': 'vertical'},     # UDL in vertical direction
            # {'type': 'LVL_peakXZERO', 'element': 1, 'magnitude': 10000},                    # Linearly varying load (peak at x=0)
            # {'type': 'LVL_peakXL', 'element': 1, 'magnitude': 10000},                       # Linearly varying load (peak at x=L)
            # {'type': 'PDL', 'element': 2, 'magnitude': 20000, 'direction': 'vertical', 'a': 1.5},  # Point load at distance 'a'
            # {'type': 'DAL', 'element': 2, 'magnitude': 15000},                              # Distributed axial load
            # {'type': 'CAL', 'element': 2, 'magnitude': 10000, 'a': 1.5},                   # Concentrated axial load at distance 'a'
        ],
        
        # NODAL FORCES: Direct forces/moments applied at nodes (in global DOF order)
        # Length must match total number of unconstrained DOFs in your system
        # Format: [Force_DOF1, Force_DOF2, ..., Moment_DOF_n] in [N] and [N⋅m]
        'normal_forces': [0, 0, 0, 0, 0, 0],  # Zero forces/moments on all unconstrained DOFs
        
        # DISPLAY OPTIONS: Control output format and visualization
        'display_magnitude': 10,        # Deflection magnification factor for plotting
        'show_matrices': True,          # True = vertical matrix display, False = compact display
        'show_deflection_labels': True  # True = show deflection values on plot
    }

    # BOUNDARY CONDITION TYPES REFERENCE:
    # 'fixed':                All DOFs constrained (ux=0, uy=0, rz=0)
    # 'pinned':               Translations constrained, rotation free (ux=0, uy=0, rz=free)
    # 'roller_x':             Y-translation constrained (ux=free, uy=0, rz=free)
    # 'roller_y':             X-translation constrained (ux=0, uy=free, rz=free)
    # 'symmetrical_roller_x': Y-translation & rotation constrained (ux=free, uy=0, rz=0)
    # 'symmetrical_roller_y': X-translation & rotation constrained (ux=0, uy=free, rz=0)
    # 'free':                 All DOFs unconstrained (ux=free, uy=free, rz=free)

    # LOAD TYPES REFERENCE:
    # 'UDL':           Uniformly distributed load [N/m]
    #                  - 'direction': 'local' (perpendicular to element) or 'vertical' (gravity direction)
    # 'LVL_peakXZERO': Linearly varying load, peak at start of element [N/m]
    # 'LVL_peakXL':    Linearly varying load, peak at end of element [N/m]  
    # 'PDL':           Point load [N]
    #                  - 'direction': 'local' or 'vertical'
    #                  - 'a': distance from start node [m] (optional, defaults to midspan)
    # 'DAL':           Distributed axial load [N/m]
    # 'CAL':           Concentrated axial load [N]
    #                  - 'a': distance from start node [m] (optional, defaults to midspan)

    # ASSEMBLY MATRIX GUIDE:
    # The assembly matrix maps element DOFs to global system DOFs
    # - Rows represent global DOFs (numbered by unconstrained DOFs in your system)
    # - Columns represent element DOFs: [ux_start, uy_start, rz_start, ux_end, uy_end, rz_end]
    # - Put '1' where element DOF corresponds to that global DOF, '0' elsewhere
    # - Total rows = total unconstrained DOFs across all nodes
    # - Must account for boundary conditions when numbering global DOFs

    results = analyzer.run_analysis(config)

if __name__ == "__main__":
    main()