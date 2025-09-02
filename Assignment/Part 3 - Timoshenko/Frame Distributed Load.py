import numpy as np
from frame_analysis import FrameAnalysis

def main():
    analyzer = FrameAnalysis()

    # Material and section properties
    E = 200e9  # Young's modulus [Pa] - typical steel: 200 GPa
    G = 77e9
    I = 1.36288e-5   # Second moment of area [m^4]
    A = 0.004496   # Cross-sectional area [m^2]
    c = 0.07
    A_s = 0.14 * 0.012
    
    
    config = {
        # NODES: Define all nodes in your structure
        'nodes': {
            1: {'coords': np.array([0, 0]),     'type': 'fixed'},        # Base left - pinned
            2: {'coords': np.array([1, 0]),     'type': 'fixed'},        # Base right - pinned
            3: {'coords': np.array([0, 2.5]),   'type': 'free'},          # Level 1 left
            4: {'coords': np.array([1, 2.5]),   'type': 'free'},          # Level 1 right
            5: {'coords': np.array([0, 5]),     'type': 'free'},          # Level 2 left
            6: {'coords': np.array([1, 5]),     'type': 'free'},          # Level 2 right
            7: {'coords': np.array([0, 7.5]),   'type': 'free'},          # Top
        },

        # ELEMENTS: Define all frame elements connecting your nodes
        'elements': [
            {   # Element 1: Node 1 → 3 (vertical column)
                'start_node': 1, 'end_node': 3,
                'length': 2.5,
                'angle': 90,  # Vertical up
                'assembly_matrix': np.array([
                    [0, 0, 0, 1, 0, 0],
                    [0, 0, 0, 0, 1, 0],
                    [0, 0, 0, 0, 0, 1],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                ]),
                'second_moment_of_area': I, 'area': A, 'youngs_modulus': E, 'shear_modulus': G, 'c': c, 'A_s': A_s,
            },
            {   # Element 2: Node 2 → 3 (diagonal)
                'start_node': 2, 'end_node': 3,
                'length': np.sqrt(2.5**2 + 1**2),
                'angle': 111.8,
                'assembly_matrix': np.array([
                    [0, 0, 0, 1, 0, 0],
                    [0, 0, 0, 0, 1, 0],
                    [0, 0, 0, 0, 0, 1],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                ]),
                'second_moment_of_area': I, 'area': A, 'youngs_modulus': E, 'shear_modulus': G, 'c': c, 'A_s': A_s,
            },
            {   # Element 3: Node 2 → 4 (vertical column)
                'start_node': 2, 'end_node': 4,
                'length': 2.5,
                'angle': 90,
                'assembly_matrix': np.array([
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 1, 0, 0],
                    [0, 0, 0, 0, 1, 0],
                    [0, 0, 0, 0, 0, 1],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                ]),
                'second_moment_of_area': I, 'area': A, 'youngs_modulus': E, 'shear_modulus': G, 'c': c, 'A_s': A_s,

            },
            {   # Element 4: Node 3 → 4 (horizontal beam)
                'start_node': 3, 'end_node': 4,
                'length': 1,
                'angle': 0,
                'assembly_matrix': np.array([
                    [1, 0, 0, 0, 0, 0],
                    [0, 1, 0, 0, 0, 0],
                    [0, 0, 1, 0, 0, 0],
                    [0, 0, 0, 1, 0, 0],
                    [0, 0, 0, 0, 1, 0],
                    [0, 0, 0, 0, 0, 1],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                ]),
                'second_moment_of_area': I, 'area': A, 'youngs_modulus': E, 'shear_modulus': G, 'c': c, 'A_s': A_s,
            },
            {   # Element 5: Node 3 → 5 (vertical column)
                'start_node': 3, 'end_node': 5,
                'length': 2.5,
                'angle': 90,
                'assembly_matrix': np.array([
                    [1, 0, 0, 0, 0, 0],
                    [0, 1, 0, 0, 0, 0],
                    [0, 0, 1, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 1, 0, 0],
                    [0, 0, 0, 0, 1, 0],
                    [0, 0, 0, 0, 0, 1],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                ]),
                'second_moment_of_area': I, 'area': A, 'youngs_modulus': E, 'shear_modulus': G, 'c': c, 'A_s': A_s,
            },
            {   # Element 6: Node 3 → 6 (diagonal)
                'start_node': 3, 'end_node': 6,
                'length': np.sqrt(2.5**2 + 1**2),
                'angle': 68.2,
                'assembly_matrix': np.array([
                    [1, 0, 0, 0, 0, 0],
                    [0, 1, 0, 0, 0, 0],
                    [0, 0, 1, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 1, 0, 0],
                    [0, 0, 0, 0, 1, 0],
                    [0, 0, 0, 0, 0, 1],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                ]),
                'second_moment_of_area': I, 'area': A, 'youngs_modulus': E, 'shear_modulus': G, 'c': c, 'A_s': A_s,
            },
            {   # Element 7: Node 4 → 6 (vertical column)
                'start_node': 4, 'end_node': 6,
                'length': 2.5,
                'angle': 90,
                'assembly_matrix': np.array([
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [1, 0, 0, 0, 0, 0],
                    [0, 1, 0, 0, 0, 0],
                    [0, 0, 1, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 1, 0, 0],
                    [0, 0, 0, 0, 1, 0],
                    [0, 0, 0, 0, 0, 1],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                ]),
                'second_moment_of_area': I, 'area': A, 'youngs_modulus': E, 'shear_modulus': G, 'c': c, 'A_s': A_s,
            },
            {   # Element 8: Node 5 → 6 (horizontal beam)
                'start_node': 5, 'end_node': 6,
                'length': 1,
                'angle': 0,
                'assembly_matrix': np.array([
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [1, 0, 0, 0, 0, 0],
                    [0, 1, 0, 0, 0, 0],
                    [0, 0, 1, 0, 0, 0],
                    [0, 0, 0, 1, 0, 0],
                    [0, 0, 0, 0, 1, 0],
                    [0, 0, 0, 0, 0, 1],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                ]),
                'second_moment_of_area': I, 'area': A, 'youngs_modulus': E, 'shear_modulus': G, 'c': c, 'A_s': A_s,
            },
            {   # Element 9: Node 5 → 7 (vertical column)
                'start_node': 5, 'end_node': 7,
                'length': 2.5,
                'angle': 90,
                'assembly_matrix': np.array([
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [1, 0, 0, 0, 0, 0],
                    [0, 1, 0, 0, 0, 0],
                    [0, 0, 1, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 1, 0, 0],
                    [0, 0, 0, 0, 1, 0],
                    [0, 0, 0, 0, 0, 1],
                ]),
                'second_moment_of_area': I, 'area': A, 'youngs_modulus': E, 'shear_modulus': G, 'c': c, 'A_s': A_s,
            },
            {   # Element 10: Node 6 → 7 (diagonal)
                'start_node': 6, 'end_node': 7,
                'length': np.sqrt(2.5**2 + 1**2),
                'angle': 111.8,
                'assembly_matrix': np.array([
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [0, 0, 0, 0, 0, 0],
                    [1, 0, 0, 0, 0, 0],
                    [0, 1, 0, 0, 0, 0],
                    [0, 0, 1, 0, 0, 0],
                    [0, 0, 0, 1, 0, 0],
                    [0, 0, 0, 0, 1, 0],
                    [0, 0, 0, 0, 0, 1],
                ]),
                'second_moment_of_area': I, 'area': A, 'youngs_modulus': E, 'shear_modulus': G, 'c': c, 'A_s': A_s,
            }
        ],

        # LOADS: Define distributed and point loads on elements
        'loads': [
            # Example loads - uncomment and modify as needed:

            {'type': 'UDL', 'element': 5, 'magnitude': 103948.2, 'direction': 'vertical'},  # 5 kN/m downward on beam 3-4
            {'type': 'UDL', 'element': 9, 'magnitude': 103948.2, 'direction': 'vertical'},  # 5 kN/m downward on beam 3-4

            # {'type': 'PDL', 'element': 8, 'magnitude': -10000, 'direction': 'vertical', 'a': 0.5},  # 10 kN downward at midspan of beam 5-6
        ],
        
        # NODAL FORCES: Direct forces/moments applied at nodes (in global DOF order)
        # Order: [3ux, 3uy, 3rz, 4ux, 4uy, 4rz, 5ux, 5uy, 5rz, 6ux, 6uy, 6rz, 7ux, 7uy, 7rz]
        'normal_forces': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        
        # DISPLAY OPTIONS: Control output format and visualization
        'display_magnitude': 100,        # Deflection magnification factor for plotting
        'show_matrices': False,          # True = vertical matrix display, False = compact display
        'show_deflection_labels': True,   # True = show deflection values on plot
        'show_plot': True,               # True = display plot, False = skip plot

        # NEW: Design parameters for maximum UDL calculation
        'yield_stress': 350e6,     # Material yield stress [Pa]
        'safety_factor': 2.5,      # Factor of safety for design
        'udl_elements': [5, 9],    # Elements to apply UDL for max load analysis
        'wind_depth': 8.0,         # Wind load depth [m]
        'wind_pressure_coeff': 0.6 # Wind pressure coefficient [Pa⋅s²/m²]
    }

    results = analyzer.run_analysis(config)

    # Separate optimization analysis
    # max_udl, wind_speed, critical_location = analyzer.find_maximum_udl(config)

if __name__ == "__main__":
    main()