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
                    # Global DOF mapping: [3ux, 3uy, 3rz, 4ux, 4uy, 4rz, 5ux, 5uy, 5rz, 6ux, 6uy, 6rz, 7ux, 7uy, 7rz]
                    # Element DOFs:       [1ux, 1uy, 1rz, 3ux, 3uy, 3rz]
                    [0, 0, 0, 1, 0, 0],  # 3ux ← 3ux (element DOF 4)
                    [0, 0, 0, 0, 1, 0],  # 3uy ← 3uy (element DOF 5)
                    [0, 0, 0, 0, 0, 1],  # 3rz ← 3rz (element DOF 6)
                    [0, 0, 0, 0, 0, 0],  # 4ux
                    [0, 0, 0, 0, 0, 0],  # 4uy
                    [0, 0, 0, 0, 0, 0],  # 4rz
                    [0, 0, 0, 0, 0, 0],  # 5ux
                    [0, 0, 0, 0, 0, 0],  # 5uy
                    [0, 0, 0, 0, 0, 0],  # 5rz
                    [0, 0, 0, 0, 0, 0],  # 6ux
                    [0, 0, 0, 0, 0, 0],  # 6uy
                    [0, 0, 0, 0, 0, 0],  # 6rz
                    [0, 0, 0, 0, 0, 0],  # 7ux
                    [0, 0, 0, 0, 0, 0],  # 7uy
                    [0, 0, 0, 0, 0, 0],  # 7rz
                ]),
                'second_moment_of_area': I, 'area': A, 'youngs_modulus': E, 'c': c,
            },
            {   # Element 2: Node 2 → 3 (diagonal)
                'start_node': 2, 'end_node': 3,
                'length': np.sqrt(2.5**2 + 1**2),
                'angle': 111.8,  # From node 2 to node 3
                'assembly_matrix': np.array([
                    [0, 0, 0, 1, 0, 0],  # 3ux ← 3ux (element DOF 4)
                    [0, 0, 0, 0, 1, 0],  # 3uy ← 3uy (element DOF 5)
                    [0, 0, 0, 0, 0, 1],  # 3rz ← 3rz (element DOF 6)
                    [0, 0, 0, 0, 0, 0],  # 4ux
                    [0, 0, 0, 0, 0, 0],  # 4uy
                    [0, 0, 0, 0, 0, 0],  # 4rz
                    [0, 0, 0, 0, 0, 0],  # 5ux
                    [0, 0, 0, 0, 0, 0],  # 5uy
                    [0, 0, 0, 0, 0, 0],  # 5rz
                    [0, 0, 0, 0, 0, 0],  # 6ux
                    [0, 0, 0, 0, 0, 0],  # 6uy
                    [0, 0, 0, 0, 0, 0],  # 6rz
                    [0, 0, 0, 0, 0, 0],  # 7ux
                    [0, 0, 0, 0, 0, 0],  # 7uy
                    [0, 0, 0, 0, 0, 0],  # 7rz
                ]),
                'second_moment_of_area': I, 'area': A, 'youngs_modulus': E, 'c': c,
            },
            {   # Element 3: Node 2 → 4 (vertical column)
                'start_node': 2, 'end_node': 4,
                'length': 2.5,
                'angle': 90,  # Vertical up
                'assembly_matrix': np.array([
                    [0, 0, 0, 0, 0, 0],  # 3ux
                    [0, 0, 0, 0, 0, 0],  # 3uy
                    [0, 0, 0, 0, 0, 0],  # 3rz
                    [0, 0, 0, 1, 0, 0],  # 4ux ← 4ux (element DOF 4)
                    [0, 0, 0, 0, 1, 0],  # 4uy ← 4uy (element DOF 5)
                    [0, 0, 0, 0, 0, 1],  # 4rz ← 4rz (element DOF 6)
                    [0, 0, 0, 0, 0, 0],  # 5ux
                    [0, 0, 0, 0, 0, 0],  # 5uy
                    [0, 0, 0, 0, 0, 0],  # 5rz
                    [0, 0, 0, 0, 0, 0],  # 6ux
                    [0, 0, 0, 0, 0, 0],  # 6uy
                    [0, 0, 0, 0, 0, 0],  # 6rz
                    [0, 0, 0, 0, 0, 0],  # 7ux
                    [0, 0, 0, 0, 0, 0],  # 7uy
                    [0, 0, 0, 0, 0, 0],  # 7rz
                ]),
                'second_moment_of_area': I, 'area': A, 'youngs_modulus': E, 'c': c,
            },
            {   # Element 4: Node 3 → 4 (horizontal beam)
                'start_node': 3, 'end_node': 4,
                'length': 1,
                'angle': 0,  # Horizontal right
                'assembly_matrix': np.array([
                    [1, 0, 0, 0, 0, 0],  # 3ux ← 3ux (element DOF 1)
                    [0, 1, 0, 0, 0, 0],  # 3uy ← 3uy (element DOF 2)
                    [0, 0, 1, 0, 0, 0],  # 3rz ← 3rz (element DOF 3)
                    [0, 0, 0, 1, 0, 0],  # 4ux ← 4ux (element DOF 4)
                    [0, 0, 0, 0, 1, 0],  # 4uy ← 4uy (element DOF 5)
                    [0, 0, 0, 0, 0, 1],  # 4rz ← 4rz (element DOF 6)
                    [0, 0, 0, 0, 0, 0],  # 5ux
                    [0, 0, 0, 0, 0, 0],  # 5uy
                    [0, 0, 0, 0, 0, 0],  # 5rz
                    [0, 0, 0, 0, 0, 0],  # 6ux
                    [0, 0, 0, 0, 0, 0],  # 6uy
                    [0, 0, 0, 0, 0, 0],  # 6rz
                    [0, 0, 0, 0, 0, 0],  # 7ux
                    [0, 0, 0, 0, 0, 0],  # 7uy
                    [0, 0, 0, 0, 0, 0],  # 7rz
                ]),
                'second_moment_of_area': I, 'area': A, 'youngs_modulus': E, 'c': c,
            },
            {   # Element 5: Node 3 → 5 (vertical column)
                'start_node': 3, 'end_node': 5,
                'length': 2.5,
                'angle': 90,  # Vertical up
                'assembly_matrix': np.array([
                    [1, 0, 0, 0, 0, 0],  # 3ux ← 3ux (element DOF 1)
                    [0, 1, 0, 0, 0, 0],  # 3uy ← 3uy (element DOF 2)
                    [0, 0, 1, 0, 0, 0],  # 3rz ← 3rz (element DOF 3)
                    [0, 0, 0, 0, 0, 0],  # 4ux
                    [0, 0, 0, 0, 0, 0],  # 4uy
                    [0, 0, 0, 0, 0, 0],  # 4rz
                    [0, 0, 0, 1, 0, 0],  # 5ux ← 5ux (element DOF 4)
                    [0, 0, 0, 0, 1, 0],  # 5uy ← 5uy (element DOF 5)
                    [0, 0, 0, 0, 0, 1],  # 5rz ← 5rz (element DOF 6)
                    [0, 0, 0, 0, 0, 0],  # 6ux
                    [0, 0, 0, 0, 0, 0],  # 6uy
                    [0, 0, 0, 0, 0, 0],  # 6rz
                    [0, 0, 0, 0, 0, 0],  # 7ux
                    [0, 0, 0, 0, 0, 0],  # 7uy
                    [0, 0, 0, 0, 0, 0],  # 7rz
                ]),
                'second_moment_of_area': I, 'area': A, 'youngs_modulus': E, 'c': c,
            },
            {   # Element 6: Node 3 → 6 (diagonal)
                'start_node': 3, 'end_node': 6,
                'length': np.sqrt(2.5**2 + 1**2),
                'angle': 68.2,  # From node 3 to node 6 (corrected angle)
                'assembly_matrix': np.array([
                    [1, 0, 0, 0, 0, 0],  # 3ux ← 3ux (element DOF 1)
                    [0, 1, 0, 0, 0, 0],  # 3uy ← 3uy (element DOF 2)
                    [0, 0, 1, 0, 0, 0],  # 3rz ← 3rz (element DOF 3)
                    [0, 0, 0, 0, 0, 0],  # 4ux
                    [0, 0, 0, 0, 0, 0],  # 4uy
                    [0, 0, 0, 0, 0, 0],  # 4rz
                    [0, 0, 0, 0, 0, 0],  # 5ux
                    [0, 0, 0, 0, 0, 0],  # 5uy
                    [0, 0, 0, 0, 0, 0],  # 5rz
                    [0, 0, 0, 1, 0, 0],  # 6ux ← 6ux (element DOF 4)
                    [0, 0, 0, 0, 1, 0],  # 6uy ← 6uy (element DOF 5)
                    [0, 0, 0, 0, 0, 1],  # 6rz ← 6rz (element DOF 6)
                    [0, 0, 0, 0, 0, 0],  # 7ux
                    [0, 0, 0, 0, 0, 0],  # 7uy
                    [0, 0, 0, 0, 0, 0],  # 7rz
                ]),
                'second_moment_of_area': I, 'area': A, 'youngs_modulus': E, 'c': c,
            },
            {   # Element 7: Node 4 → 6 (vertical column)
                'start_node': 4, 'end_node': 6,
                'length': 2.5,
                'angle': 90,  # Vertical up
                'assembly_matrix': np.array([
                    [0, 0, 0, 0, 0, 0],  # 3ux
                    [0, 0, 0, 0, 0, 0],  # 3uy
                    [0, 0, 0, 0, 0, 0],  # 3rz
                    [1, 0, 0, 0, 0, 0],  # 4ux ← 4ux (element DOF 1)
                    [0, 1, 0, 0, 0, 0],  # 4uy ← 4uy (element DOF 2)
                    [0, 0, 1, 0, 0, 0],  # 4rz ← 4rz (element DOF 3)
                    [0, 0, 0, 0, 0, 0],  # 5ux
                    [0, 0, 0, 0, 0, 0],  # 5uy
                    [0, 0, 0, 0, 0, 0],  # 5rz
                    [0, 0, 0, 1, 0, 0],  # 6ux ← 6ux (element DOF 4)
                    [0, 0, 0, 0, 1, 0],  # 6uy ← 6uy (element DOF 5)
                    [0, 0, 0, 0, 0, 1],  # 6rz ← 6rz (element DOF 6)
                    [0, 0, 0, 0, 0, 0],  # 7ux
                    [0, 0, 0, 0, 0, 0],  # 7uy
                    [0, 0, 0, 0, 0, 0],  # 7rz
                ]),
                'second_moment_of_area': I, 'area': A, 'youngs_modulus': E, 'c': c,
            },
            {   # Element 8: Node 5 → 6 (horizontal beam)
                'start_node': 5, 'end_node': 6,
                'length': 1,
                'angle': 0,  # Horizontal right
                'assembly_matrix': np.array([
                    [0, 0, 0, 0, 0, 0],  # 3ux
                    [0, 0, 0, 0, 0, 0],  # 3uy
                    [0, 0, 0, 0, 0, 0],  # 3rz
                    [0, 0, 0, 0, 0, 0],  # 4ux
                    [0, 0, 0, 0, 0, 0],  # 4uy
                    [0, 0, 0, 0, 0, 0],  # 4rz
                    [1, 0, 0, 0, 0, 0],  # 5ux ← 5ux (element DOF 1)
                    [0, 1, 0, 0, 0, 0],  # 5uy ← 5uy (element DOF 2)
                    [0, 0, 1, 0, 0, 0],  # 5rz ← 5rz (element DOF 3)
                    [0, 0, 0, 1, 0, 0],  # 6ux ← 6ux (element DOF 4)
                    [0, 0, 0, 0, 1, 0],  # 6uy ← 6uy (element DOF 5)
                    [0, 0, 0, 0, 0, 1],  # 6rz ← 6rz (element DOF 6)
                    [0, 0, 0, 0, 0, 0],  # 7ux
                    [0, 0, 0, 0, 0, 0],  # 7uy
                    [0, 0, 0, 0, 0, 0],  # 7rz
                ]),
                'second_moment_of_area': I, 'area': A, 'youngs_modulus': E, 'c': c,
            },
            {   # Element 9: Node 5 → 7 (vertical column)
                'start_node': 5, 'end_node': 7,
                'length': 2.5,  # Corrected length
                'angle': 90,    # Vertical up
                'assembly_matrix': np.array([
                    [0, 0, 0, 0, 0, 0],  # 3ux
                    [0, 0, 0, 0, 0, 0],  # 3uy
                    [0, 0, 0, 0, 0, 0],  # 3rz
                    [0, 0, 0, 0, 0, 0],  # 4ux
                    [0, 0, 0, 0, 0, 0],  # 4uy
                    [0, 0, 0, 0, 0, 0],  # 4rz
                    [1, 0, 0, 0, 0, 0],  # 5ux ← 5ux (element DOF 1)
                    [0, 1, 0, 0, 0, 0],  # 5uy ← 5uy (element DOF 2)
                    [0, 0, 1, 0, 0, 0],  # 5rz ← 5rz (element DOF 3)
                    [0, 0, 0, 0, 0, 0],  # 6ux
                    [0, 0, 0, 0, 0, 0],  # 6uy
                    [0, 0, 0, 0, 0, 0],  # 6rz
                    [0, 0, 0, 1, 0, 0],  # 7ux ← 7ux (element DOF 4)
                    [0, 0, 0, 0, 1, 0],  # 7uy ← 7uy (element DOF 5)
                    [0, 0, 0, 0, 0, 1],  # 7rz ← 7rz (element DOF 6)
                ]),
                'second_moment_of_area': I, 'area': A, 'youngs_modulus': E, 'c': c,
            },
            {   # Element 10: Node 6 → 7 (diagonal)
                'start_node': 6, 'end_node': 7,
                'length': np.sqrt(2.5**2 + 1**2),
                'angle': 111.8,  # From node 6 to node 7
                'assembly_matrix': np.array([
                    [0, 0, 0, 0, 0, 0],  # 3ux
                    [0, 0, 0, 0, 0, 0],  # 3uy
                    [0, 0, 0, 0, 0, 0],  # 3rz
                    [0, 0, 0, 0, 0, 0],  # 4ux
                    [0, 0, 0, 0, 0, 0],  # 4uy
                    [0, 0, 0, 0, 0, 0],  # 4rz
                    [0, 0, 0, 0, 0, 0],  # 5ux
                    [0, 0, 0, 0, 0, 0],  # 5uy
                    [0, 0, 0, 0, 0, 0],  # 5rz
                    [1, 0, 0, 0, 0, 0],  # 6ux ← 6ux (element DOF 1)
                    [0, 1, 0, 0, 0, 0],  # 6uy ← 6uy (element DOF 2)
                    [0, 0, 1, 0, 0, 0],  # 6rz ← 6rz (element DOF 3)
                    [0, 0, 0, 1, 0, 0],  # 7ux ← 7ux (element DOF 4)
                    [0, 0, 0, 0, 1, 0],  # 7uy ← 7uy (element DOF 5)
                    [0, 0, 0, 0, 0, 1],  # 7rz ← 7rz (element DOF 6)
                ]),
                'second_moment_of_area': I, 'area': A, 'youngs_modulus': E, 'c': c,
            }
        ],

        # LOADS: Define distributed and point loads on elements
        'loads': [
            # Example loads - uncomment and modify as needed:

            # {'type': 'UDL', 'element': 5, 'magnitude': -19600.7, 'direction': 'local'},  # 5 kN/m downward on beam 3-4
            # {'type': 'UDL', 'element': 9, 'magnitude': -19600.7, 'direction': 'local'},  # 5 kN/m downward on beam 3-4

            {'type': 'UDL', 'element': 5, 'magnitude': -4880, 'direction': 'local'},  # 5 kN/m downward on beam 3-4
            {'type': 'UDL', 'element': 9, 'magnitude': -4880, 'direction': 'local'},  # 5 kN/m downward on beam 3-4

            # {'type': 'PDL', 'element': 8, 'magnitude': -10000, 'direction': 'vertical', 'a': 0.5},  # 10 kN downward at midspan of beam 5-6
        ],
        
        # NODAL FORCES: Direct forces/moments applied at nodes (in global DOF order)
        # Order: [3ux, 3uy, 3rz, 4ux, 4uy, 4rz, 5ux, 5uy, 5rz, 6ux, 6uy, 6rz, 7ux, 7uy, 7rz]
        'normal_forces': [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0],
        
        # DISPLAY OPTIONS: Control output format and visualization
        'display_magnitude': 100,        # Deflection magnification factor for plotting
        'show_matrices': False,          # True = vertical matrix display, False = compact display
        'show_deflection_labels': True   # True = show deflection values on plot
    }

    results = analyzer.run_analysis(config)

if __name__ == "__main__":
    main()