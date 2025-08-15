# CODE FOR ENME302 TEST
# This code will be used for the ENME302 Test - August 19th 2025

# The test will use frame elements
# There will be no more than four frame elements per question. Each problem will have either 2, 3, or 4 frame elements.

# The structure of the code will be similar to the tutorial code. I will build methods for 2, 3, and 4 frame elements,
# and whichever one is needed will be used in the test. Each file will ask for the element properties, angles, and assembly matrices.
# The subsequent calculations and plots will be completed in each file.

# Some notes:


"""
LINEAR ALGEBRA IN FINITE ELEMENT ANALYSIS (FEA)
================================================

1. COORDINATE SYSTEMS
---------------------
- Local coordinate: Element-specific coordinate system.
- Global coordinate: System-wide coordinate system.

Key Concept:
    Local → Transformation Matrix (Aᵉ) → Global
2. VECTORS
----------
- Local vector: Defined in the element's local coordinate system.
- Global vector: Defined in the global coordinate system.
    - A node is a point where force acts.

3. STIFFNESS MATRICES
---------------------
- Kᵉ: Local stiffness matrix (element level).
- K̂ᵉ: Transformed stiffness matrix in global coordinates.
    - K̂ᵉ = (Aᵉ)ᵀ · Kᵉ · Aᵉ

- K_G: Global stiffness matrix for the entire structure.
    - K_G = Σ (Aᵉ)ᵀ · Kᵉ · Aᵉ

4. ASSEMBLY MATRIX
------------------
- Aᶜ: Assembly matrix
    - Maps element matrices into the global system.

5. DISPLACEMENT VECTORS
- D: Global displacement vector.
- dᵉ: Local displacement vector.
    - dᵉ = Aᵉ · D

5. FORCE VECTORS
----------------
- Q: External force vector (global)
    - Q = K_G · q

- F′: Element forces at each node (global)
    - F′ = K · (A′)ᵀ · q

- d′: Element displacements (global)
    - d′ = A′ · (A′)ᵀ · q

- f′: Element forces at each node (local)
    - f′ = K · d′
"""






"""
ENME302 - Finite Element Analysis (FEA) - Course Notes
University of Canterbury, Semester 2, 2025
Instructor: Prof. Geoff Rodgers

===========================================================
1. INTRODUCTION TO FEA
===========================================================

Finite Element Analysis (FEA) is a numerical method for solving complex structural problems by discretizing them into simpler elements.

Key Steps:
-----------
A. Choose element type:
    - Truss (axial only)
    - Beam (bending and shear)
    - Frame (axial + bending + shear)

B. Identify boundary conditions:
    - Fixed, Free, Constrained

C. Apply loads:
    - Point forces, Moments, Distributed loads

D. Define Degrees of Freedom (DOFs):
    - Locations where displacements are solved

===========================================================
2. SYMBOLS AND DEFINITIONS
===========================================================

- f: Vector of applied forces
- d: Vector of nodal displacements
- K: Stiffness matrix
- E: Young’s modulus
- A: Cross-sectional area
- L: Element length
- I: Second moment of area
- σ: Stress
- ε: Strain
- u(x): Displacement function
- ψ(x): Shape function
- θ: Angle between global and local axes
- Λ: Transformation matrix
- q: Global DOF vector
- Q: Global force vector

===========================================================
3. BAR (TRUSS) ELEMENT FORMULATION
===========================================================

Displacement:
    u(x) = ψ₁(x)·d₁ + ψ₂(x)·d₂
    ψ₁(x) = 1 - x/L
    ψ₂(x) = x/L

Strain:
    ε(x) = (d₂ - d₁) / L

Stress:
    σ(x) = E * ε(x)

Stiffness Matrix:
    K = (EA / L) * [[1, -1], [-1, 1]]

===========================================================
4. ENERGY METHODS
===========================================================

Work-Energy Method:
    U = W = (1/2)·P·δ

Virtual Work Principle:
    1·Δ = Σ(fᵢ·Fᵢ·Lᵢ) / (Aᵢ·Eᵢ)

===========================================================
5. TRANSFORMATION BETWEEN COORDINATES
===========================================================

Transformation Matrix (Λ):
    Λ = [[cosθ, sinθ, 0, 0],
         [0, 0, cosθ, sinθ]]

Transformations:
    d_local = Λ · d_global
    f_global = Λᵀ · f_local
    K_global = Λᵀ · K_local · Λ

===========================================================
6. GLOBAL SYSTEM ASSEMBLY
===========================================================

Assembly Matrix (Aᵉ): Maps local DOFs to global DOFs

Global Stiffness Matrix:
    K_G = Σ Aᵉᵀ · K̂ᵉ · Aᵉ

Global Force Vector:
    Q = Σ Aᵉᵀ · fᵉ

Solve:
    Q = K_G · q

===========================================================
7. BEAM ELEMENTS
===========================================================

Assumptions:
    - Slender beams (L >> cross-section)
    - Small deflections
    - No axial deformation

Shape Functions (Hermite):
    v(x) = H₁(x)·d₁ + H₂(x)·θ₁ + H₃(x)·d₂ + H₄(x)·θ₂

Stiffness Matrix:
    K_beam = (EI / L³) * [
        [12, 6L, -12, 6L],
        [6L, 4L², -6L, 2L²],
        [-12, -6L, 12, -6L],
        [6L, 2L², -6L, 4L²]
    ]

===========================================================
8. FRAME ELEMENTS
===========================================================

Combines axial and flexural behavior:
    - 6 DOF: [u₁, v₁, θ₁, u₂, v₂, θ₂]

Frame Stiffness Matrix:
    K_frame = K_bar + K_beam (superposition)

Transformation:
    K̂_frame = Λᵀ · K_frame · Λ

===========================================================
9. PYTHON IMPLEMENTATION STRATEGY
===========================================================

Functions:
    - local_bar(E, A, L)
    - global_bar(K_local, θ)
    - local_beam(E, I, L)
    - global_beam(K_local, θ)
    - assemble_global(K_list, A_list)
    - solve_displacements(KG, Q)
    - compute_internal_forces(K_local, d_local)

===========================================================
10. FINAL REMARKS
===========================================================

- FEA is powerful for solving complex structures.
- Understanding the formulation helps debug and validate software results.
- Always check boundary conditions and constraints.
- For large deformations, iterative nonlinear analysis is required.
"""
