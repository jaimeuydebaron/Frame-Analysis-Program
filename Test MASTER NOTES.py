"""
ENME302 - Finite Element Analysis (FEA) Test Notes
University of Canterbury, Semester 2, 2025
Instructor: Prof. Geoff Rodgers

===========================================================
1. INTRODUCTION TO FEA
===========================================================

Finite Element Analysis (FEA) is a numerical method for solving
structural problems by discretizing them into smaller elements.

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

- f   : Vector of applied forces
- d   : Vector of nodal displacements
- K   : Stiffness matrix
- E   : Young’s modulus
- A   : Cross-sectional area
- L   : Element length
- I   : Second moment of area
- σ   : Stress
- ε   : Strain
- u(x): Displacement function
- ψ(x): Shape function
- θ   : Angle between global and local axes
- Λ   : Transformation matrix
- q   : Global DOF vector
- Q   : Global force vector
- Aᵉ  : Transformation matrix (local → global)
- Aᶜ  : Assembly matrix (element → global)

===========================================================
3. COORDINATE SYSTEMS & TRANSFORMATIONS
===========================================================

- Local coordinate: Element-specific system.
- Global coordinate: Structure-wide system.

Transformations:
    d_local  = Aᵉ · d_global
    f_global = (Aᵉ)ᵀ · f_local
    K_global = (Aᵉ)ᵀ · K_local · Aᵉ

Angle definition:
    - Start from the global x-axis (Xg).
    - Rotate anticlockwise to the element local x-axis (Xe).

===========================================================
4. ELEMENT FORMULATIONS
===========================================================

(A) Bar (Truss) Element:
------------------------
Displacement:
    u(x) = ψ₁(x)·d₁ + ψ₂(x)·d₂
    ψ₁(x) = 1 - x/L
    ψ₂(x) = x/L

Strain:
    ε(x) = (d₂ - d₁) / L

Stress:
    σ(x) = E · ε(x)

Stiffness Matrix:
    K_bar = (EA / L) * [[1, -1], [-1, 1]]

(B) Beam Element:
-----------------
Assumptions:
    - Slender beams (L >> cross-section)
    - Small deflections
    - No axial deformation

Shape Functions (Hermite):
    v(x) = H₁(x)·d₁ + H₂(x)·θ₁ + H₃(x)·d₂ + H₄(x)·θ₂

Stiffness Matrix:
    K_beam = (EI / L³) * [
        [12,   6L,  -12,   6L],
        [6L,  4L²,  -6L,  2L²],
        [-12, -6L,   12,  -6L],
        [6L,  2L²,  -6L,  4L²]
    ]

(C) Frame Element:
------------------
- Combines axial and flexural behavior.
- 6 DOF per element: [u₁, v₁, θ₁, u₂, v₂, θ₂]

Frame Stiffness:
    K_frame = K_bar + K_beam
    K̂_frame = (Aᵉ)ᵀ · K_frame · Aᵉ

===========================================================
5. LOADING TYPES
===========================================================

(A) Uniformly Distributed Load (UDL):
    Q_UDL = [0, wL/2, wL²/12, 0, wL/2, -wL²/12]

(B) Linearly Varying Load (LVL):
    Peak at x = L:
        Q_LVL = [0, 3wL/20, wL²/30, 0, 7wL/20, -wL²/20]
    Peak at x = 0:
        Q_LVL = [0, 7wL/20, wL²/20, 0, 3wL/20, -wL²/30]

(C) Mid-span Point Load:
    Q_PL = [0, w/2, wL/8, 0, w/2, -wL/8]

(D) Distributed Axial Load
    Q_DAL = [wL/2, 0, 0, -wL/2, 0, 0]

(E) Concentrated Axial Load
    Q_CAL = [0, 0, 0, 0, 0, 0]

===========================================================
6. GLOBAL SYSTEM ASSEMBLY
===========================================================

- Kᵉ   : Local element stiffness matrix
- K̂ᵉ  : Element stiffness matrix in global coordinates
          K̂ᵉ = (Aᵉ)ᵀ · Kᵉ · Aᵉ

- K_G  : Global stiffness matrix for structure
          K_G = Σ (Aᵉ)ᵀ · Kᵉ · Aᵉ

- Aᶜ   : Assembly matrix, maps element matrices into global system.

Vectors:
    - D  : Global displacement vector
    - dᵉ : Local displacement vector  →  dᵉ = Aᵉ · D
    - Q  : External global force vector  →  Q = K_G · q
    - f′ : Element forces in local coords
    - F′ : Element forces in global coords

===========================================================
7. NODE & ELEMENT NUMBERING
===========================================================

1. Node numbering:
   - Top → bottom, then left → right.

2. Element numbering:
   - Choose consistent scheme (left → right or top → bottom).
   - Element local axis always goes Node 1 → Node 2.

3. DOF Alignment:
   - Each qi aligns with the same DOF across all connected elements.

===========================================================
8. ENERGY METHODS
===========================================================

Work-Energy Method:
    U = W = (1/2)·P·δ

Virtual Work Principle:
    1·Δ = Σ(fᵢ·Fᵢ·Lᵢ) / (Aᵢ·Eᵢ)

===========================================================
9. PYTHON IMPLEMENTATION STRATEGY
===========================================================

Functions to build:
    - local_bar(E, A, L)
    - global_bar(K_local, θ)
    - local_beam(E, I, L)
    - global_beam(K_local, θ)
    - assemble_global(K_list, A_list)
    - solve_displacements(KG, Q)
    - compute_internal_forces(K_local, d_local)

===========================================================
10. QUICK REMINDERS
===========================================================

- "q" is always on the left-hand side of matrix equations (displacements).
- "D" is associated with element DOFs (roof of equation).
- Align DOFs across elements at shared nodes.
- Check numbering consistency for nodes, elements, and axes.
- For nonlinear problems (large deflections), iterative methods required.
"""
