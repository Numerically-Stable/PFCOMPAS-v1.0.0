# PFCOMPAS-v1.0.0: Phase-Field Modeling of Fracture

PFCOMPAS-v1.0.0 is a self-contained, high-performance MATLAB software package for quasi-static brittle fracture based on the regularized AT2 functional. Designed for structural engineers and computational mechanicians, the framework bypasses standard interpreted-language limitations via a parallel map-reduce assembly engine utilizing precomputed sparsity patterns, achieving near-compiled execution rates.

## Core Features
* **Constitutive Asymmetry:** Natively incorporates isotropic, Amor volumetric-deviatoric, and Miehe spectral tension-compression splits.
* **Algorithmic Rigor:** Deploys hyperbolically smoothed Macaulay brackets with dynamic scale invariance to provide exact, $C^{1}$-continuous consistent tangents. 
* **Cohesive Zone Library:** Integrates 10 distinct degradation functions to artificially delay crack nucleation and simulate threshold-based cohesive zone behaviors.
* **Spatial Diagnostics:** Features an explicit spatial diagnostic suite based on Verfürth a posteriori error estimators, which evaluate localized traction and damage gradient interface flux jumps across interior edges.
* **Adaptive Controller:** An adaptive predictive load-stepping heuristic dynamically scales load increments to bypass manual tuning across stiff pre-cracking and unstable post-peak softening regimes.

```mermaid
flowchart TD

    %% ================================================================
    %% GLOBAL INITIALIZATION
    %% ================================================================
    A([Start]) --> B[Preprocessing and<br/>Initialization]

    B --> C[Initialize Solution State<br/>
            u₀, φ₀, H₀]

    %% ================================================================
    %% OUTER ADAPTIVE LOAD-STEP LOOP
    %% ================================================================
    C --> D[Accepted State<br/>
            uₙ, φₙ, Hₙ, λₙ]

    D --> E[Propose Trial Load Step<br/>
            λₙ₊₁ᵗʳ = λₙ + Δλ]

    %% ================================================================
    %% STAGGERED SOLVER
    %% ================================================================
    E --> F[Staggered Solver]

    F --> G[Displacement<br/>Newton Solve]

    G --> H[History Field<br/>Update]

    H --> I[Phase-Field<br/>Newton Solve]

    I --> J{Staggered<br/>Converged?}

    %% ================================================================
    %% STAGGERED NON-CONVERGENCE
    %% ================================================================
    J -- No --> K[Reduce Δλ<br/>and Retry]

    K -->|Δλ = Δλ / 2| D

    %% ================================================================
    %% POST-PROCESSING
    %% ================================================================
    J -- Yes --> L[Trial Post-Processing<br/>
                    E, RF and Diagnostics]

    L --> M{max Δφ ≤ Δφₘₐₓ?}

    %% ================================================================
    %% EXCESSIVE PHASE-FIELD EVOLUTION
    %% ================================================================
    M -- No --> K

    %% ================================================================
    %% ACCEPTANCE
    %% ================================================================
    M -- Yes --> N[Accept Trial State<br/>
                    uₙ₊₁, φₙ₊₁, Hₙ₊₁]

    N --> O[Adaptive Update of Δλ<br/>
            Based on Solver Difficulty<br/>
            and max Δφ]

    %% ================================================================
    %% OUTPUT / DIAGNOSTICS
    %% ================================================================
    O --> P[Store Accepted-Step Data<br/>
            Stress, Error Indicators,<br/>
            Energies, Energy Balance<br/>
            and ParaView Output]

    P --> Q[Next Load Step]

    %% ================================================================
    %% OUTER LOAD-STEP RETURN
    %% ================================================================
    Q --> D

    %% ================================================================
    %% STYLING
    %% ================================================================
    classDef startEnd fill:#E8F5E9,stroke:#2E7D32,stroke-width:2px,color:#1B1B1B;
    classDef initialization fill:#E3F2FD,stroke:#1565C0,stroke-width:1.5px,color:#1B1B1B;
    classDef accepted fill:#E8F5E9,stroke:#388E3C,stroke-width:2px,color:#1B1B1B;
    classDef trial fill:#FFF8E1,stroke:#F9A825,stroke-width:1.5px,color:#1B1B1B;
    classDef solver fill:#EDE7F6,stroke:#6A1B9A,stroke-width:1.5px,color:#1B1B1B;
    classDef decision fill:#FFF3E0,stroke:#EF6C00,stroke-width:2px,color:#1B1B1B;
    classDef rejection fill:#FFEBEE,stroke:#C62828,stroke-width:2px,color:#1B1B1B;
    classDef output fill:#E0F2F1,stroke:#00695C,stroke-width:1.5px,color:#1B1B1B;

    class A startEnd;
    class B,C initialization;
    class D,N accepted;
    class E trial;
    class F,G,H,I solver;
    class J,M decision;
    class K rejection;
    class L,O,P,Q output;
```
## Requirements
* MATLAB (> R2021a)
* Parallel Computing Toolbox

## Installation & Quick Start
1. Clone the repository to your local machine:
   `git clone https://github.com/NUmerically-Stable/PF-COMPAS.git`
2. Add the repository directory and its subfolders to your MATLAB path.
3. **Run the Minimal Working Example (MWE):**
   Open MATLAB and execute the Single Edge Notched Tension (SENT) benchmark:
   `main.m`
   *This script autonomously loads the mesh (make sure to add it to path), executes the adaptive staggered active-set Newton solver, and streams output data to the disk.*

## Visualization
State variables are streamed directly to disk as XML Unstructured Grid files (`.vtu`) and registered within a master temporal collection wrapper (`.pvd`), preserving workstation memory. Open the `.pvd` file in ParaView to analyze the chronological evolution of the damage field and Verfürth error indicators.

## Documentation
Comprehensive theoretical formulations, Jacobian derivations, and algorithmic flowcharts are available in the [PFCOMPAS-v1.0.0 Documentation Report](docs/PFCOMPAS-v1.0.0_Documentation.pdf).

## License
This project is licensed under the MIT License - see the `LICENSE` file for details. 

## Contact
For questions or support, contact Pranjal Saxena (pranjals21@iitk.ac.in or contactpranjal02@gmail.com).
