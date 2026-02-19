# RSTE_budget

A comprehensive OpenFOAM function object library for calculating and analysing the Reynolds Stress Transport Equation (RSTE) budget in Large Eddy Simulation (LES).

## Overview

This repository provides a validated implementation of all terms in the Reynolds Stress Transport Equation budget for OpenFOAM-based LES simulations. The library enables detailed turbulence analysis by computing convection, production, turbulent transport, pressure-diffusion, pressure-strain, viscous diffusion, and dissipation terms, along with subgrid-scale (SGS) contributions.

### Key Features

- **Complete RSTE Budget**: Implementation of all terms in the Reynolds Stress Transport Equation
- **Validated Implementation**: Rigorously tested against DNS data for canonical flows
- **OpenFOAM Integration**: Seamless integration with OpenFOAM's function object framework
- **Cartesian Coordinates**: Full implementation in Cartesian coordinate system
- **SGS Contributions**: Subgrid-scale term calculations for LES - used as a residual term value
- **DNS Compatibility**: Can be run in DNS mode by setting turbulence model to "laminar"
- **Memory Efficient**: Optimised field calculations using OpenFOAM utilities
- **Parallel Support**: Full parallel computation support with HPC considerations

## Table of Contents

- [Installation](#installation)
- [Quick Start](#quick-start)
- [Project Structure](#project-structure)
- [Usage](#usage)
- [Validation](#validation)
- [Implementation Details](#implementation-details)
- [Examples](#examples)
- [Contributing](#contributing)
- [License](#license)

## Installation

### Prerequisites

- OpenFOAM (tested with recent version OpenFOAM-v2506)
- Git (for version control)
- Python 3.x (for parameter setup scripts)

### Setup Instructions

1. **Clone the repository**:
   ```bash
   git clone https://github.com/AUfluids/RSTE_budget.git
   cd RSTE_budget
   ```

2. **Copy the RSTE calculation codes**:
   ```bash
   # Copy the calcRSTEcode directory to your OpenFOAM case's system directory
   cp -r calcRSTEcode/ /path/to/your/case/system/
   ```

3. **Configure your OpenFOAM case**:
   - Add the RSTE function objects to your `system/controlDict`
   - Ensure proper field averaging setup
   - Configure parallel execution if needed

## Quick Start

### Running the Example Case

1. **Navigate to the example**:
   ```bash
   cd examples/channel_flow_Re_tau_180
   ```
   
2. **Run the simulation**:
   ```bash
   # Make the script executable
   chmod +x Allrun
   
   # Run the case
   ./Allrun
   ```

3. **View results**:
   - Check the `logs/` directory for simulation output
   - RSTE budget terms are written to the case directory
   - Use ParaView or other post-processing tools to visualise results

### Basic Usage

The RSTE budget calculation is integrated into your OpenFOAM case through function objects. Add the following to your `system/controlDict`:

```cpp
functions
{
    #include "calcRSTEcode/fieldAverageRSTE"
    #include "calcRSTEcode/calcProduction"
    #include "calcRSTEcode/calcConvection"
    #include "calcRSTEcode/calcTurbulentTransport"
    #include "calcRSTEcode/calcPressureDiffusion"
    #include "calcRSTEcode/calcPressureStrain"
    #include "calcRSTEcode/calcViscousDiffusion"
    #include "calcRSTEcode/calcDissipation"
    #include "calcRSTEcode/calcSGS"
}
```

## Project Structure

```
RSTE_budget/
├── calcRSTEcode/                    # Main RSTE calculation codes
│   ├── auxiliaryFunctionsRSTE      # Helper functions
│   ├── calcConvection              # Convective term calculation
│   ├── calcDissipation             # Dissipation term calculation
│   ├── calcPressureDiffusion       # Pressure diffusion term
│   ├── calcPressureStrain          # Pressure-strain correlation
│   ├── calcProduction              # Production term calculation
│   ├── calcSGS                     # Subgrid-scale term calculation
│   ├── calcTurbulentTransport      # Turbulent transport term
│   ├── calcViscousDiffusion        # Viscous diffusion term
│   └── fieldAverageRSTE            # Field averaging configuration
├── examples/
│   └── channel_flow_Re_tau_180/    # Example channel flow case
│       ├── 0.orig/                 # Initial conditions
│       ├── constant/               # Physical properties
│       ├── system/                 # OpenFOAM configuration
│       │   ├── calcRSTEcode/       # Case-specific RSTE codes
│       │   ├── controlDict         # Solver control
│       │   ├── fvSchemes           # Numerical schemes
│       │   └── fvSolution          # Solver settings
│       ├── pythonScripts/          # Parameter setup scripts
│       ├── Allrun                  # Run script
│       └── Allclean                # Cleanup script
├── publication_description.tex     # Scientific publication details
├── LICENSE                         # License information
└── README.md                       # This file
```

## Usage

### Integration with Your Case

1. **Copy RSTE codes**:
```bash
   cp -r calcRSTEcode/ /path/to/your/case/system/
   ```

2. **Modify controlDict**:
   Add the RSTE function objects to your `system/controlDict`:

   ```cpp
   functions
   {
       #include "calcRSTEcode/fieldAverageRSTE"
       // Add other terms as needed
       #include "calcRSTEcode/calcProduction"
       #include "calcRSTEcode/calcConvection"
       // ... etc
   }
   ```

3. **Configure field averaging**:
   The `fieldAverageRSTE` function object sets up the necessary field averaging for:
   - Mean velocity components
   - Reynolds stress components
   - Higher-order correlations
   - Pressure-velocity correlations

4. **Run your simulation**:
   ```bash
   # Serial execution
   pimpleFoam
   
   # Parallel execution
   mpirun -np 4 pimpleFoam -parallel
   ```

### DNS Mode

The RSTE budget calculation can also be run in DNS mode. To run in DNS mode:

1. **Set turbulence model to laminar**:
   ```cpp
   // In constant/turbulenceProperties
   simulationType  laminar;
   ```

2. **Run the simulation**:
   The same execution commands apply, but the SGS terms will be zero since no subgrid-scale model is active.

**Note**: DNS mode is primarily useful for validation purposes or when comparing with DNS data, but requires significantly higher computational resources and finer meshes.

### HPC Considerations

For complex HPC systems, add the following to your `system/controlDict` to avoid compilation issues:

```cpp
runTimeModifiable true;

OptimisationSwitches
{
    fileModificationSkew 10; // Increase from default 1-2
    // maxFileModificationPolls 100; // Optional: increase from default 50
}
```

### Output Fields

The library generates the following output fields:

- **Reynolds Stress Components**: `RMean` (symmetric tensor)
- **Production Term**: `P_ij_RSTE` (symmetric tensor)
- **Convective Term**: `C_ij_RSTE` (symmetric tensor)
- **Turbulent Transport**: `T_ij_RSTE` (symmetric tensor)
- **Pressure Diffusion**: `Dp_ij_RSTE` (symmetric tensor)
- **Pressure-Strain**: `Phi_ij_RSTE` (symmetric tensor)
- **Viscous Diffusion**: `D_ij_RSTE` (symmetric tensor)
- **Dissipation**: `eps_ij_RSTE` (symmetric tensor)
- **SGS Term**: `SGS_ij_RSTE` (symmetric tensor)

## Validation

The implementation has been rigorously validated against high-fidelity DNS data for canonical wall-bounded turbulent flows:

### Validation Cases

1. **Channel Flow** (Re_τ = 180)
   - Comparison with DNS data from Lee & Moser (2015)
   - Mesh refinement study (L₁ to L₄)
   - All RSTE budget terms validated

2. **Pipe Flow** (Re_τ = 180) - Not included in this repository - Available in publication
   - Comparison with DNS data from El Khoury et al. (2013)
   - TKE budget validation (coordinate-invariant)
   - Multiple mesh resolutions

### Key Findings

- **Systematic Convergence**: Results converge systematically towards DNS data with mesh refinement
- **Physical Accuracy**: All budget terms show correct physical behaviour
- **Numerical Robustness**: Implementation handles near-wall regions correctly
- **SGS Validation**: Subgrid-scale contributions diminish appropriately with mesh refinement

## Implementation Details

### Theoretical Background

The Reynolds Stress Transport Equation in Cartesian coordinates:

```
∂R_ij/∂t = C_ij + P_ij + T_ij + D_p,ij + Φ_ij + D_ν,ij + ε_ij + Π_ij*
```

Where:
- **C_ij**: Convective term
- **P_ij**: Production term
- **T_ij**: Turbulent transport term
- **D_p,ij**: Pressure diffusion term
- **Φ_ij**: Pressure-strain correlation
- **D_ν,ij**: Viscous diffusion term
- **ε_ij**: Viscous dissipation term
- **Π_ij***: Subgrid-scale term (LES only)

### Numerical Implementation

- **Field Expansion**: Higher-order correlations expanded for on-the-fly averaging
- **Conservative Form**: Convective terms implemented in conservative form
- **Gradient Calculations**: Uses OpenFOAM's `grad` and `div` functions
- **Boundary Conditions**: Proper handling of wall boundary conditions
- **Memory Management**: Efficient field creation and management

### Coordinate Systems

- **Current**: Cartesian coordinates (x, y, z)
- **Future**: Extension to cylindrical and spherical coordinates planned
- **TKE Budget**: Coordinate-invariant (trace of RSTE budget)

## Examples

### Channel Flow Example

The included channel flow example demonstrates:

- **Setup**: Complete OpenFOAM case configuration
- **Mesh**: Structured hexahedral mesh with wall refinement
- **Solver**: pimpleFoam with WALE SGS model
- **Initial Conditions**: Python script with automatic initial conditions
- **Analysis**: Complete RSTE budget calculation

**Important Note**: Since the simulation starts from plug-flow conditions, a certain initial perturbation of the flow is necessary to onset turbulence and obtain a fully-developed turbulent channel flow. The example case DOES NOT includes appropriate initial perturbations to ensure proper turbulent development nor a script to post-process the RSTE budget data.

### Running the Example

```bash
cd examples/channel_flow_Re_tau_180

# Run the simulation
./Allrun

# Clean up (optional)
./Allclean
```

### Custom Cases

To adapt for your own cases:

1. **Modify mesh**: Update `system/blockMeshDict`
2. **Adjust parameters**: Edit `pythonScripts/paramSetup.py`
3. **Configure solver**: Modify `system/fvSolution` and `system/fvSchemes`
4. **Set initial conditions**: Update `0.orig/` directory

## Contributing

### Getting Started

1. Fork the repository
2. Create a feature branch: `git checkout -b feature/new-feature`
3. Make your changes
4. Test with the example case
5. Submit a pull request

### Development Guidelines

- **Code Style**: Follow OpenFOAM coding conventions
- **Documentation**: Update comments and documentation
- **Testing**: Validate against known solutions
- **Performance**: Consider memory and computational efficiency

### Areas for Contribution

- **Coordinate Systems**: Extension to cylindrical/spherical coordinates
- **Memory Optimisation**: Reduced memory footprint
- **Additional Terms**: Implementation of additional turbulence terms
- **Validation**: New validation cases and comparisons
- **Documentation**: Improved user guides and examples
- **OpenFOAM Library**: Future implementation as a proper OpenFOAM library to avoid per-case compilation

## License

This project is licensed under the Creative Commons License - see the [LICENSE](LICENSE) file for details.

## Authors

- **Mario Javier Rincón** - *Lead Developer* - Aarhus University
- **Christoffer Hansen** - *Lead-Developer* - Aarhus University

## Acknowledgments

- **Aarhus University** - Mechanical and Production Engineering Department
- **Kamstrup A/S** - Fluid Mechanics R&D Group in Quality and Sustainability Department
- **OpenFOAM Community** - For the excellent CFD framework
- **DNS Data Providers** - Lee & Moser, El Khoury et al. for validation data
- **Scientific Community** - For turbulence modelling research

## References

- Rincón, M. J., Hansen, C., Reclari, M., & Abkar, M. (2025). Implementation and verification of the Reynolds stress transport equations in OpenFOAM. *OpenFOAM Journal*, UNDER REVIEW.
- Lee, M. & Moser, R.D. (2015). Direct numerical simulation of turbulent channel flow up to Re_τ ≈ 5200. *Journal of Fluid Mechanics*, 774, 395-415.
- El Khoury, G.K., Schlatter, P., Noorani, A., Fischer, P.F., & Brethouwer, G. (2013). Direct numerical simulation of turbulent pipe flow at moderately high Reynolds numbers. *Flow, Turbulence and Combustion*, 91(3), 475-495.
- Pope, S.B. (2000). *Turbulent Flows*. Cambridge University Press.

## Related Work

This implementation supports research in:
- **Turbulence Modelling**: Development and validation of RANS models
- **LES Analysis**: Understanding subgrid-scale effects
- **Wall-Bounded Flows**: Analysis of near-wall turbulence physics
- **Educational Purposes**: Teaching turbulence physics and CFD

---

**Note**: This repository provides a validated tool for RSTE budget analysis in OpenFOAM. For questions or issues, please open an issue on GitHub or contact the authors.