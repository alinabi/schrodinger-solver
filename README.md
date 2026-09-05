# Schrödinger Equation Solver

A comprehensive Haskell library for solving the time-dependent Schrödinger equation for systems of N charged particles with arbitrary masses, evolving in time-dependent electromagnetic fields with full particle-particle Coulomb interactions.

## Mathematical Background

The time-dependent Schrödinger equation:

```
iℏ ∂ψ/∂t = Ĥ(t) ψ
```

where the Hamiltonian for N particles is:

```
Ĥ(t) = Σᵢ [-ℏ²/(2mᵢ) ∇ᵢ² + qᵢ φ(rᵢ, t) - qᵢ/c (vᵢ · A(rᵢ, t))]
     + Σᵢ<ⱼ [qᵢ qⱼ / (4πε₀|rᵢ - rⱼ|)]
```

### Key Components

1. **Kinetic Energy**: `-ℏ²/(2mᵢ) ∇ᵢ²` — governs particle momentum dynamics
2. **External Potential**: `qᵢ φ(rᵢ, t)` — interaction with applied electric field
3. **Magnetic Coupling**: `-qᵢ/c (vᵢ · A(rᵢ, t))` — interaction with magnetic field
4. **Coulomb Interaction**: `qᵢ qⱼ / (4πε₀|rᵢ - rⱼ|)` — particle-particle electrostatic repulsion/attraction

## Features

- ✅ **Multi-particle systems**: Arbitrary number of particles with independent masses and charges
- ✅ **Full Coulomb interactions**: All pair-wise electrostatic interactions included
- ✅ **Time-dependent fields**: Support for arbitrary time-dependent EM field configurations
- ✅ **Multiple integrators**:
  - Split-step Fourier method (good for linear kinetic energy term)
  - Runge-Kutta 4th order (accurate for smooth potentials)
- ✅ **Observable computation**: Energy, probability density, expectation values
- ✅ **Physical constants**: Built-in fundamental constants (ℏ, ε₀, c, etc.)

## Physical Constants

The library uses SI units throughout:

- ℏ = 1.054571817 × 10⁻³⁴ J·s
- ε₀ = 8.854187812 × 10⁻¹² F/m
- c = 299,792,458 m/s
- e = 1.602176634 × 10⁻¹⁹ C
- mₑ = 9.109383701 × 10⁻³¹ kg
- mₚ = 1.672621900 × 10⁻²⁷ kg

## Core Types

### Particle
Represents a single quantum particle:
```haskell
data Particle = Particle
  { particleId   :: Int              -- Unique ID
  , mass         :: Double           -- Mass (kg)
  , charge       :: Double           -- Charge (C)
  , position     :: Vector Double    -- Position (m)
  , momentum     :: Vector Double    -- Momentum (kg·m/s)
  }
```

### System
Describes the complete quantum system:
```haskell
data System = System
  { particles    :: [Particle]       -- List of particles
  , gridPoints   :: Int              -- Grid resolution per dimension
  , boxSize      :: Double           -- Simulation box size (m)
  , currentTime  :: Double           -- Current simulation time (s)
  }
```

### EMField
Time-dependent electromagnetic field specification:
```haskell
type EMField = (Double, Vector Double, Vector Double) -> (Double, Vector Double)
-- Takes (time, position, velocity) and returns (φ, A)
-- where φ is scalar potential, A is vector potential
```

### Wavefunction
The quantum state:
```haskell
type Wavefunction = Matrix (Complex Double)
```

## Examples

### 1. Two-Electron System (Helium-like ion)

```haskell
let e = 1.602176634e-19  -- Elementary charge
    me = 9.1093837015e-31  -- Electron mass
    
    electron1 = Particle 1 me (-e) (fromList [-1e-10, 0, 0]) (fromList [0, 1e-24, 0])
    electron2 = Particle 2 me (-e) (fromList [1e-10, 0, 0]) (fromList [0, -1e-24, 0])
    
    system = System [electron1, electron2] 64 1e-9 0

-- Static electric field: E = 1e8 V/m
let eField (t, r, _) = (1e8 * (r ! 2), fromList [0, 0, 0])

-- Evolve for 10 fs
let (finalSys, finalPsi) = evolveSystem system initialPsi eField 1e-17 1e-14
```

### 2. Hydrogen Atom in Laser Field

```haskell
let electron = Particle 1 me (-e) (fromList [0, 0, 5.29e-11]) (fromList [0, 0, 0])
    system = System [electron] 128 5e-10 0

-- Laser field: linearly polarized, 200 nm wavelength, 10 MV/cm
let laserField (t, r, _) =
      let omega = 2 * pi * 5e15
          e_amp = 1e10
      in (-(e_amp / omega) * sin (omega * t) * (r ! 2),
          fromList [0, 0, e_amp / omega * cos (omega * t)])

let (finalSys, finalPsi) = evolveSystem system psi laserField 1e-18 1e-13
```

### 3. Coulomb Explosion

```haskell
let mp = 1.67262192e-27  -- Proton mass
    
    proton1 = Particle 1 mp e (fromList [-1e-15, 0, 0]) (fromList [0, 0, 0])
    proton2 = Particle 2 mp e (fromList [1e-15, 0, 0]) (fromList [0, 0, 0])
    
    system = System [proton1, proton2] 32 1e-13 0

-- No external field - pure Coulomb repulsion
let noField _ = (0, fromList [0, 0, 0])

let (finalSys, finalPsi) = evolveSystem system psi noField 1e-20 1e-15
```

## Usage

### Installation

```bash
cabal build
cabal install
```

### Running Examples

```bash
cabal run schrodinger-example
```

### Using as a Library

Add to your `cabal` file:
```
build-depends: schrodinger-solver
```

Then import:
```haskell
import Physics.Schrodinger
import Physics.Examples
```

## API Reference

### Time Evolution

```haskell
-- Evolve system to target time
evolveSystem :: System -> Wavefunction -> EMField -> TimeStep -> Double -> (System, Wavefunction)

-- Generate N evolution steps
evolutionSteps :: System -> Wavefunction -> EMField -> TimeStep -> Int -> [(System, Wavefunction)]
```

### Observable Computation

```haskell
-- Probability density |ψ(r, t)|²
probability :: Wavefunction -> Matrix Double

-- Expectation value ⟨O⟩ = ⟨ψ|O|ψ⟩
expectationValue :: Matrix (Complex Double) -> Wavefunction -> Complex Double

-- Total kinetic energy
kineticEnergy :: System -> Wavefunction -> Double

-- Total potential energy (Coulomb + external)
potentialEnergy :: System -> EMField -> Wavefunction -> Double
```

## Implementation Notes

### Numerical Methods

1. **Split-Step Fourier Method** (`splitStepEvolution`):
   - Alternates between position and momentum space
   - Separates kinetic and potential evolution
   - Excellent for periodic boundary conditions
   - Error: O(dt³)

2. **Runge-Kutta 4th Order** (`rk4Evolution`):
   - Direct integration in position space
   - Better for smooth, slowly-varying potentials
   - Handles boundary conditions naturally
   - Error: O(dt⁵)

### Grid Representation

The wavefunction is represented on an N-dimensional grid with:
- `gridPoints` per dimension (typically 64-256)
- Grid spacing: `dx = boxSize / gridPoints`
- Spatial domain: [-boxSize/2, boxSize/2]³

### Computational Complexity

- **Time evolution**: O(N · M³) per step (where N = particles, M = gridPoints)
- **FFT overhead**: O(M³ log M) with split-step method
- **Coulomb interaction**: O(N²) per timestep
- **Memory**: O(M³) for 3D wavefunction on M³ grid

## Physical Applications

1. **Atomic Ionization**: Study electron escape in strong fields
2. **Rydberg Dynamics**: High-energy state evolution in external fields
3. **Molecular Dynamics**: Track electron distribution in molecules
4. **Quantum Computing**: Model qubit systems with environmental interactions
5. **Plasma Physics**: Collective dynamics of charged particles
6. **Laser-Matter Interaction**: Strong-field physics and high-harmonic generation

## Limitations and Future Work

- [ ] Spin effects (add spin degree of freedom)
- [ ] Relativistic corrections (Dirac equation)
- [ ] Dissipation and decoherence (open quantum systems)
- [ ] GPU acceleration for large grid sizes
- [ ] Adaptive mesh refinement
- [ ] Better handling of singularities in Coulomb potential
- [ ] Periodic boundary conditions
- [ ] Absorbing boundary conditions (prevent reflection)

## References

1. Tannor, D. J. (2007). *Introduction to Quantum Mechanics: A Time-Dependent Perspective*
2. Thaller, B. (1992). *The Dirac Equation* (Springer)
3. Suzuki, M. (1976). "Decomposition formulas of exponential operators and Lie exponentials"
4. Kosloff, R. (1988). "Time-dependent quantum-mechanical methods for molecular dynamics"
5. Feit, M. D., Fleck, J. A., & Steiger, A. (1982). "Solution of the Schrödinger equation by a spectral method"

## License

MIT License - See LICENSE file for details

## Author

Quantum Computing Group
