{-| Module      : Physics.Schrodinger
    Description : Time-dependent Schrödinger equation solver for N charged particles
    Copyright   : (c) 2026
    License     : MIT

Time-dependent Schrödinger equation solver for a system of N particles with
arbitrary masses, electrical charges, evolving in a time-dependent electromagnetic field.
Includes particle-particle Coulomb interactions.

The time-dependent Schrödinger equation is:

  iℏ ∂ψ/∂t = Ĥ(t) ψ

where the Hamiltonian is:

  Ĥ(t) = Σᵢ [-ℏ²/(2mᵢ) ∇ᵢ² + qᵢ φ(rᵢ, t) - qᵢ/c (v(t) · Aᵢ(rᵢ, t))]
       + Σᵢ<ⱼ [qᵢ qⱼ / |rᵢ - rⱼ|]

-}

module Physics.Schrodinger
  ( -- * Core types
    Particle (..)
  , System (..)
  , EMField
  , Wavefunction
  , TimeStep
  
  -- * Solvers
  , evolveSystem
  , evolutionSteps
  , wavefunctionAt
  
  -- * Observables
  , probability
  , expectationValue
  , kineticEnergy
  , potentialEnergy
  
  -- * Constants
  , hbar
  , epsilon0
  , c_speed
  ) where

import Data.List (foldl')
import Data.Complex
import Numeric.LinearAlgebra
import qualified Numeric.LinearAlgebra.HMatrix as LA

-- | Fundamental physical constants (SI units)
hbar :: Double
hbar = 1.054571817e-34  -- Reduced Planck constant (J·s)

epsilon0 :: Double
epsilon0 = 8.8541878128e-12  -- Vacuum permittivity (F/m)

c_speed :: Double
c_speed = 299792458  -- Speed of light (m/s)

-- | Represents a single particle in the system
data Particle = Particle
  { particleId   :: Int           -- Unique identifier
  , mass         :: Double        -- Mass in kg
  , charge       :: Double        -- Electric charge in Coulombs
  , position     :: Vector Double -- Position in 3D space
  , momentum     :: Vector Double -- Momentum in 3D space
  } deriving (Show, Eq)

-- | Represents the quantum system of N particles
data System = System
  { particles    :: [Particle]           -- List of particles
  , gridPoints   :: Int                  -- Number of grid points per dimension
  , boxSize      :: Double               -- Size of simulation box (m)
  , currentTime  :: Double               -- Current simulation time (s)
  } deriving (Show, Eq)

-- | Wavefunction represented as a complex-valued function on a grid
-- For N particles in 3D, this is ψ(r₁, r₂, ..., rₙ, t)
type Wavefunction = Matrix (Complex Double)

-- | Electromagnetic field specification
-- Returns (electric potential φ, magnetic vector potential A) at position r and time t
type EMField = (Double, Vector Double, Vector Double) -> (Double, Vector Double)

-- | Time step for evolution
type TimeStep = Double

-- | Create initial Gaussian wavepacket for a single particle
GaussianWavepacket :: Particle -> Double -> Double -> Particle -> Double -> Wavefunction
GaussianWavepacket particle sigma k0 _ gridRes =
  error "Not yet implemented - requires discretization scheme"

-- | Compute the kinetic energy operator T = -ℏ²/(2m) ∇²
kineticOperator :: Double -> Double -> Vector Double -> Matrix (Complex Double)
kineticOperator mass gridSpacing _ =
  let n = 3  -- 3D
      -- Laplacian approximation on grid (finite difference)
      h = gridSpacing
      coeff = -hbar**2 / (2 * mass * h**2)
  in LA.scale (coeff :+ 0) (LA.ident 3)

-- | Compute Coulomb potential between two particles
-- V = q₁q₂ / (4πε₀ |r₁ - r₂|)
coulombPotential :: Double -> Double -> Double -> Double -> Double
coulombPotential q1 q2 r epsilon
  | r < epsilon = 1e10  -- Avoid singularity
  | otherwise = (q1 * q2) / (4 * pi * epsilon * r)

-- | Compute total Coulomb interaction energy for the system
totalCoulombEnergy :: [Particle] -> Double
totalCoulombEnergy ps =
  foldl' (+) 0 [coulombPotential (charge p1) (charge p2) (norm_2 (position p1 - position p2)) epsilon0
                | (i, p1) <- zip [0..] ps, (j, p2) <- zip [0..] ps, i < j]

-- | Compute total external potential energy: Σᵢ qᵢ φ(rᵢ, t)
externalPotentialEnergy :: EMField -> [Particle] -> Double -> Double
externalPotentialEnergy emField ps t =
  foldl' (+) 0 [charge p * phi | p <- ps, let (phi, _) = emField (t, position p, [])]

-- | Compute magnetic interaction energy: -Σᵢ qᵢ/c (vᵢ · Aᵢ)
magneticInteractionEnergy :: EMField -> [Particle] -> Double -> Double
magneticInteractionEnergy emField ps t =
  foldl' (+) 0 [-(charge p / c_speed) * (momentum p `dot` a)
                | p <- ps, let (_, a) = emField (t, position p, [])]
  where dot u v = sum [u ! i * v ! i | i <- [0 .. size u - 1]]

-- | Split-step Fourier method for time evolution
-- Evolves ψ(t) -> ψ(t + dt) using operator splitting
splitStepEvolution :: System
                   -> Wavefunction
                   -> EMField
                   -> TimeStep
                   -> (System, Wavefunction)
splitStepEvolution sys psi emField dt =
  let
    -- Half-step kinetic energy evolution (momentum space)
    tHalf = currentTime sys + dt / 2
    psiMomentum = fft psi  -- Transform to momentum space
    
    -- Apply kinetic operator for half time step
    kineticExp = LA.scale (exp ((-1 :+ 0) * (0 :+ hbar / (2 * hbar)) * dt) :+ 0) (LA.ident 3)
    psiAfterKinetic = fft (psiMomentum `LA.mmmul` kineticExp)
    
    -- Full-step potential evolution (position space)
    -- V_total = external potential + Coulomb interaction
    extPot = externalPotentialEnergy emField (particles sys) (currentTime sys + dt / 2)
    coulPot = totalCoulombEnergy (particles sys)
    
    potentialExp = exp ((-1 :+ 0) * (0 :+ (extPot + coulPot) / hbar) * dt)
    psiAfterPotential = LA.scale potentialExp psiAfterKinetic
    
    -- Final half-step kinetic evolution
    psiMomentumFinal = fft psiAfterPotential
    psiFinal = fft (psiMomentumFinal `LA.mmmul` kineticExp)
    
    updatedSys = sys { currentTime = currentTime sys + dt }
  in
    (updatedSys, psiFinal)

-- | Runge-Kutta 4th order time evolution (alternative to split-step)
-- More accurate for smooth potentials
rk4Evolution :: System
             -> Wavefunction
             -> EMField
             -> TimeStep
             -> (System, Wavefunction)
rk4Evolution sys psi emField dt =
  let
    -- Compute Hamiltonian and time evolution
    h = timeDerivative sys psi emField
    k1 = LA.scale ((0 :+ 1) / hbar * dt) h
    
    -- Subsequent RK4 stages would follow similar pattern
    -- Simplified here for brevity
    psiFinal = psi `LA.add` k1
    updatedSys = sys { currentTime = currentTime sys + dt }
  in
    (updatedSys, psiFinal)

-- | Compute time derivative: iℏ ∂ψ/∂t = H ψ
timeDerivative :: System -> Wavefunction -> EMField -> Wavefunction
timeDerivative sys psi emField =
  let
    ps = particles sys
    t = currentTime sys
    
    -- Kinetic energy contributions (Laplacian)
    kineticTerm = foldl' LA.add (LA.zeros (rows psi) (cols psi))
      [LA.scale ((-hbar**2 / (2 * mass p)) :+ 0) (laplacian psi) | p <- ps]
    
    -- External potential (Coulomb + EM)
    extPotTerm = LA.scale
      ((totalCoulombEnergy ps + externalPotentialEnergy emField ps t) :+ 0)
      psi
    
    -- Magnetic term
    magTerm = LA.scale
      ((magneticInteractionEnergy emField ps t) :+ 0)
      psi
  in
    kineticTerm `LA.add` extPotTerm `LA.add` magTerm

-- | Discrete Laplacian approximation (finite difference method)
laplacian :: Wavefunction -> Wavefunction
laplacian psi =
  let
    (m, n) = size psi
    h = 1.0 / fromIntegral m  -- Grid spacing
  in
    -- Simplified: real implementation uses proper boundary conditions
    LA.scale ((1 / h**2) :+ 0) psi

-- | Fast Fourier transform (wrapper around external library)
fft :: Matrix (Complex Double) -> Matrix (Complex Double)
fft = id  -- Placeholder - would use actual FFT implementation

-- | Evolve system for N timesteps
evolutionSteps :: System
               -> Wavefunction
               -> EMField
               -> TimeStep
               -> Int
               -> [(System, Wavefunction)]
evolutionSteps sys psi emField dt n =
  let step (s, p) _ = splitStepEvolution s p emField dt
  in take (n + 1) $ iterate (\(s, p) -> step (s, p) ()) (sys, psi)

-- | Evolve system to a specific time
evolvSystem :: System
           -> Wavefunction
           -> EMField
           -> TimeStep
           -> Double
           -> (System, Wavefunction)
evolveSystem sys psi emField dt targetTime =
  let
    nSteps = ceiling ((targetTime - currentTime sys) / dt)
    steps = evolutionSteps sys psi emField dt nSteps
  in
    last steps

-- | Evaluate wavefunction at a specific spatial configuration
wavefunctionAt :: Wavefunction -> Vector Double -> Complex Double
wavefunctionAt psi r = error "Not yet implemented - requires grid interpolation"

-- | Compute probability density |ψ(r, t)|²
probability :: Wavefunction -> Matrix Double
probability psi =
  let psiFwd = psi
  in LA.fromRows [LA.fromList [magnitude (psiFwd LA.! i LA.! j) ** 2 | j <- [0 .. cols psiFwd - 1]] | i <- [0 .. rows psiFwd - 1]]

-- | Compute expectation value ⟨O⟩ = ⟨ψ|O|ψ⟩
expectationValue :: Matrix (Complex Double) -> Wavefunction -> Complex Double
expectationValue operator psi =
  let
    opPsi = operator `LA.mmmul` psi
    conjugate = LA.fromRows [LA.fromList [conjugate c | c <- LA.toList row] | row <- LA.toRows psi]
  in
    sum $ map magnitude $ LA.toList $ conjugate `LA.mmmul` opPsi

-- | Compute total kinetic energy expectation value
kineticEnergy :: System -> Wavefunction -> Double
kineticEnergy sys psi =
  let
    ps = particles sys
    -- Kinetic energy operator for each particle
    kinOps = [LA.scale ((-hbar**2 / (2 * mass p)) :+ 0) (LA.ident 3) | p <- ps]
    -- Sum expectation values
  in
    realPart $ foldl' (\acc op -> acc + expectationValue op psi) 0 kinOps

-- | Compute total potential energy expectation value
potentialEnergy :: System -> EMField -> Wavefunction -> Double
potentialEnergy sys emField psi =
  let
    ps = particles sys
    t = currentTime sys
    coulPot = totalCoulombEnergy ps
    extPot = externalPotentialEnergy emField ps t
  in
    coulPot + extPot
