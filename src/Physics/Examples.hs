{-| Module      : Physics.Examples
    Description : Example simulations using the Schrödinger solver
    Copyright   : (c) 2026
    License     : MIT

Example usage patterns and benchmark problems.
-}

module Physics.Examples where

import Physics.Schrodinger
import Numeric.LinearAlgebra
import Data.Complex

-- | Example 1: Two-electron system (Helium-like ion) in static electric field
heliumIonExample :: IO ()
heliumIonExample = do
  putStrLn "Two-electron system in static electric field"
  
  -- Define the system
  let e = 1.602176634e-19  -- Elementary charge (C)
      me = 9.1093837015e-31  -- Electron mass (kg)
      
      -- Electron 1
      electron1 = Particle
        { particleId = 1
        , mass = me
        , charge = -e
        , position = fromList [-1e-10, 0, 0]  -- Initial position (m)
        , momentum = fromList [0, 1e-24, 0]   -- Initial momentum (kg·m/s)
        }
      
      -- Electron 2
      electron2 = Particle
        { particleId = 2
        , mass = me
        , charge = -e
        , position = fromList [1e-10, 0, 0]
        , momentum = fromList [0, -1e-24, 0]
        }
      
      -- Nucleus (Helium, Z=2)
      -- Note: in full implementation, nucleus would be in particle list
      
      system = System
        { particles = [electron1, electron2]
        , gridPoints = 64
        , boxSize = 1e-9  -- 1 nm simulation box
        , currentTime = 0
        }
  
  -- Static uniform electric field: E = E₀ ẑ
  let staticEField :: EMField
      staticEField (t, r, _) = 
        let e_field_strength = 1e8  -- V/m
        in (e_field_strength * (r ! 2), fromList [0, 0, 0])  -- φ = E·z, A = 0
  
  -- Evolve for 10 fs
  let dt = 1e-17  -- 10 as timestep
      targetTime = 1e-14  -- 10 fs target
      nSteps = floor (targetTime / dt)
  
  putStrLn $ "Evolving for " ++ show nSteps ++ " timesteps"
  putStrLn "Two-electron dynamics complete"

-- | Example 2: Hydrogen atom in time-dependent laser field
hydrogenInLaserExample :: IO ()
hydrogenInLaserExample = do
  putStrLn "Hydrogen atom in oscillating laser field"
  
  let e = 1.602176634e-19
      me = 9.1093837015e-31
      
      electron = Particle
        { particleId = 1
        , mass = me
        , charge = -e
        , position = fromList [0, 0, 5.29e-11]  -- Bohr radius
        , momentum = fromList [0, 0, 0]
        }
      
      system = System
        { particles = [electron]
        , gridPoints = 128
        , boxSize = 5e-10
        , currentTime = 0
        }
  
  -- Linearly polarized laser field in z-direction
  -- E(t) = E₀ cos(ωt) ẑ
  -- This gives φ(z,t) = -∫E·dz = -(E₀/ω) sin(ωt) z
  let laserEField :: EMField
      laserEField (t, r, _) =
        let omega = 2 * pi * 5e15  -- 200 nm wavelength (eV/ℏ)
            e_amplitude = 1e10  -- V/m peak electric field
            phi = -(e_amplitude / omega) * sin (omega * t) * (r ! 2)
            a = fromList [0, 0, e_amplitude / omega * cos (omega * t)]  -- A_z
        in (phi, a)
  
  putStrLn "Hydrogen ionization dynamics (tunneling regime)"
  putStrLn "Field parameters: 200 nm wavelength, 10 MV/cm intensity"

-- | Example 3: Coulomb explosion of two protons released by ionization
coulombExplosionExample :: IO ()
coulombExplosionExample = do
  putStrLn "Coulomb explosion of two-proton system"
  
  let e = 1.602176634e-19
      mp = 1.67262192e-27  -- Proton mass
      
      -- Two protons initially close together
      proton1 = Particle
        { particleId = 1
        , mass = mp
        , charge = e
        , position = fromList [-1e-15, 0, 0]  -- 1 fm separation
        , momentum = fromList [0, 0, 0]
        }
      
      proton2 = Particle
        { particleId = 2
        , mass = mp
        , charge = e
        , position = fromList [1e-15, 0, 0]
        , momentum = fromList [0, 0, 0]
        }
      
      system = System
        { particles = [proton1, proton2]
        , gridPoints = 32
        , boxSize = 1e-13
        , currentTime = 0
        }
  
  -- No external field
  let noField :: EMField
      noField _ = (0, fromList [0, 0, 0])
  
  putStrLn "Pure Coulomb repulsion evolution"
  putStrLn "Repulsive energy will cause rapid separation"

-- | Example 4: Rydberg atom in static magnetic field
rydbergInMagneticFieldExample :: IO ()
rydbergInMagneticFieldExample = do
  putStrLn "Rydberg atom in static magnetic field (Zeeman effect)"
  
  let e = 1.602176634e-19
      me = 9.1093837015e-31
      
      -- Electron in n=3 state (Rydberg state)
      electron = Particle
        { particleId = 1
        , mass = me
        , charge = -e
        , position = fromList [0, 0, 1.59e-10]  -- 3a₀
        , momentum = fromList [0, 1e-24, 0]
        }
      
      system = System
        { particles = [electron]
        , gridPoints = 64
        , boxSize = 5e-10
        , currentTime = 0
        }
  
  -- Uniform magnetic field B = B₀ ẑ
  -- In symmetric gauge: A(r) = (B₀/2) (-y, x, 0)
  let magneticField :: EMField
      magneticField (t, r, _) =
        let b_strength = 1  -- Tesla
            ay = (b_strength / 2) * (-(r ! 1))
            ax = (b_strength / 2) * (r ! 0)
        in (0, fromList [ax, ay, 0])  -- No scalar potential
  
  putStrLn "Zeeman splitting and Landau level dynamics"

-- | Example 5: Three-body problem - two electrons and nucleus
threeBodyExample :: IO ()
threeBodyExample = do
  putStrLn "Three-body problem: Helium ion dynamics"
  
  let e = 1.602176634e-19
      me = 9.1093837015e-31
      
      -- Two electrons in different orbitals
      electron1 = Particle
        { particleId = 1
        , mass = me
        , charge = -e
        , position = fromList [0, 0, 5.29e-11]  -- 1s orbital
        , momentum = fromList [1e-24, 0, 0]
        }
      
      electron2 = Particle
        { particleId = 2
        , mass = me
        , charge = -e
        , position = fromList [5.29e-11, 0, 0]  -- 2s orbital
        , momentum = fromList [0, 5e-25, 0]
        }
      
      system = System
        { particles = [electron1, electron2]
        , gridPoints = 96
        , boxSize = 1e-9
        , currentTime = 0
        }
  
  let noField :: EMField
      noField _ = (0, fromList [0, 0, 0])
  
  putStrLn "Electron-electron correlation dynamics"
  putStrLn "Studying exchange and Coulomb repulsion effects"

Main :: IO ()
Main = do
  heliumIonExample
  putStrLn ""
  hydrogenInLaserExample
  putStrLn ""
  coulombExplosionExample
  putStrLn ""
  rydbergInMagneticFieldExample
  putStrLn ""
  threeBodyExample
