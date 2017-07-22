{-# LANGUAGE FlexibleContexts #-}

-- |
-- Module: Math.Integrators.StormerVerlet
--
--
-- Störmer-Verlet is an order 2 symplectic method. This means it will
-- preserve the Hamiltonian for the system the differential equations
-- describe, for example, important for modelling planetary motion;
-- the application of something like the much-loved Runge-Kutta 4th
-- order method would either model the planet spiralling toward or
-- away from the Sun!
--
-- Here's a diagram showing the orbit of Jupiter around the Sun.
-- 
-- <<diagrams/src_Math_Integrators_StormerVerlet_jupiterOrbit.svg#diagram=jupiterOrbit&height=400&width=500>>
--
-- To create this, consider the \(n\)-body problem. The Hamiltonian is
--
-- \[
-- {\mathbb H} = \frac{1}{2}\sum_{i=0}^n \frac{p_i^\top p_i}{m_i} - \frac{G}{2}\sum_{i=0}^n\sum_{j \neq i} \frac{m_i m_j}{\|q_i - q_j\|}
-- \]
--
-- Apply Hamilton's equations will gives \(2n\) first order
-- equations. To use 'stormerVerlet2' this needs to be \(n\) second order
-- equations. In this case, the Lagrangian is easy
--
-- \[
-- {\mathcal{L}} = \frac{1}{2}\sum_{i=0}^n \frac{p_i^\top p_i}{m_i} + \frac{G}{2}\sum_{i=0}^n\sum_{j \neq i} \frac{m_i m_j}{\|q_i - q_j\|}
-- \]
--
-- Applying Lagrange's equation
--
-- \[
-- \frac{\mathrm{d}}{\mathrm{d}t}\bigg(\frac{\partial{\mathcal{L}}}{\partial\dot{q}_j}\bigg) = \frac{\partial{\mathcal{L}}}{\partial{q}_j}
-- \]
--
-- gives
--
-- \[
-- m_j\ddot{q}_j = G\sum_{k \neq j}m_k m_j \frac{q_k - q_j}{\|q_k - q_j\|^3}
-- \]
--
-- For \(n = 2\) this gives
--
-- \[
-- \begin{aligned}
-- \ddot{q}_1 &= m_2G\frac{q_1 - q_2}{\|q_1 - q_2\|^3} \\
-- \ddot{q}_2 &= m_1G\frac{q_2 - q_1}{\|q_2 - q_1\|^3}
-- \end{aligned}
-- \]
--
-- > {-# LANGUAGE NegativeLiterals      #-}
-- > {-# LANGUAGE TypeFamilies          #-}
-- > {-# LANGUAGE FlexibleContexts      #-}
-- > {-# LANGUAGE MultiParamTypeClasses #-}
-- > 
-- > import qualified Data.Vector as V
-- > import Control.Monad.ST
-- > 
-- > import Math.Integrators.StormerVerlet
-- > import Math.Integrators
-- > 
-- > import qualified Linear as L
-- > import Linear.V
-- > import Data.Maybe ( fromJust )
-- > 
-- > import Diagrams.Prelude
-- > 
-- > import Control.Monad
-- > import Control.Monad.State.Class
-- > 
-- > import Plots
-- >
-- > -- First some constants describing the system
-- >
-- > gConst :: Double
-- > gConst = 6.67384e-11
-- > 
-- > nStepsTwoPlanets :: Int
-- > nStepsTwoPlanets = 44
-- >
-- > -- A step size of 100 days!
-- >
-- > stepTwoPlanets :: Double
-- > stepTwoPlanets = 24 * 60 * 60 * 100
-- > 
-- > sunMass, jupiterMass :: Double
-- > sunMass     = 1.9889e30
-- > jupiterMass = 1.8986e27
-- > 
-- > jupiterPerihelion :: Double
-- > jupiterPerihelion = 7.405736e11
-- > 
-- > jupiterV :: L.V3 Double
-- > jupiterV = L.V3 (-1.0965244901087316e02) (-1.3710001990210707e04) 0.0
-- > 
-- > jupiterQ :: L.V3 Double
-- > jupiterQ = L.V3 (-jupiterPerihelion) 0.0 0.0
-- > 
-- > sunV :: L.V3 Double
-- > sunV = L.V3 0.0 0.0 0.0
-- > 
-- > sunQ :: L.V3 Double
-- > sunQ =  L.V3 0.0 0.0 0.0
-- >
-- > -- The right hand side of the second order differential equation system.
-- >
-- > kepler :: L.V2 (L.V3 Double) -> L.V2 (L.V3 Double)
-- > kepler (L.V2 q1 q2) =
-- >     let r  = q2 L.^-^ q1
-- >         ri = r `L.dot` r
-- >         rr = ri * (sqrt ri)
-- >         q1' = pure gConst * r / pure rr
-- >         q2' = negate q1'
-- >         q1'' = q1' * pure sunMass
-- >         q2'' = q2' * pure jupiterMass
-- >     in L.V2 q1'' q2''
-- >
-- > -- Initial values
-- >
-- > initPQs :: L.V2 (L.V2 (L.V3 Double))
-- > initPQs = L.V2 (L.V2 jupiterV sunV) (L.V2 jupiterQ sunQ)
-- >
-- > -- Steps at which to evolve the system
-- >
-- > tm :: V.Vector Double
-- > tm = V.enumFromStepN 0 stepTwoPlanets nStepsTwoPlanets
-- >
-- > -- The results
-- >
-- > result1 :: V.Vector (L.V2 (L.V2 (L.V3 Double)))
-- > result1 = runST $ integrateV (\h -> stormerVerlet2 kepler (pure h)) initPQs tm
-- > 
-- > preMorePts :: [(Double, Double)]
-- > preMorePts = map (\(L.V2 _ (L.V2 (L.V3 x y _z) _)) -> (x,y))  (V.toList result1)
-- > 
-- > morePts :: [P2 Double]
-- > morePts = map p2 $ preMorePts
-- >
-- > -- Finally plot the results
-- >
-- > addPoint :: (Plotable (Diagram B) b, MonadState (Axis b V2 Double) m) =>
-- >             Double -> (Double, Double) -> m ()
-- > addPoint o (x, y) = addPlotable'
-- >                     ((circle 1e11 :: Diagram B) #
-- >                      fc brown #
-- >                      opacity o #
-- >                      translate (r2 (x, y)))
-- > 
-- > jSaxis :: Axis B V2 Double
-- > jSaxis = r2Axis &~ do
-- >   addPlotable' ((circle 1e11 :: Diagram B) # fc yellow)
-- >   let l = length preMorePts
-- >   let os = [0.05,0.1..]
-- >   let ps = take (l `div` 4) [0,4..]
-- >   zipWithM_ addPoint os (map (preMorePts!!) ps)
-- >   linePlot' $ map unp2 $ take 200 morePts
-- > 
-- > jupiterOrbit = renderAxis jSaxis # bg ivory
--

module Math.Integrators.StormerVerlet
    ( stormerVerlet2H
    , stormerVerlet2
    ) where

import Linear
import Control.Lens

-- | Störmer-Verlet integration scheme for systems of the form
-- \(\mathbb{H}(p,q) = T(p,q) + V(p,q)\)
stormerVerlet2H :: (Applicative f, Num (f a), Fractional a) =>
              a            -- ^ Step size
           -> (f a -> f a) -- ^ \(\frac{\partial H}{\partial q}\)
           -> (f a -> f a) -- ^ \(\frac{\partial H}{\partial p}\)
           -> V2 (f a)     -- ^ Current \((p, q)\) as a 2-dimensional vector
           -> V2 (f a)     -- ^ New \((p, q)\) as a 2-dimensional vector
stormerVerlet2H hh nablaQ nablaP prev = V2 qNew pNew
  where
    h2   = hh / 2
    hhs  = pure hh
    hh2s = pure h2
    qsPrev = prev ^. _x
    psPrev = prev ^. _y
    pp2  = psPrev - hh2s * nablaQ qsPrev
    qNew = qsPrev + hhs * nablaP pp2
    pNew = pp2 - hh2s * nablaQ qNew

-- | Störmer-Verlet integration scheme for system: \(\ddot{\mathbf{q}} = f(\mathbf{q})\)
stormerVerlet2 :: (Applicative f, Num (f a), Fractional a)
               => (f a -> f a)              -- ^ \(f\)
               -> a                         -- ^ Step size
               -> V2 (f a)                  -- ^ Current \((p, q)\) as a 2-dimensional vector
               -> V2 (f a)                  -- ^ New \((p, q)\) as a 2-dimensional vector
stormerVerlet2 f h prev =
    let h'  = h
        h2' = 0.5 * h
        p1  = prev ^. _x + pure h2' * (f (prev ^. _y))
        q'  = prev ^. _y + pure h' * p1
        p'  = p1 + pure h2' * (f q')
    in V2 p' q'

