{-# LANGUAGE FlexibleContexts #-}
-- |
-- Module: Math.Integrators.StormerVerlet
--
-- Stormer-Verlet has following properties
--
--  * the method is of order 2
--  * this is sympletic method
--  * method is symmetric
--  * this method excatly conserves quadratic first integrals \(p^T C\) (angular momentum of \(N\)-body system)
--
module Math.Integrators.StormerVerlet
    ( oneStepH98
    , stormerVerlet2
    ) where

import Linear
import Control.Lens

-- | Stormer-Verlet integration scheme for systems of the form
oneStepH98 :: (Applicative f, Num (f a), Fractional a) =>
              a            -- ^ Step size
           -> (f a -> f a) -- ^ \(\frac{\partial H}{\partial q}\)
           -> (f a -> f a) -- ^ \(\frac{\partial H}{\partial p}\)
           -> V2 (f a)     -- ^ Current \((p, q)\) as a 2-dimensional vector
           -> V2 (f a)     -- ^ New \((p, q)\) as a 2-dimensional vector
oneStepH98 hh nablaQ nablaP prev = V2 qNew pNew
  where
    h2   = hh / 2
    hhs  = pure hh
    hh2s = pure h2
    qsPrev = prev ^. _x
    psPrev = prev ^. _y
    pp2  = psPrev - hh2s * nablaQ qsPrev
    qNew = qsPrev + hhs * nablaP pp2
    pNew = pp2 - hh2s * nablaQ qNew

-- | Stormer-Verlet integration scheme for system: 
--      
--   \[   q'' = f(q) \]
--
-- <<diagrams/src_Math_Integrators_StormerVerlet_mySquare.svg#diagram=mySquare&height=300&width=200>>
--
-- > mySquare = square 1 # fc blue # myTransf
-- > myTransf = rotateBy (1/27)
--
-- For example, for the \(n\)-body problem, the Hamiltonian is
--
-- \[
-- {\mathbb H} = \frac{1}{2}\sum_{i=0}^n \frac{p_i^\top p_i}{m_i} - \frac{G}{2}\sum_{i=0}^n\sum_{j \neq i} \frac{m_i m_j}{\|q_i - q_j\|}
-- \]
--
-- Apply Hamilton's equations will give us \(2n\) first order
-- equations. To use 'stormerVerlet2' we need \(n\) second order
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
-- m_1\ddot{q}_1 &= G\frac{q_1 - q_2}{\|q_1 - q_2\|^3} \\
-- m_2\ddot{q}_2 &= G\frac{q_2 - q_1}{\|q_2 - q_1\|^3}
-- \end{aligned}
-- \]
--
-- > kepler :: L.V2 (L.V3 Double) -> L.V2 (L.V3 Double)
-- > kepler (L.V2 q1 q2) =
-- >     let r  = q2 L.^-^ q1          -- q2 - q1
-- >         ri = r `L.dot` r          -- ||q2-q1||^2
-- >         rr = ri * (sqrt ri)
-- >         q1' = r / pure rr
-- >         q2' = negate q1'
-- >     in L.V2 q1' q2'
--
stormerVerlet2 :: (Applicative f, Num (f a), Fractional a)
               => (f a -> f a)              -- ^ function f
               -> a
               -> V2 (f a)
               -> V2 (f a)
stormerVerlet2 f h prev =
    let h'  = h
        h2' = 0.5 *h
        p1  = prev ^. _x + pure h2' * (f (prev ^. _y))
        q'  = prev ^. _y + pure h' * p1
        p'  = p1 + pure h2' * (f q')
    in V2 p' q'

