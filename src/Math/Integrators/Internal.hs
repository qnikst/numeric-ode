module Math.Integrators.Internal
    where

{- | Integrator function
 -   \( \Phi [h] : y_0 \rightarrow y_1\)
 -}
type Integrator a = Double  -- ^ Step
                  -> a      -- ^ Initial value
                  -> a      -- ^ Next value

