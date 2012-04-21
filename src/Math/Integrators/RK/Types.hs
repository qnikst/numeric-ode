module Math.Integrators.RK.Types
    where


-- | type implicit solver
data ImplicitRkType a = FixedPoint (Int -> a -> a -> Bool) | NewtonIteration
