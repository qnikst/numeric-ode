module Math.Integrators.RK.Types
    where


-- | type implicit solver
data ImplicitRkType = FixedPoint | NewtonIteration
                    deriving (Eq,Show)


