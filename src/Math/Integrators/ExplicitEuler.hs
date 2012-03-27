module Math.Integrators.ExplicitEuler

explicitEuler :: Double -> (a -> a) -> a -> a
explicitEuler h f y = y+h * (f y)
