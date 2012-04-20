module Math.Integrators.RK.Internal
    ( MExp(..)
    )
    where


-- | Internal type that users by
data MExp = Delimeter | Row (Maybe Double,[Double]) deriving (Show,Eq)

isExplicit :: [MExp] -> Bool
isExplicit  = (all (\(i,(Row (_,x))) -> i> length x)) . (zip [1..]) . top
    where 
        top = takeWhile (/= Delimeter) 

-- | helpers to work with data
rowRhs (Row (_,ls)) = ls
rowRhs _ = error "no a row"
