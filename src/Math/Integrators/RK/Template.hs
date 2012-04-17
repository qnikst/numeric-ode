{-# LANGUAGE TemplateHaskell #-}
module Math.Integrators.RK.Template
    where

import Data.Maybe

import Language.Haskell.TH.Quote
import Language.Haskell.TH


import Text.Parsec
import qualified Text.Parsec.Token as P
import Text.Parsec.Language (haskellDef) 
import Text.Parsec.Expr
import Text.Parsec.String
import qualified GHC.Num

import Control.Monad

import Debug.Trace

lexer       = P.makeTokenParser haskellDef { P.reservedOpNames = ["*","/","+","-","sqrt","sin","cos"] }

whiteSpace= P.whiteSpace lexer
lexeme    = P.lexeme lexer
symbol    = P.symbol lexer
float     = P.float lexer
parens    = P.parens lexer
semi      = P.semi lexer
natural   = P.natural lexer
identifier= P.identifier lexer
reserved  = P.reserved lexer
reservedOp= P.reservedOp lexer

expr    :: Parser Double
expr    = buildExpressionParser table factor
        <?> "expression"
factor  =   parens expr
        <|> try float
        <|> fmap realToFrac natural
        <?> "simple expression"
table   = [ [prefix "-" negate]
          , [prefix "sqrt" sqrt,prefix "sin" sin,prefix "cos" cos]
          , [op "*" (*) AssocLeft, op "/" (/) AssocLeft]
          , [op "+" (+) AssocLeft, op "-" (-) AssocLeft]
          ]          
          where
             op s f assoc = Infix (do{ reservedOp s; return f} <?> "operator") assoc
             prefix s f   = Prefix (do { reservedOp s; return f} <?> "prefix")


erun :: Parser Double -> String -> Double
erun p input = erun' (do { whiteSpace ; x <- p ; eof; return x})
    where
        erun' p' = case (parse p' "" input) of
                    Left err -> error $ "Parse error at "++(show err)
                    Right x  -> x

data MExp = Delimeter | Row (Maybe Double,[Double]) deriving (Show,Eq)

rowRhs (Row (_,ls)) = ls
rowRhs _ = error "no a row"


valuesFromString :: String -> [MExp]
valuesFromString = 
    mapMaybe go . lines
    where
        go ('-':s) = Just Delimeter
        go ('#':s) = Nothing
        go (ls) | null.filter (==' ')$ ls = Nothing
                | otherwise = 
            let (lhs,_:rhs) = span (/='|') ls
                l = case filter (/= ' ') lhs of
                        "" -> Nothing
                        ls -> Just $! erun expr ls
            in Just $ Row (l, map (erun expr) $! grp '&' rhs)
        grp c s = case dropWhile (==c) s of
            "" -> []
            s' -> if any (/=' ') w then w : grp c c'' else grp c c''
                  where (w,c'') = break (==c) s'


qrk  :: QuasiQuoter
qrk  =  QuasiQuoter {quoteExp = x}
    where
        x s = case () of
                _ | trace (show $ valuesFromString s) False -> undefined
                _ -> rk $ valuesFromString s

jv = Just . VarE
jld = Just . LitE. {-DoublePrimL-} RationalL . toRational
ld    = LitE . RationalL . toRational
plus  = VarE $! mkName "+"
vplus = VarE $! mkName "^+^"
vmult = VarE $! mkName "*^"
mult  = VarE $! mkName "*"



rk :: [MExp] -> Q Exp
rk mExp = do
    let (ab,_:c:[]) = break (==Delimeter) mExp
        lenA     = length ab
    let f = mkName "f"
        h = mkName "h"
        t = mkName "t"
        y = mkName "y"
    kn <- forM [1..lenA] (\_ -> newName "k")
    let kvv = zip kn ab
    ks' <- forM kvv $ \(k,r) -> do
            t <- rkt1 r kn
            return $ ValD (VarP k) (NormalB t) []
    y' <- rkt2 c kn
    return $ LamE [VarP f, VarP h, TupP [VarP t,VarP y]] $ LetE ks' y'
    where 
        rkt1 (Row (Just c,ls)) ks = do
            let f = mkName "f"
                h = mkName "h"
                t = mkName "t"
                y = mkName "y"
            let ft = InfixE (jv t) plus (Just $ InfixE (jld c) mult (jv h))
                st = if null ls 
                        then (VarE y)
                        else
                            InfixE (jv y) vplus
                                (Just $ InfixE (Just $ AppE (VarE $ mkName "realToFrac") (VarE h)) (vmult)
                                    (Just $ foldl1 (\x y -> InfixE (Just x) vplus (Just y)) $
                                                zipWith (\k l -> InfixE (Just $ AppE (VarE $ mkName "realToFrac") (ld l)) vmult (jv k)) ks ls
                                    )
                                )
            return $ AppE ( AppE (VarE f) ft) st
        rkt2 (Row (_,ls)) ks = do
            let f = mkName "f"
                y = mkName "y"
                h = mkName "h"
                t = mkName "t"
            return $ TupE 
                        [ InfixE (jv t) plus (jv h)
                        , InfixE (jv y) vplus
                          (Just $ InfixE (Just $ AppE (VarE $ mkName "realToFrac") (VarE h)) vmult
                            (Just $ foldl1 (\x y -> InfixE (Just x) vplus (Just y)) $
                                        zipWith (\k l -> InfixE (Just $ AppE (VarE $ mkName "realToFrac") (ld l)) vmult (jv k)) ks ls
                            )
                          )
                        ]

{-
implicit = \f implicit (t,y) -> 
    let (z1,z2,...zN) = implicit
                            (\(z1,z2..zN) ->
                                let z1' = if (z1 > eps) then h * sum (a_{ij} f (y+z1)) else z1
                                    z2' = if (z2 > eps) then h * sum (................) else z2
                                    ....
                            in tuple (z1',z2',z3'..zN'))
                            (z1,z2,z3,...,zN')
    in y ^+^ h *^ sum ( b_1 * (y ^+^ z1))
-}

