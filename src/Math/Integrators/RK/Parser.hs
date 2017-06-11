{-# OPTIONS_GHC -Wwarn #-} -- We need this option, because we want to remove this module in future
module Math.Integrators.RK.Parser
    ( readMatrixTable
    )
    where

import Data.Maybe

-- Parsec stuff
import Text.Parsec
import qualified Text.Parsec.Token as P
import Text.Parsec.Language (haskellDef) 
import Text.Parsec.Expr
import Text.Parsec.String

import Math.Integrators.RK.Internal

readMatrixTable :: String -> [MExp]
readMatrixTable  = 
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

lexer     = P.makeTokenParser haskellDef { P.reservedOpNames = ["*","/","+","-","sqrt","sin","cos"] }

whiteSpace= P.whiteSpace lexer
lexeme    = P.lexeme lexer
symbol    = P.symbol lexer
float     = P.float lexer
parens    = P.parens lexer
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

