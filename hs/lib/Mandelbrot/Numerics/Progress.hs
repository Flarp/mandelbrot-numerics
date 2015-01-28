module Mandelbrot.Numerics.Progress
  ( Progress(..)
  , skip
  ) where

data Progress a = Failed !a | Continue !a (Progress a) | Done !a
  deriving (Eq, Ord, Read, Show)

skip :: Int -> Progress a -> Progress a
skip n p
  | n <= 0 = p
  | otherwise = case p of
      Continue _ p' -> skip (n - 1) p'
      _ -> p
