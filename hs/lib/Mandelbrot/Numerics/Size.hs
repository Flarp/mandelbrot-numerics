{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE FlexibleContexts #-}
module Mandelbrot.Numerics.Size
  ( size
  ) where

import Mandelbrot.Numerics.Complex

size
  :: (RealFloat r, Square r, Square (Complex r))
  => Int -> Complex r -> Complex r
size p c = go 1 0 1 1
  where
    go q !z !b !l
      | q == p = recip (b * sqr l)
      | otherwise = go (q + 1) z' b' l'
      where
        z' = sqr z + c
        l' = double (z' * l)
        b' = b + recip l'
