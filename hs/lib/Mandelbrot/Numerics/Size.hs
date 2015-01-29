{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE FlexibleContexts #-}
module Mandelbrot.Numerics.Size
  ( size
  ) where

import Mandelbrot.Numerics.Complex

size
  :: (RealFloat r, Square r, Square (Complex r))
  => Int -> Complex r -> Complex r
{-# SPECIALIZE size :: Int -> Complex Double -> Complex Double #-}
size !p !c
  | p == 1 = 1
  | otherwise = go 2 c 0 (double c)
  where
    go !q !z !b !l
      | q == p = recip ((b + recip l) * sqr l)
      | otherwise = go (q + 1) (sqr z + c) (b + recip l) (double (z * l))
