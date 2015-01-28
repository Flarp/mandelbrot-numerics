{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE FlexibleContexts #-}
module Mandelbrot.Numerics.Interior
  ( interior
  ) where

import Data.Strict.Tuple
  ( Pair(..)
  )

import Mandelbrot.Numerics.Complex

interior
  :: (RealFloat r, Square r, Square (Complex r))
  => Int -> Complex r -> Complex r -> Complex r -> [Pair (Complex r) (Complex r)]
interior p i z0 c0
  | notFiniteC z0 || notFiniteC c0 = []
  | otherwise = go 0 z0 1 0 0 0
  where
    go q !z !dz !dc !dzdz !dcdz
      | q >= p = (z' :!: c') : interior p i z' c'
      | otherwise =
          go (q + 1) (sqr z + c0) (double (z * dz)) (double (z * dc) + 1)
            (double (z * dzdz + sqr dz)) (double (z * dcdz + dc * dz))
      where
        det = (dz - 1) * dcdz - dc * dzdz;
        z' = z0 - (dcdz * (z - z0) - dc * (dz - i)) / det;
        c' = c0 - ((dz - 1) * (dz - i) - dzdz * (z - z0)) / det;
