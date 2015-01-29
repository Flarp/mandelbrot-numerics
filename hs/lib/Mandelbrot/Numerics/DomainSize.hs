{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE FlexibleContexts #-}
module Mandelbrot.Numerics.DomainSize
  ( domainSize
  ) where

import Mandelbrot.Numerics.Complex

domainSize
  :: (RealFloat r, Square r, Square (Complex r))
  => Int -> Complex r -> r
{-# SPECIALIZE domainSize :: Int -> Complex Double -> Double #-}
domainSize p c = go 1 (magnitudeSquared c) c 1
  where
    go q !zq2 !z !dc
      | q >= p = sqrt zq2 / magnitude dc
      | otherwise = go (q + 1) zq2' (sqr z + c) (double (z * dc) + 1)
      where
        zq2' = zq2 `min` magnitudeSquared z
