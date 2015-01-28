{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE FlexibleContexts #-}
module Mandelbrot.Numerics.Nucleus
  ( nucleus
  ) where

import Mandelbrot.Numerics.Complex

nucleus
  :: (RealFloat r, Square r, Square (Complex r))
  => Int -> Complex r -> [Complex r]
nucleus !p !c0
  | p < 1 || notFiniteC c0 = []
  | otherwise = go 0 0 0
  where
    go i z dc
      | i == p = c' : nucleus p c'
      | otherwise = go (i + 1) (sqr z + c0) (double (z * dc) + 1)
      where
        c' = c0 - z / dc
