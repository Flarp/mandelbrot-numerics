{-# LANGUAGE BangPatterns #-}
{-# LANGUAGE FlexibleContexts #-}
module Mandelbrot.Numerics.Child
  ( Child(..)
  , child
  , children
  , childSize
  ) where

import Data.Ratio
  ( (%)
  , denominator
  )
import Data.Strict.Tuple
  ( Pair(..)
  )

import Mandelbrot.Numerics.Complex
import Mandelbrot.Numerics.Interior
import Mandelbrot.Numerics.Nucleus
import Mandelbrot.Numerics.Shape
import Mandelbrot.Numerics.Size

data Child r = Child !Int !(Complex r)
  deriving (Eq, Read, Show)

childSize
  :: (RealFloat r, Square r, Square (Complex r), Approx r, Approx (Complex r))
  => Int -> Complex r -> Rational -> Maybe r
{-# SPECIALIZE childSize :: Int -> Complex Double -> Rational -> Maybe Double #-}
childSize !p !c0 !t = case shape p c0 of
  Just Cardioid | finite s -> Just $! s * sin (pi * fromRational t) * 2
  Just Circle | finite s -> Just $! s
  _ -> Nothing
  where
    !s = magnitude (size p c0) / (q * q)
    !q = fromInteger (denominator t)

child
  :: (RealFloat r, Square r, Square (Complex r), Approx r, Approx (Complex r))
  => Int -> Complex r -> Rational -> Maybe (Child r)
{-# SPECIALIZE child :: Int -> Complex Double -> Rational -> Maybe (Child Double) #-}
child p c0 t = case childSize p c0 t of
  Just s -> case converge 2 64 $ interior p (l * b) c0 c0 of
    Just (_:!:root1) -> case converge 2 64 $ interior p b c0 c0 of
      Just (_:!:root0) ->
         let n = root0 + (0.5 * s :+ 0) * signum (root0 - root1)
         in  case converge 2 64 $ nucleus p1 n of
              Just c1 -> Just (Child p1 c1)
              _ -> Nothing
      _ -> Nothing
    _ -> Nothing
  _ -> Nothing
  where
    l = fromRational $ (denominator t - 1) % (denominator t)
    b = cis (2 * pi * fromRational t)
    p1 = p * fromInteger (denominator t)

children
  :: (RealFloat r, Square r, Square (Complex r), Approx r, Approx (Complex r))
  => Int -> Complex r -> [Rational] -> [Child r]
{-# SPECIALIZE children :: Int -> Complex Double -> [Rational] -> [Child Double] #-}
children _ _ [] = []
children p0 c0 (t:ts) = case child p0 c0 t of
  Just c@(Child p1 c1) -> c : children p1 c1 ts
  _ -> []
