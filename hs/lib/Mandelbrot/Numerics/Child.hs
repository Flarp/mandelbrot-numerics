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
  :: (RealFloat r, Square r, Square (Complex r))
  => Int -> Complex r -> Rational -> Maybe r
childSize p c0 t = case shape p c0 of
  Just Cardioid | finite s -> Just $ s * sin (pi * fromRational t) * 2
  Just Circle | finite s -> Just $ s
  _ -> Nothing
  where
    s = magnitude (size p c0) / (q * q)
    q = fromInteger (denominator t)

child
  :: (RealFloat r, Square r, Square (Complex r))
  => Int -> Complex r -> Rational -> Maybe (Child r)
child p c0 t = case childSize p c0 t of
  Just s -> case interior' p (l * b) c0 c0 of
    Just (_:!:root1) -> case interior' p b c0 c0 of
      Just (_:!:root0) -> case nucleus' p1 (root0 + 0.5 * (s :+ 0) * normalize (root0 - root1)) of
        Just c1 -> Just (Child p1 c1)
        _ -> Nothing
      _ -> Nothing
    _ -> Nothing
  _ -> Nothing
  where
    l = fromRational $ (denominator t - 1) % (denominator t)
    b = cis (2 * pi * fromRational t)
    p1 = p * fromInteger (denominator t)
    interior' q i z c = converge2 $ interior q i z c
    nucleus' q c = converge $ nucleus q c
    converge = lastMay . filter finiteC . take 64
    converge2 = lastMay . filter finiteC2 . take 64
    lastMay [] = Nothing
    lastMay xs = Just (last xs)

normalize c = c / (magnitude c :+ 0)

finiteC2 (a:!:b) = finiteC a && finiteC b

children
  :: (RealFloat r, Square r, Square (Complex r))
  => Int -> Complex r -> [Rational] -> [Child r]
children _ _ [] = []
children p0 c0 (t:ts) = case child p0 c0 t of
  Just c@(Child p1 c1) -> c : children p1 c1 ts
  _ -> []
