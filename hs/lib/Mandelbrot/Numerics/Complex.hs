{- |
Complex numbers with operations useful in Mandelbrot set numeric algorithms.
-}
module Mandelbrot.Numerics.Complex
  ( module Data.Complex
  , magnitudeSquared
  , finite
  , notFinite
  , finiteC
  , notFiniteC
  , Square(..)
  ) where

import Data.Complex

magnitudeSquared :: Num r => Complex r -> r
magnitudeSquared (x :+ y) = x * x + y * y

notFiniteC :: RealFloat r => Complex r -> Bool
notFiniteC = not . finiteC

finiteC :: RealFloat r => Complex r -> Bool
finiteC (x :+ y) = finite x && finite y

finite :: RealFloat r => r -> Bool
finite = not . notFinite

notFinite :: RealFloat r => r -> Bool
notFinite x = isNaN x || isInfinite x

class Num t => Square t where
  sqr :: t -> t
  sqr x = x * x
  double :: t -> t
  double x = 2 * x

instance (RealFloat t, Square t) => Square (Complex t) where
  sqr (x:+y) = (sqr x - sqr y) :+ double (x * y)
  double (x :+ y) = double x :+ double y

instance Square Double
