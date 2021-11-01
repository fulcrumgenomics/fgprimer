/*
 * The MIT License
 *
 * Copyright (c) 2021 Fulcrum Genomics LLC
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in
 * all copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
 * THE SOFTWARE.
 *
 */

package com.fulcrumgenomics.primerdesign

import com.fulcrumgenomics.commons.CommonsDef

import java.nio.file.Path

/** `Predef` like object that is intended to be imported with `import PrimerDesignDef._`
  * in order to import commonly used types and functions for primer design.
  */
object PrimerDesignDef extends CommonsDef {
  type PathToPrimers = Path
  type PathToPrimerPairs = Path

  /** An implicit to convert a [[Tuple3]] to a [[MinOptMax]]. */
  implicit def toMinOptMax[T](tuple: (T, T, T))(implicit numeric: Numeric[T]): MinOptMax[T] = {
    MinOptMax(tuple._1, tuple._2, tuple._3)
  }

  /** A case class that stores a minimum, optimal, and maximum.  Optimal should be zero if it should be ignored. */
  case class MinOptMax[T](min: T, opt: T, max: T)(implicit numeric: Numeric[T]) extends Iterable[T] {
    require(numeric.lteq(min, max))
    require(numeric.zero == opt || (numeric.lteq(min, opt) && numeric.lteq(opt, max)))
    def tuple: (T, T, T) = (min, opt, max)

    def range: T = numeric.plus(numeric.minus(max, min), numeric.fromInt(1))

    def toDouble: MinOptMax[Double] = {
      new MinOptMax[Double](
        min = numeric.toDouble(this.min),
        opt = numeric.toDouble(this.opt),
        max = numeric.toDouble(this.max)
      )
    }

    def iterator: Iterator[T] = Iterator(min, opt, max)

    override def toString: String = f"(${numeric.toDouble(min)}%.2f, ${numeric.toDouble(opt)}%.2f, ${numeric.toDouble(max)}%.2f)"
  }

  /** Splits a string into a pair of strings at the first index of the split character. */
  private[primerdesign] def splitInTwo(s: String, ch: Char) : (String, String) = {
    s.indexOf(ch) match {
      case -1 => throw new IllegalArgumentException(s"Value '$s' did not contain split character '$ch'.")
      case i  => (s.take(i), s.drop(i+1))
    }
  }
}
