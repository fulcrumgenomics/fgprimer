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
 */

package com.fulcrumgenomics.primerdesign

import com.fulcrumgenomics.commons.util.LazyLogging
import com.fulcrumgenomics.primerdesign.primer3.NtThermoAlign

import java.io.Closeable

/**
  * Class for checking for primer dimers that wraps ntthal.
  */
class DimerChecker(private val ntthal: NtThermoAlign = new NtThermoAlign()) extends LazyLogging with Closeable {

  /** Calculates the Tm of the worst duplex formed by the two sequences. */
  def tmOf(s1: String, s2: String): Double = {
    val pair = if (s1 < s2) (s1, s2) else (s2, s1)
    ntthal.duplexTm(s1, s2)
  }

  /** Returns the count of dimers formed between the query primer and target primers. */
  def countDimers(query: String, targets: Iterable[String], minTm: Double): Int = {
    targets.map(t => tmOf(query, t)).count(_ >= minTm)
  }

  override def close(): Unit = {
    this.ntthal.close()
  }
}
