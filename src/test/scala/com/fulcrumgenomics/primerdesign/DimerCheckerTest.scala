/*
 * The MIT License
 *
 * Copyright (c) 2021 Fulcrum Genomics
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

import com.fulcrumgenomics.primerdesign.primer3.NtThermoAlign
import com.fulcrumgenomics.primerdesign.testing.UnitSpec

class DimerCheckerTest extends UnitSpec {
  "DimerChecker" should "find some dimers" in {
    // Four primers - first two are just rotations of each other; second two are reverse complements of each other
    val primers = IndexedSeq("AGGAGGAGGAGGAGGAGGAGG", "GAGGAGGAGGAGGAGGAGGAG", "GGAGCTGATCGCTAGCTGATA", "TATCAGCTAGCGATCAGCTCC")
    val checker = new DimerChecker(ntthal=new NtThermoAlign())

    checker.countDimers(primers.head, primers, minTm=30) shouldBe 0
    checker.countDimers(primers(1),   primers, minTm=30) shouldBe 0
    checker.countDimers("CCTCCTCCTCCTCCTCCTCCT", primers, minTm=30) shouldBe 2

    checker.countDimers(primers(2), primers, minTm=30) shouldBe 1
    checker.countDimers(primers(3), primers, minTm=30) shouldBe 1
  }
}
