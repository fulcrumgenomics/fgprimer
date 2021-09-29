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

package com.fulcrumgenomics.primerdesign.primer3

import com.fulcrumgenomics.primerdesign.testing.UnitSpec
import htsjdk.samtools.util.SequenceUtil

class NtThermAlignTest extends UnitSpec {
  "NtThermoAlign" should "calculate a very low duplex Tm for non-matching sequences" in {
    new NtThermoAlign().duplexTm("AAAAAAAAAAAAAAAAAAAA", "CCCCCCCCCCCCCCCCCCCC") shouldBe 0
  }

  it should "calculate a high Tm for two sequences that are RCs of one another" in {
    val s1 = "CTGACTGACTTGAGTTCGCTA"
    val s2 = SequenceUtil.reverseComplement(s1)
    new NtThermoAlign().duplexTm(s1, s2) shouldBe 51.634492 +- 0.0001
  }

  it should "handle two sequences of quite different lengths" in {
    val s1 = "CTGACTGACTTGAGTTCGCTA"
    val s2 = SequenceUtil.reverseComplement(s1.substring(10))
    new NtThermoAlign().duplexTm(s1, s2) shouldBe 27.204056 +- 0.0001
  }
}
