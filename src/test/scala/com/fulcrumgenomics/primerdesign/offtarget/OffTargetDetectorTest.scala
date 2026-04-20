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

package com.fulcrumgenomics.primerdesign.offtarget

import com.fulcrumgenomics.primerdesign.api.{Mapping, Primer, PrimerPair, Strand}
import com.fulcrumgenomics.primerdesign.testing.UnitSpec
import htsjdk.samtools.util.SequenceUtil

/**
  * Tests for the PCR Simulator
  */
class OffTargetDetectorTest extends UnitSpec {

  private val ref = testResource("miniref.fa", classOf[OffTargetDetector].getPackage)

  def pp(left: String, right: String) : PrimerPair = {
    val lp = Primer(left,  tm=65, penalty=0, Mapping("chr1", 500, 500+left.length-1,  Strand.Positive))
    val rp = Primer(right, tm=65, penalty=0, Mapping("chr1", 750, 750+right.length-1, Strand.Negative))
    PrimerPair(left=lp, right=rp, amplicon=Mapping("chr1", 500, 750+right.length-1), "", 85, penalty=0)
  }

  "OffTargetDetector" should "find a single mapping for a primer pair" in {
    val pcr = new OffTargetDetector(
      ref                             = ref,
      maxMismatches                   = 2,
      maxPrimerHits                   = 100,
      maxPrimerPairHits               = 1,
      threePrimeRegionLength          = 10,
      maxMismatchesInThreePrimeRegion = 2,
      maxAmpliconSize                 = 450
    )
    val primerPair = pp("GGCTAGAGTGCAGTGGTGCGATCT", SequenceUtil.reverseComplement("TACCGTGCCTGGCTGATTGCCT"))
    val actual     = pcr.check(primerPair)
    val expected   = OffTargetResult(primerPair, passes=true, mappings=Seq(Mapping("chr1", 781, 1042)))
    actual shouldBe expected
  }
}
