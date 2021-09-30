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

package com.fulcrumgenomics.primerdesign.api

import com.fulcrumgenomics.primerdesign.testing.UnitSpec
import com.fulcrumgenomics.util.{Io, Sequences}
import org.scalatest.OptionValues

class PrimerPairTest extends UnitSpec with OptionValues {

  private val dict = makeSequenceDictionary

  private val leftPrimers = Seq(
    Primer(bases="GATTACA", tm=12.34, penalty=56.78, mapping=Mapping("chr1", 1, 7)),
    Primer(bases="TGTAATC", tm=87.65, penalty=43.21, mapping=Mapping("chr22", 100, 106)),
    Primer(bases="",        tm=12.34, penalty=56.78, mapping=Mapping("chr1", 1, 40)),
    Primer(bases="GGGGGGG", tm=12.34, penalty=56.78, mapping=Mapping("chr1", 1, 7)),
  )

  private val rightPrimers = leftPrimers.map { left =>
    val rightMapping = left.mapping.copy(
      start  = left.mapping.start + 50,
      end    = left.mapping.end + 50,
      strand = Strand.Negative
    )
    left.copy(
      bases = Sequences.revcomp(left.bases),
      tm    = left.tm + 10,
      penalty = left.penalty + 10,
      mapping = rightMapping
    )
  }

  private def toAmpliconSequence(left: String, right: String, length: Int): String = {
    if (left.isEmpty || right.isEmpty) ""
    else {
      // FIXME: does not add the leftPrimerMappings/rightPrimerMappings sequences to the amplicon
      val bases = left.zip(right).map { case (l, r) => s"$l$r" }.mkString.take(length)
      val times = (length / bases.length) + 1
      (bases * times).take(length)
    }
  }

  private val primerPairs = leftPrimers.zip(rightPrimers).map { case (left, right) =>
    val amplicon = left.mapping.copy(end = right.mapping.end)
    val ampliconSequence = toAmpliconSequence(left.bases, right.bases, amplicon.length)
    new PrimerPair(
      left             = left,
      right            = right,
      amplicon         = amplicon,
      ampliconSequence = ampliconSequence,
      tm               = left.tm + right.tm,
      penalty          = left.penalty + right.penalty
    )
  }

  "PrimerPair" should "require that the amplicon sequence is empty or that its length matches the amplicon mapping length" in {
    this.primerPairs.headOption.foreach { pp =>
      // OK: amplicon sequence is empty
      pp.copy(ampliconSequence = "")
      // OK: amplicon sequence length is the same as the amplicon mapping length.  Kind of circular check, but OK.
      val left  = Primer(bases="A", tm=0, penalty=0, mapping=Mapping("chr1", 1, 1))
      val right = Primer(bases="A", tm=0, penalty=0, mapping=Mapping("chr1", 3, 3))
      PrimerPair(left=left, right=right, amplicon=Mapping("chr1", 1, 3), ampliconSequence="AAA", tm=0, penalty=0)
      // NOK: amplicon sequence length differs from the amplicon mapping length
      an[Exception] should be thrownBy pp.copy(ampliconSequence = pp.ampliconSequence + "A")
      an[Exception] should be thrownBy pp.copy(amplicon = pp.amplicon.copy(start=pp.amplicon.start+1))
      // NOK: copy validates
      an[Exception] should be thrownBy pp.copy(ampliconSequence = pp.ampliconSequence + "A")
    }
  }

  it should "require that the amplicon is the leftPrimerMappings-most leftPrimerMappings-primer base to rightPrimerMappings-most rightPrimerMappings-primer base" in {
    this.primerPairs.headOption.foreach { pp: PrimerPair =>
      an[Exception] should be thrownBy pp.copy(right = pp.right.copy(mapping=pp.right.mapping.copy(end=pp.right.mapping.end+1)))
      an[Exception] should be thrownBy pp.copy(amplicon = pp.amplicon.copy(end = pp.amplicon.end + 1))
    }
  }

  it should "require that the primers are on the same reference" in {
    this.primerPairs.headOption.foreach { pp: PrimerPair =>
      an[Exception] should be thrownBy pp.copy(amplicon = pp.amplicon.copy(refName = "no-name"))
      an[Exception] should be thrownBy pp.copy(left = pp.left.copy(mapping = pp.left.mapping.copy(refName = "no-name")))
      an[Exception] should be thrownBy pp.copy(right = pp.right.copy(mapping = pp.left.mapping.copy(refName = "no-name")))
    }
  }

  it should "give the length" in {
    val expectedValues: Seq[Int] = Seq(57, 57, 90, 57)
    this.primerPairs should have length expectedValues.length
    this.primerPairs.map(_.length) should contain theSameElementsInOrderAs expectedValues
  }

  it should "give the inner mapping" in {
    val expectedMappings = Seq(
      Mapping("chr1",  8, 50),
      Mapping("chr22", 107, 149),
      Mapping("chr1",  41, 50),
      Mapping("chr1",  8, 50)
    )
    this.primerPairs should have length expectedMappings.length
    this.primerPairs.map(_.inner) should contain theSameElementsInOrderAs expectedMappings
  }

  it should "give the gc for the amplicon sequence" in {
    val expectedValues: Seq[Double] = Seq(29.82, 28.07, 0, 100.0)
    this.primerPairs should have length expectedValues.length
    this.primerPairs.map(_.gc).zip(expectedValues).foreach { case (actual, expected) =>
      actual shouldBe expected +- 0.01
    }
  }

  it should "give a detailed BED representation given a custom name" in {
    val expectedStrings = Seq(
      s"chr1\t0\t57\tfulcrum_chr1_1_57_F\t500\t+\t0\t57\t100,100,100\t3\t7,43,7\t0,7,50",
      s"chr22\t99\t156\tfulcrum_chr22_100_156_F\t500\t+\t99\t156\t100,100,100\t3\t7,43,7\t0,7,50",
      s"chr1\t0\t90\tfulcrum_chr1_1_90_F\t500\t+\t0\t90\t100,100,100\t3\t40,10,40\t0,40,50",
      s"chr1\t0\t57\tfulcrum_chr1_1_57_F\t500\t+\t0\t57\t100,100,100\t3\t7,43,7\t0,7,50"
    )
    this.primerPairs should have length expectedStrings.length

    this.primerPairs.map(_.copy(namePrefix=Some("fulcrum")).toBed12Row).zip(expectedStrings).foreach { case (act, exp) =>
      act shouldBe exp
    }
  }

  it should "give a concise string representation" in {
    val expectedStrings = Seq(
      "GATTACA\t12.34\t56.78\tchr1:1-7:+\tTGTAATC\t22.34\t66.78\tchr1:51-57:-\tGTAGTTTAAACTACGTAGTTTAAACTACGTAGTTTAAACTACGTAGTTTAAACTACG\t34.68\t123.56",
      "TGTAATC\t87.65\t43.21\tchr22:100-106:+\tGATTACA\t97.65\t53.21\tchr22:150-156:-\tTGGATTATAATCCATGGATTATAATCCATGGATTATAATCCATGGATTATAATCCAT\t185.3\t96.42",
      "*\t12.34\t56.78\tchr1:1-40:+\t*\t22.34\t66.78\tchr1:51-90:-\t*\t34.68\t123.56",
      "GGGGGGG\t12.34\t56.78\tchr1:1-7:+\tCCCCCCC\t22.34\t66.78\tchr1:51-57:-\tGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCGCG\t34.68\t123.56"
    )
    this.primerPairs should have length expectedStrings.length
    this.primerPairs.map(_.toString).zip(expectedStrings).foreach { case (act, exp) =>
      act shouldBe exp
    }
  }

  it should "compare based on the amplicon's mapping" in {
    this.primerPairs.combinations(2).foreach { case Seq(left, right) =>
      PrimerPair.compare(left, right, dict, byAmplicon=true) shouldBe Mapping.compare(left.amplicon, right.amplicon, this.dict)
    }
  }

  it should "compare based on the primer pair's leftPrimerMappings and rightPrimerMappings primers" in {
    this.primerPairs.combinations(2).foreach { case Seq(left, right) =>
      val leftCompare = Mapping.compare(left.left.mapping, right.left.mapping, this.dict)
      val expectedResult = if (leftCompare == 0) Mapping.compare(left.right.mapping, right.right.mapping, this.dict) else leftCompare
      PrimerPair.compare(left, right, dict, byAmplicon=false) shouldBe expectedResult
    }
  }

  it should "write primer pairs to a BED file" in {
    val path = makeTempFile("primers.", ".sam")
    PrimerPair.toBed(path, this.primerPairs)
    val lines = Io.readLines(path).toSeq
    lines.length shouldBe this.primerPairs.length + 1
    lines.headOption.value should include ("track")
    this.primerPairs.zip(lines.drop(1)).foreach { case (primer, line) =>
      line shouldBe primer.toBed12Row
      line.split('\t') should have length 12
    }
  }
}
