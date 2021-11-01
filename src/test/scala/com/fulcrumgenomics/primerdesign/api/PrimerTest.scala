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
import com.fulcrumgenomics.util.Io
import org.scalatest.OptionValues

class PrimerTest extends UnitSpec with OptionValues {

  private val primers = Seq(
    Primer(bases="GATTACA", tm=12.34, penalty=56.78, mapping=Mapping("chr1", 1, 7)),
    Primer(bases="TGTAATC", tm=87.65, penalty=43.21, mapping=Mapping("chr22", 100, 106)),
    Primer(bases="",        tm=12.34, penalty=56.78, mapping=Mapping("chr1", 1, 1000)),
    Primer(bases="GGGGGGG", tm=12.34, penalty=56.78, mapping=Mapping("chr1", 1, 7, Strand.Negative)),
  )

  private val header = makeSamHeader
  private val dict   = makeSequenceDictionary

  "PrimerLike.lt" should "return a function that correctly orders primers" in {
    val lt     = Primer.lt(dict)
    val sorted = this.primers.sortWith(lt)
    sorted.map(_.bases) should contain theSameElementsInOrderAs Seq("GATTACA", "GGGGGGG", "", "TGTAATC")
  }

  "Primer" should "not be built when there are bases and the bases' length does not match the mapping length" in {
    an[Exception] should be thrownBy Primer(bases="GATTACA", tm=0, penalty=0, mapping=Mapping("chr1", 1, 1000))
  }

  it should "give the gc percentage" in {
    val expectedValues: Seq[Double] = Seq(28.57, 28.57, 0, 100.0)
    this.primers should have length expectedValues.length
    this.primers.map(_.gc).zip(expectedValues).foreach { case (actual, expected) =>
        actual shouldBe expected +- 0.01
    }
  }

  it should "give the length" in {
    val expectedValues: Seq[Int] = Seq(7, 7, 1000, 7)
    this.primers should have length expectedValues.length
    this.primers.map(_.length) should contain theSameElementsInOrderAs expectedValues
  }

  it should "give the maximum homopolymer length" in {
    val expectedValues: Seq[Int] = Seq(2, 2, 0, 7)
    this.primers should have length expectedValues.length
    this.primers.map(_.longestHomopolymer) should contain theSameElementsInOrderAs expectedValues
  }

  it should "give a concise string representation" in {
    val expectedStrings = Seq(
      "GATTACA\t12.34\t56.78\tchr1:1-7:+",
      "TGTAATC\t87.65\t43.21\tchr22:100-106:+",
      "*\t12.34\t56.78\tchr1:1-1000:+",
      "GGGGGGG\t12.34\t56.78\tchr1:1-7:-"
    )
    this.primers should have length expectedStrings.length
    this.primers.map(_.toString) should contain theSameElementsInOrderAs expectedStrings
  }

  it should "give a detailed BED representation given a custom name" in {
    val expectedStrings = Seq(
      s"chr1\t0\t7\tfulcrum_chr1_1_7_F\t500\t+\t0\t7\t100,100,100\t1\t7\t0",
      s"chr22\t99\t106\tfulcrum_chr22_100_106_F\t500\t+\t99\t106\t100,100,100\t1\t7\t0",
      s"chr1\t0\t1000\tfulcrum_chr1_1_1000_F\t500\t+\t0\t1000\t100,100,100\t1\t1000\t0",
      s"chr1\t0\t7\tfulcrum_chr1_1_7_R\t500\t-\t0\t7\t100,100,100\t1\t7\t0"
    )
    this.primers should have length expectedStrings.length

    this.primers.map(_.copy(namePrefix=Some("fulcrum")).toBed12Row).zip(expectedStrings).foreach { case (act, exp) =>
      act shouldBe exp
    }
  }

  it should "compare primers based on the mapping" in {
    this.primers.combinations(2).foreach { case Seq(left, right) =>
        Primer.compare(left, right, dict) shouldBe Mapping.compare(left.mapping, right.mapping, this.dict)
    }
  }

  it should "write primers to a BED file" in {
    val path = makeTempFile("primers.", ".sam")
    Primer.toBed(path, this.primers)
    val lines = Io.readLines(path).toSeq
    lines.length shouldBe this.primers.length + 1
    lines.headOption.value should include ("track")
    this.primers.zip(lines.drop(1)).foreach { case (primer, line) =>
        line shouldBe primer.toBed12Row
        line.split('\t') should have length 12
    }
  }
}
