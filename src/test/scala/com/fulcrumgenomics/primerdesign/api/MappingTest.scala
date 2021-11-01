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

import com.fulcrumgenomics.fasta.{SequenceDictionary, SequenceMetadata}
import com.fulcrumgenomics.primerdesign.testing.UnitSpec
import htsjdk.samtools.SAMSequenceRecord
import htsjdk.samtools.util.CoordMath

/**
  * Tests for the mapping class.
  */
class MappingTest extends UnitSpec {
  "Mapping" should "accept valid coordinate ranges and reject invalid ones" in {
    an[AssertionError] shouldBe thrownBy { Mapping("chr1", -1, 10, Strand.Positive) }
    an[AssertionError] shouldBe thrownBy { Mapping("chr1", 1, -1, Strand.Positive) }
    Mapping("chr1", 1, 10, Strand.Negative)
    Mapping("chr2", 100, 99, Strand.Positive) // allowed to represent zero-width intervals
  }

  it should "return the length" in {
    Mapping("chr1", 1, 1).length shouldBe 1
    Mapping("chr1", 50, 100).length shouldBe 51
  }

  it should "make valid sub-mappings" in {
    val m = Mapping("chr1", 5000, 5999, Strand.Positive)
    m.resolve(start=1, length=10)  shouldBe Mapping("chr1", 5000, 5009, Strand.Positive)
    m.resolve(start=1, length=501) shouldBe Mapping("chr1", 5000, 5500, Strand.Positive)
    m.resolve(start=501, length=1) shouldBe Mapping("chr1", 5500, 5500, Strand.Positive)
  }

  it should "disallow invalid submappings" in {
    val m = Mapping("chr1", 5000, 5999, Strand.Positive)
    an[AssertionError] shouldBe thrownBy { m.resolve(start = -1,     length=10) }
    an[AssertionError] shouldBe thrownBy { m.resolve(start = 10000,  length=10) }
    an[AssertionError] shouldBe thrownBy { m.resolve(start = 950,    length=100) }
  }

  it should "correctly project positions into a containing mapping" in {
    val m = Mapping("chr1", 1000, 1999, Strand.Positive)
    m.project(1000) shouldBe 1
    m.project(1999) shouldBe 1000
    m.project(1100) shouldBe 101
  }

  it should "refuse to project positions into a non-containing mapping" in {
    val m = Mapping("chr1", 1000, 1999, Strand.Positive)
    an[AssertionError] shouldBe thrownBy { m.project(-1) }
    an[AssertionError] shouldBe thrownBy { m.project(1) }
    an[AssertionError] shouldBe thrownBy { m.project(500) }
    an[AssertionError] shouldBe thrownBy { m.project(999) }
    an[AssertionError] shouldBe thrownBy { m.project(2000) }
    an[AssertionError] shouldBe thrownBy { m.project(10000) }
  }

  it should "correctly determine mapping that do and don't overlap" in {
    Mapping("chr1", 500, 1000).overlaps(Mapping("chr1",  500,  1000)) shouldBe true
    Mapping("chr1", 500, 1000).overlaps(Mapping("chr1",  600,   700)) shouldBe true
    Mapping("chr1", 500, 1000).overlaps(Mapping("chr1",  400,  2000)) shouldBe true
    Mapping("chr1", 500, 1000).overlaps(Mapping("chr1",  750,   850)) shouldBe true
    Mapping("chr1", 500, 1000).overlaps(Mapping("chr1", 1000,  1001)) shouldBe true

    Mapping("chr1", 500, 1000).overlaps(Mapping("chr2",  500,  1000)) shouldBe false
    Mapping("chr1", 500, 1000).overlaps(Mapping("chr1", 5000, 10000)) shouldBe false
    Mapping("chr1", 500, 1000).overlaps(Mapping("chr1",    1,   499)) shouldBe false
    Mapping("chr1", 500, 1000).overlaps(Mapping("chr1", 1001,  1100)) shouldBe false
  }

  it should "correctly compute overlap" in {
    Mapping("chr1", 500, 1000).lengthOfOverlapWith(Mapping("chr1",  500,  1000)) shouldBe CoordMath.getLength(500, 1000)
    Mapping("chr1", 500, 1000).lengthOfOverlapWith(Mapping("chr1",  600,   700)) shouldBe CoordMath.getLength(600, 700)
    Mapping("chr1", 500, 1000).lengthOfOverlapWith(Mapping("chr1",  400,  2000)) shouldBe CoordMath.getLength(500, 1000)
    Mapping("chr1", 500, 1000).lengthOfOverlapWith(Mapping("chr1",  750,   850)) shouldBe CoordMath.getLength(750, 850)
    Mapping("chr1", 500, 1000).lengthOfOverlapWith(Mapping("chr1", 1000,  1001)) shouldBe CoordMath.getLength(1000, 1000)

    Mapping("chr1", 500, 1000).lengthOfOverlapWith(Mapping("chr2",  500,  1000)) shouldBe 0
    Mapping("chr1", 500, 1000).lengthOfOverlapWith(Mapping("chr1", 5000, 10000)) shouldBe 0
    Mapping("chr1", 500, 1000).lengthOfOverlapWith(Mapping("chr1",    1,   499)) shouldBe 0
    Mapping("chr1", 500, 1000).lengthOfOverlapWith(Mapping("chr1", 1001,  1100)) shouldBe 0
  }

  it should "correctly determine if one mapping fully contains the other" in {
    Mapping("chr1", 500, 1000).contains(Mapping("chr1",  500,  1000)) shouldBe true
    Mapping("chr1", 500, 1000).contains(Mapping("chr1",  600,  1000)) shouldBe true
    Mapping("chr1", 500, 1000).contains(Mapping("chr1",  500,  900)) shouldBe true
    Mapping("chr1", 500, 1000).contains(Mapping("chr1",  600,  900)) shouldBe true

    Mapping("chr1", 500, 1000).contains(Mapping("chr2",  500,  1000)) shouldBe false
    Mapping("chr1", 500, 1000).contains(Mapping("chr1",  499,  1000)) shouldBe false
    Mapping("chr1", 500, 1000).contains(Mapping("chr2",  500,  1001)) shouldBe false
    Mapping("chr1", 500, 1000).contains(Mapping("chr2",  499,  1001)) shouldBe false
  }

  it should "union two mappings" in {
    Seq(
      (Mapping("chr1", 500, 1000), Mapping("chr1", 500, 1000), Mapping("chr1", 500, 1000)), // the same mapping
      (Mapping("chr1", 500, 1000), Mapping("chr1", 700, 800),  Mapping("chr1", 500, 1000)), // a contains b
      (Mapping("chr1", 700, 800),  Mapping("chr1", 500, 1000), Mapping("chr1", 500, 1000)), // b contains a
      (Mapping("chr1", 500, 800),  Mapping("chr1", 700, 1000), Mapping("chr1", 500, 1000)),  // overlap but do not contain
      (Mapping("chr1", 500, 800),  Mapping("chr1", 801, 1000), Mapping("chr1", 500, 1000))  // abuts
    ).foreach { case (left, right, expected) =>
        left.union(right) shouldBe expected
        right.union(left) shouldBe expected
    }
  }

  it should "throw an exception when unioning mappings on difference references" in {
    an[Exception] should be thrownBy Mapping("chr1", 500, 1000).union(Mapping("chr2", 500, 1000))
    an[Exception] should be thrownBy Mapping("chr2", 500, 1000).union(Mapping("chr1", 500, 1000))
  }

  it should "throw an exception when unioning mappings that do not overlap" in {
    an[Exception] should be thrownBy Mapping("chr1", 500, 1000).union(Mapping("chr1", 1500, 2000))
    an[Exception] should be thrownBy Mapping("chr1", 1500, 2000).union(Mapping("chr1", 500, 1000))
  }

  it should "shift mappings" in {
    Mapping("chr1", 500, 1000).shift(10) shouldBe Mapping("chr1", 510, 1010)
    Mapping("chr1", 500, 1000).shift(-10) shouldBe Mapping("chr1", 490, 990)
  }

  it should "throw an exception if the mapping is shifted before the beginning of the reference" in {
    an[Exception] should be thrownBy Mapping("chr1", 500, 1000).shift(-500)
    an[Exception] should be thrownBy Mapping("chr1", 500, 1000).shift(-1000)
  }

  it should "return the five prime position in the primer" in {
    Mapping("chr1", 500, 1000, Strand.Positive).fivePrimePosition shouldBe 500
    Mapping("chr1", 500, 1000, Strand.Negative).fivePrimePosition shouldBe 1000
  }

  it should "return a string representation" in {
    Mapping("chr1", 500, 1000, Strand.Positive).toString shouldBe "chr1:500-1000"
    Mapping("chr1", 500, 1000, Strand.Negative).toString shouldBe "chr1:500-1000"
    Mapping("chr1", 1000, 1000, Strand.Negative).toString shouldBe "chr1:1000"
    Mapping("chr1", 1000, 999, Strand.Negative).toString shouldBe "chr1:1000-999"
  }

  it should "return a string representation (with strand)" in {
    Mapping("chr1", 500, 1000, Strand.Positive).toStringWithStrand shouldBe "chr1:500-1000:+"
    Mapping("chr1", 500, 1000, Strand.Negative).toStringWithStrand shouldBe "chr1:500-1000:-"
    Mapping("chr1", 1000, 1000, Strand.Negative).toStringWithStrand shouldBe "chr1:1000:-"
    Mapping("chr1", 1000, 999, Strand.Negative).toStringWithStrand shouldBe "chr1:1000-999:-"
  }

  it should "return a mapping in BED format" in {
    Mapping("chr1", 500, 1000, Strand.Positive).toBed shouldBe "chr1\t499\t1000\tchr1:500-1000,501\t500\t+"
    Mapping("chr1", 500, 1000, Strand.Negative).toBed shouldBe "chr1\t499\t1000\tchr1:500-1000,501\t500\t-"
  }

  private val mappingTestCases: Seq[(Mapping, Mapping, Int)] = Seq(
    (Mapping("chr1", 1, 2, Strand.Positive), Mapping("chr1", 1, 2, Strand.Positive), 0),  // same
    (Mapping("chr1", 1, 2, Strand.Negative), Mapping("chr1", 1, 2, Strand.Negative), 0),  // same
    (Mapping("chr1", 1, 2, Strand.Positive), Mapping("chr1", 2, 2, Strand.Positive), -1), // start is later
    (Mapping("chr1", 1, 2, Strand.Positive), Mapping("chr1", 1, 3, Strand.Positive), -1), // end is later
    (Mapping("chr1", 1, 2, Strand.Positive), Mapping("chr1", 1, 2, Strand.Negative), -1), // strand is later
  )

  it should "be comparable within a reference sequence (without a sequence dictionary)" in {
    this.mappingTestCases.foreach { case (left, right, result) =>
      left.compare(right) shouldBe result
      right.compare(left) shouldBe -result
    }
  }

  it should "throw an exception when comparing across reference sequences (without a sequence dictionary)" in {
    an[Exception] should be thrownBy Mapping("chr1", 1, 2, Strand.Positive).compare(Mapping("chr2", 1, 2, Strand.Positive))
  }

  it should "be comparable across a reference sequence (without a sequence dictionary)" in {
    val dict = makeSequenceDictionary

    val extraTestCases = Seq(
      (Mapping("chr1", 1, 2, Strand.Positive), Mapping("chr2", 1, 2, Strand.Positive), -1),
      (Mapping("chr1", 10, 12, Strand.Positive), Mapping("chr2", 1, 2, Strand.Positive), -1)
    )
    (this.mappingTestCases ++ extraTestCases).foreach { case (left, right, result) =>
      Mapping.compare(left, right, dict) shouldBe result
      Mapping.compare(right, left, dict) shouldBe -result
    }
  }
}
