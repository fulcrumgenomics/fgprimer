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

import com.fulcrumgenomics.primerdesign.PrimerDesignDef._
import com.fulcrumgenomics.primerdesign.VariantLookup
import com.fulcrumgenomics.primerdesign.api.{Mapping, Strand}
import com.fulcrumgenomics.primerdesign.primer3.Primer3FailureReason._
import com.fulcrumgenomics.primerdesign.testing.UnitSpec
import com.fulcrumgenomics.util.Sequences
import htsjdk.samtools.reference.ReferenceSequenceFileFactory

/**
  * Tests for the class that drives interaction with Primer3. This test requires that
  * 'primer3_core' in on the Path, otherwise it will fail.
  *
  * Uses the test 'miniref.fa', which contains two "chromosomes":
  * name   length
  * chr1      577
  * chr2    9,821
  */
class Primer3Test extends UnitSpec {
  private val genome = testResource("miniref.fa", classOf[Primer3].getPackage)
  private val vcf    = testResource("miniref.variants.vcf", classOf[Primer3].getPackage)
  private val target = Mapping("chr1", 201, 250, Strand.Positive)
  private val examplePrimer3Input = Primer3Parameters(
    howMany        = 10,
    ampliconSizes  = (100, 150, 200),
    ampliconTms    = (55, 0, 100).toDouble,
    primerSizes    = (20, 25, 30),
    primerTms      = (55.0, 65.0, 75.0),
    primerGcs      = (30, 45, 65),
    primerMaxPolyx = 4
  )

  {
    "Primer3.designPrimers" should "design some primer pairs" in {
      val ref    = ReferenceSequenceFileFactory.getReferenceSequenceFile(genome)
      val p3     = new Primer3(genome = genome, variantLookup = VariantLookup.empty)
      val target = Mapping("chr1", 201, 250, Strand.Positive)
      val params = Primer3Parameters(
        howMany        = 10,
        ampliconSizes  = (100, 150, 200),
        ampliconTms    = (55, 0, 100).toDouble,
        primerSizes    = (20, 25, 30),
        primerTms      = (55.0, 65.0, 75.0),
        primerGcs      = (30, 45, 65),
        primerMaxPolyx = 4
      )
      val result = p3.designPrimers(Primer3Input(target, DesignPrimerPairsTask, params))
      p3.close()

      val pairs = result.primerPairs

      pairs.size shouldBe 10
      pairs.foreach(pair => {
        pair.amplicon.length should be <= 200

        pair.primers.foreach(primer => {
          primer.length should (be >= 20 and be <= 30)
          primer.tm should (be >= 55.0 and be <= 75.0)
          primer.gc should (be >= 30.0 and be <= 65.0)
        })

        val leftFromRef = new String(ref.getSubsequenceAt(pair.left.mapping.refName, pair.left.mapping.start, pair.left.mapping.end).getBases)
        val rightFromRef = Sequences.revcomp(new String(ref.getSubsequenceAt(pair.right.mapping.refName, pair.right.mapping.start, pair.right.mapping.end).getBases))

        pair.left.bases.toUpperCase shouldBe leftFromRef.toUpperCase
        pair.right.bases.toUpperCase shouldBe rightFromRef.toUpperCase
      })
    }

    it should "design some leftPrimerMappings primers" in {
      val ref    = ReferenceSequenceFileFactory.getReferenceSequenceFile(genome)
      val p3     = new Primer3(genome = genome, variantLookup = VariantLookup.empty)
      val target = Mapping("chr1", 201, 250, Strand.Positive)
      val params = Primer3Parameters(
        howMany        = 1000,
        ampliconSizes  = (100, 200, 250),
        ampliconTms    = (55, 0, 100).toDouble,
        primerSizes    = (29, 30, 31),
        primerTms      = (63.0, 65.0, 67.0),
        primerGcs      = (30, 45, 65),
        primerMaxPolyx = 4
      )

      val result = p3.designPrimers(Primer3Input(target, DesignLeftPrimersTask, params))
      p3.close()

      val primers = result.primers
      primers.foreach(primer => {
        primer.length should (be >= 29 and be <= 31)
        primer.tm should (be >= 63.0 and be <= 67.0)
        primer.gc should (be >= 30.0 and be <= 65.0)
        primer.mapping.start should be < primer.mapping.end
        primer.mapping.end should be < 201

        val fromRef = new String(ref.getSubsequenceAt(primer.mapping.refName, primer.mapping.start, primer.mapping.end).getBases)
        primer.bases.toUpperCase shouldBe fromRef.toUpperCase
      })
    }

    it should "design some rightPrimerMappings primers" in {
      val ref    = ReferenceSequenceFileFactory.getReferenceSequenceFile(genome)
      val p3     = new Primer3(genome = genome, variantLookup = VariantLookup.empty)
      val target = Mapping("chr1", 201, 250, Strand.Positive)
      val params = Primer3Parameters(
        howMany        = 1000,
        ampliconSizes  = (100, 200, 250),
        ampliconTms    = (55, 0, 100).toDouble,
        primerSizes    = (29, 30, 31),
        primerTms      = (63.0, 65.0, 67.0),
        primerGcs      = (30, 45, 65),
        primerMaxPolyx = 4
      )
      val result = p3.designPrimers(Primer3Input(target, DesignRightPrimersTask, params))
      p3.close()

      val primers = result.primers
      primers.foreach(primer => {
        primer.length should (be >= 29 and be <= 31)
        primer.tm should (be >= 63.0 and be <= 67.0)
        primer.gc should (be >= 30.0 and be <= 65.0)
        primer.mapping.end   should be > primer.mapping.start
        primer.mapping.start should be > 250

        val fromRef = Sequences.revcomp(new String(ref.getSubsequenceAt(primer.mapping.refName, primer.mapping.start, primer.mapping.end).getBases))
        primer.bases.toUpperCase shouldBe fromRef.toUpperCase
      })
    }
  }

  /**
    * Complicated test that ensures that masking of common variants in the design sequence is working as intended.
    * The VCF contains a number of both common and rare SNPs, as well as some common insertions, deletions and
    * mixed events.
    */
  "Primer3.designSequences" should "correctly pick the variants to mask when constructing design sequence" in {
    val region = Mapping(refName="chr2", start=9000, end=9110, strand=Strand.Positive)
    //                      9000      9010      9020      9030      9040      9050      9060      9070      9080      9090      9100      9110
    val expectedUnmasked = "AATATTCTTGCTGCTTATGCAGCTGACATTGTTGCCCTCCCTAAAGCAACCAAGTAGCCTTTATTTCCCACAGTGAAAGAAAACGCTGGCCTATCAGTTACATTACAAAAG"
    val expectedMasked   = "AATATTCTTGNTGCTTATGCNGCTGACATTGTTGCCCTCCCTAAAGCAACNAAGTAGCCTNTATTTCCCANAGTGAAAGANNACGCTGGCCNNTCAGTTANNNTACAAAAG"
    val p3 = new Primer3(genome = genome, variantLookup=VariantLookup.cached(vcf))

    val (actualUnmasked, actualMasked) = p3.designSequences(region=region)

    actualUnmasked shouldBe expectedUnmasked
    actualMasked   shouldBe expectedMasked
  }

  "Primer3.buildFailures" should "parse a single failure string and return the elements as failures" in {
    val p3 = new Primer3(genome = genome, variantLookup = VariantLookup.empty)
    p3.buildFailures(examplePrimer3Input, Nil)("considered 3576, lowercase masking of 3' end 3576, ok 0") shouldBe IndexedSeq(Primer3Failure(LowercaseMasking, 3576))

    p3.buildFailures(examplePrimer3Input, Nil)("considered 1998, GC content failed 10, GC clamp failed 105, low tm 181, high tm 179, lowercase masking of 3' end 1120, ok 393") shouldBe
      IndexedSeq(Primer3Failure(LowercaseMasking, 1120),
        Primer3Failure(LowTm, 181),
        Primer3Failure(HighTm, 179),
        Primer3Failure(GcClamp, 105),
        Primer3Failure(GcContent, 10))

    p3.buildFailures(examplePrimer3Input, Nil)("considered 1000, ok 1000") shouldBe IndexedSeq()
  }

  it should "parse a string with an unknown failure reason and ignore it" in {
    val p3 = new Primer3(genome = genome, variantLookup = VariantLookup.empty)
    p3.buildFailures(examplePrimer3Input, Nil)("considered 1000, wib-wobbled 100, ok 900") shouldBe IndexedSeq()
  }

  it should "parse a multiple failure strings and combine them appropriately" in {
    val p3 = new Primer3(genome = genome, variantLookup = VariantLookup.empty)
    val failures = p3.buildFailures(examplePrimer3Input, Nil)("considered 3285, GC clamp failed 16, low tm 24, long poly-x seq 12, lowercase masking of 3' end 3208, ok 25",
      "considered 2992, GC clamp failed 26, low tm 28, high tm 32, long poly-x seq 13, lowercase masking of 3' end 2824, ok 61")

    failures shouldBe IndexedSeq(
      Primer3Failure(LowercaseMasking, 2824+3208),
      Primer3Failure(LowTm, 24+28),
      Primer3Failure(GcClamp, 16+26),
      Primer3Failure(HighTm, 32),
      Primer3Failure(LongPolyX, 12+13)
    )
  }
}
