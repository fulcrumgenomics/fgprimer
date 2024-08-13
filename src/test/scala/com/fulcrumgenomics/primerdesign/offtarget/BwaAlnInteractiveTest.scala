package com.fulcrumgenomics.primerdesign.offtarget

import com.fulcrumgenomics.primerdesign.offtarget.BwaAlnInteractive.{Hit, Query}
import com.fulcrumgenomics.primerdesign.testing.UnitSpec

class BwaAlnInteractiveTest extends UnitSpec {

  private val ref = testResource("miniref.fa", classOf[OffTargetDetector].getPackage)

  "Hit.apply" should "should create a new hit" in {
    // only the "negative" flag should be changed when reverse complementing
    val fwd = Hit(chrom = "chr1", start=20, negative=true, cigar="30M1D20M", edits=3, rc=false)
    val rc = Hit(chrom = "chr1", start=20, negative=true, cigar="30M1D20M", edits=3, rc=true)

    rc.negative shouldBe false
    rc.copy(negative=true) shouldBe fwd
  }

  "BwaAlnInteractive" should "map a single hit to the mini reference" in {
    val bases = "GGCTAGGTGCAGTGGTGCGATCT"
    val query = Query(id="test", bases=bases)

    Seq(true, false).foreach { reverseComplement =>
      val bwa = new BwaAlnInteractive(
        ref                 = ref,
        seedLength          = 10,
        maxHits             = 100,
        maxMismatches       = 1,
        maxMismatchesInSeed = 1,
        maxGapOpens         = 1,
        reverseComplement   = reverseComplement
      )
      val result = bwa.map(query = query.bases, id = query.id)

      result.query shouldBe query
      result.hitCount shouldBe 1
      val hit = result.hits.head
      hit.chrom shouldBe "chr1"
      hit.start shouldBe 781
      hit.negative shouldBe false
      hit.cigar.toString shouldBe "6M1D17M"
      hit.edits shouldBe 1
      hit.end shouldBe (781 + hit.cigar.getReferenceLength - 1)
      hit.mismatches shouldBe 0

      bwa.close()
    }
  }
}
