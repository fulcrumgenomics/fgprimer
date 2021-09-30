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

package com.fulcrumgenomics.primerdesign.offtarget

import com.fulcrumgenomics.FgBioDef.{FilePath, PathToFasta, unreachable}
import com.fulcrumgenomics.commons.util.Logger
import com.fulcrumgenomics.fasta.{SequenceDictionary, Topology}
import com.fulcrumgenomics.primerdesign.api.{Mapping, PrimerPair}
import com.fulcrumgenomics.primerdesign.offtarget.BwaAlnInteractive.Query
import com.fulcrumgenomics.primerdesign.offtarget.OffTargetResult.NoMappings
import com.fulcrumgenomics.util.ProgressLogger
import htsjdk.samtools.util.CoordMath

import java.io.Closeable
import scala.collection.mutable

object OffTargetResult {
  val NoMappings: Seq[Mapping] = Seq.empty
}

/**
  * Information obtained by running a single primer pair through the off-target detector.
  *
  * @param primerPair the primer pair submitted
  * @param passes true if the primer pair passes all checks, false otherwise
  * @param mappings the set of mappings of the primer pair to the genome or an empty Seq if mappings were not retained
  * @param leftPrimerMappings the set of mappings for the left primer, independent of the pair mappings, or NoMappings
  * @param rightPrimerMappings the set of mappings for the right primer, independent of the pair mappings, or NoMappings
  */
case class OffTargetResult(primerPair: PrimerPair,
                           passes: Boolean,
                           mappings: Seq[Mapping] = NoMappings,
                           leftPrimerMappings: Seq[Mapping]  = NoMappings,
                           rightPrimerMappings: Seq[Mapping] = NoMappings)

/**
  * A class for detecting off-target mappings of primers and primer pairs that uses a custom version of "bwa aln".
  * The off-target detection is faster and more sensitive than traditional isPCR and in addition can correctly detect
  * primers that are repetitive and contain many thousands or millions of mappings to the genome.
  *
  * Note that while this class invokes BWA with multiple threads, it is not itself thread-safe.  Only one thread at
  * a time should invoke methods on this class without external synchronization.
  *
  * @param ref the reference genome fasta file (must be index with BWA)
  * @param maxPrimerHits the maximum number of hits an individual primer can have in the genome before it is
  *                      considered an invalid primer, and all primer pairs containing the primer failed.
  * @param maxPrimerPairHits the maximum number of amplicons a primer pair can make and be considered passing
  * @param threePrimeRegionLength the number of bases at the 3' end of the primer in which the parameter
  *                               maxMismatchesInThreePrimeRegion is evaluated
  * @param maxMismatchesInThreePrimeRegion the maximum number of mismatches that are tolerated in the
  *                                        three prime region of each primer defined by threePrimeRegionLength
  * @param maxMismatches the maximum number of mismatches allowed in the full length primer (including any in the
  *                      three prime region)
  * @param maxAmpliconSize the maximum amplicon size to consider amplifiable
  * @param cacheResults if true, cache results for faster re-querying
  * @param threads the number of threads to use when invoking bwa
  * @param keepMappings if true, [[OffTargetResult]] objects will have amplicon mapping populated, otherwise not
  * @param keepPrimerMappings if true, [[OffTargetResult]] objects will have left and right primer mappings
  * @param logger if provided, use the given logger to produce info/debug logging during processing
  */
class OffTargetDetector(val bwaExecutable: FilePath = BwaAlnInteractive.DefaultBwaAlnInteractiveExecutable,
                        val ref: PathToFasta,
                        val maxPrimerHits: Int,
                        val maxPrimerPairHits: Int,
                        val threePrimeRegionLength: Int,
                        val maxMismatchesInThreePrimeRegion: Int,
                        val maxMismatches:Int,
                        val maxAmpliconSize: Int,
                        val cacheResults: Boolean = true,
                        val threads: Int = 16,
                        val keepMappings: Boolean = true,
                        val keepPrimerMappings: Boolean = false,
                        val logger: Option[Logger] = None
                  ) extends Closeable {

  private val dict = SequenceDictionary.extract(ref)
  private val primerCache: mutable.Map[String, BwaAlnInteractive.Result] = mutable.HashMap()
  private val primerPairCache: mutable.Map[PrimerPair, OffTargetResult] = mutable.HashMap()

  private val bwa = new BwaAlnInteractive(
    bwa                 = bwaExecutable,
    ref                 = ref,
    reverseComplement   = true,
    threads             = threads,
    seedLength          = threePrimeRegionLength,
    maxMismatchesInSeed = maxMismatchesInThreePrimeRegion,
    maxMismatches       = maxMismatches,
    maxHits             = maxPrimerHits
  )

  /** Checks a PrimerPair for off-target sites in the genome at which it might amplify. */
  def check(pp: PrimerPair): OffTargetResult = {
    check(Some(pp)).getOrElse(pp, unreachable("Where did our primer pair go!"))
  }

  /**
    * Checks a collection of primer pairs for off-target sites, returning a map of PrimerPair => Result.
    *
    * @param pps zero or more PrimerPairs
    * @return a Map containing all given primer pairs as keys with the values being the result of off-target checking
    */
  def check(pps: Iterable[PrimerPair]) : Map[PrimerPair, OffTargetResult] = {
    logger.foreach(_.info(f"Mapping ${pps.size}%,d primer pairs."))

    // Check the cache first
    val primerPairResults = mutable.HashMap[PrimerPair, OffTargetResult]()

    // Get the primer pairs to map.  If the primer pair is found in the cache, use that
    val primerPairsToMap = if (!cacheResults) pps else pps.filter { primerPair =>
      this.primerPairCache.get(primerPair) match {
        case None         => true  // Map it!
        case Some(result) =>
          primerPairResults(primerPair) = result
          false
      }
    }

    if (cacheResults) logger.foreach(_.info(f"Found ${primerPairResults.size}%,d cached results."))

    // Then handle any cache misses
    if (primerPairResults.size < pps.size) {
      // Get the primers to map
      val primersToMap = primerPairsToMap
        .flatMap { pp => Seq(pp.left.bases, pp.right.bases) }
        .toSet
        .toSeq
        .filterNot(b => cacheResults && this.primerCache.contains(b))

      // Build the unique list of queries to map with BWA
      val queries: Seq[Query] = primersToMap.map(bases => Query(id=bases, bases=bases))

      // Map the queries with BWA
      logger.foreach(_.info(f"Mapping ${queries.size}%,d queries."))
      val hitsByPrimer: Map[String, BwaAlnInteractive.Result] = this.bwa.map(queries).map(r => r.query.id -> r).toMap
      logger.foreach(_.info(f"Mapped ${hitsByPrimer.size}%,d queries."))
      if (cacheResults) this.primerCache ++= hitsByPrimer

      // Gets the mapping results for a single primer, looking in the cache first, otherwise in the results from BWA
      def getPrimerResult(bases: String): BwaAlnInteractive.Result = if (cacheResults) this.primerCache(bases) else hitsByPrimer(bases)

      // Go through all the primer pairs that need to be mapped
      val progress = logger.map(ProgressLogger(_, noun="primer pairs", unit=1e6.toInt))
      val numPassing = primerPairsToMap.count { pp =>
        // Get the mappings for the leftPrimerMappings and rightPrimerMappings primer respectively
        val p1: BwaAlnInteractive.Result = getPrimerResult(pp.left.bases)
        val p2: BwaAlnInteractive.Result = getPrimerResult(pp.right.bases)
        // Get all possible amplicons from the leftPrimerMappings and rightPrimerMappings primer hits, filtering if there are too many for either
        val result = if (p1.hitCount > maxPrimerHits || p2.hitCount > maxPrimerHits) {
          OffTargetResult(primerPair=pp, passes=false, mappings=Nil)
        }
        else {
          val amps = toAmplicons(p1.hits, p2.hits, this.maxAmpliconSize)

          OffTargetResult(
            primerPair          = pp,
            passes              = amps.lengthCompare(maxPrimerPairHits) <= 0,
            mappings            = if (keepMappings) amps else Seq.empty,
            leftPrimerMappings  = if (keepPrimerMappings) p1.hits.map(hitToMapping) else Seq.empty,
            rightPrimerMappings = if (keepPrimerMappings) p2.hits.map(hitToMapping) else Seq.empty
          )
        }

        progress.foreach(_.record())
        if (cacheResults) this.primerPairCache(pp) = result
        primerPairResults(pp) = result
        result.passes
      }
      logger.foreach(_.info(f"Found $numPassing%,d passing (on-target) primers out of ${primerPairResults.size}%,d."))
    }

    primerPairResults.toMap
  }

  /** Takes a set of hits for a one or more left primers and right primers and constructs Amplicon mappings anywhere
    * a left primer hit and a right primer hit align in F/R orientation up to `maxLen` apart on the same chromosome. */
  def toAmplicons(lefts: Iterable[BwaAlnInteractive.Hit], rights: Iterable[BwaAlnInteractive.Hit], maxLen: Int): Seq[Mapping] = {
    val amplicons = Seq.newBuilder[Mapping]
    for (h1 <- lefts; h2 <- rights) {
      if (h1.negative != h2.negative && h1.chrom == h2.chrom) {
        val (plus, minus) = if (h1.negative) (h2, h1) else (h1, h2)

        if (minus.start > plus.end && CoordMath.getLength(plus.start, minus.end) <= maxLen) {
          amplicons += Mapping(plus.chrom, plus.start, minus.end)
        }
      }
    }

    amplicons.result()
  }

  /** Converts a Bwa Hit object to a Mapping. */
  private def hitToMapping(hit: BwaAlnInteractive.Hit): Mapping = new Mapping(refName = hit.chrom, start = hit.start, end = hit.end)

  /** Terminates the underlying BWA process. */
  override def close(): Unit = this.bwa.close()
}
