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
import com.fulcrumgenomics.util.ProgressLogger
import htsjdk.samtools.util.CoordMath

import java.io.Closeable
import scala.collection.mutable

/**
  * Result type for running a single PrimerPair through the PcrSimulator.
  */
case class PrimerOffTargetDetectorResult(primerPair: PrimerPair,
                                         passes: Boolean,
                                         mappings: Seq[Mapping],
                                         left: Seq[Mapping] = Seq.empty,
                                         right: Seq[Mapping] = Seq.empty)

/**
  * A PCR Simulator that is a replacement for in-silico PCR and is both significantly faster
  * and more sensitive.
  *
  * @param ref the reference sequence (must be index with BWA)
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
  */
class PcrSimulator(val bwaExecutable: FilePath = BwaAlnInteractive.DefaultBwaAlnInteractiveExecutable,
                   val ref: PathToFasta,
                   val maxPrimerHits: Int = 1000,
                   val maxPrimerPairHits: Int = 1,
                   val threePrimeRegionLength: Int = 15,
                   val maxMismatchesInThreePrimeRegion: Int = 3,
                   val maxMismatches:Int = 3,
                   val maxAmpliconSize: Int = 600,
                   val cacheResults: Boolean = true,
                   val threads: Int = 8,
                   val keepMappings: Boolean = true,
                   val KeepPrimerMappings: Boolean = false,
                   val logger: Option[Logger] = None
                  ) extends Closeable {

  private val dict = SequenceDictionary.extract(ref)

  private val primerCache: mutable.Map[String, BwaAlnInteractive.Result] = mutable.HashMap()
  private val primerPairCache: mutable.Map[PrimerPair, PrimerOffTargetDetectorResult] = mutable.HashMap()

  private val bwa = new BwaAlnInteractive(bwa=bwaExecutable, ref=ref, reverseComplement=true, threads=threads,
    seedLength=threePrimeRegionLength, maxMismatchesInSeed=maxMismatchesInThreePrimeRegion,
    maxMismatches=maxMismatches, maxHits=maxPrimerHits)

  /** Checks a single PrimerPair for possible amplification sites in the genome. */
  def check(primerPair: PrimerPair): PrimerOffTargetDetectorResult = {
    check(Option(primerPair)).getOrElse(primerPair, unreachable("Must have results for the primer pair!"))
  }

  /**
    * Checks one or more PrimerPairs for expected amplicons within the genome.
    *
    * @param primerPairs one or more PrimerPairs
    * @return a Map that contains an entry per PrimerPair with the results of PCR simulations
    */
  def check(primerPairs: Iterable[PrimerPair]) : Map[PrimerPair, PrimerOffTargetDetectorResult] = {
    logger.foreach(_.info(f"Mapping ${primerPairs.size}%,d primer pairs."))

    // Check the cache first
    val primerPairResults = mutable.HashMap[PrimerPair, PrimerOffTargetDetectorResult]()

    // Get the primer pairs to map.  If the primer pair is found in the cache, use that
    val primerPairsToMap = if (!cacheResults) primerPairs else primerPairs.filter { primerPair =>
      this.primerPairCache.get(primerPair) match {
        case None         => true  // Map it!
        case Some(result) =>
          primerPairResults(primerPair) = result
          false
      }
    }

    if (cacheResults) logger.foreach(_.info(f"Found ${primerPairResults.size}%,d cached results."))

    // Then handle any cache misses
    if (primerPairResults.size < primerPairs.size) {

      // Get the primers to map
      val primersToMap = primerPairsToMap
        .flatMap { pp => Seq(pp.left.bases, pp.right.bases) }
        .toSet.toSeq
        .filterNot(b => cacheResults && this.primerCache.contains(b))

      // Build the unique list of queries to map with BWA
      val queries: Seq[Query] = primersToMap.map(bases => Query(id=bases, bases=bases))

      // Map the queries with BWA
      logger.foreach(_.info(f"Mapping ${queries.size}%,d queries."))
      val hitsByPrimer: Map[String, BwaAlnInteractive.Result] = this.bwa.map(queries).map(r => r.query.id -> r).toMap
      logger.foreach(_.info(f"Mapped ${hitsByPrimer.size}%,d queries."))
      // TODO handle hits to alts?

      // Gets the mapping results for a single primer, looking in the cache first, otherwise in the results from BWA
      def getPrimerResult(bases: String): BwaAlnInteractive.Result = if (!cacheResults) hitsByPrimer(bases) else {
        this.primerCache.getOrElse(bases, {
          val hit = hitsByPrimer(bases)
          this.primerCache(bases) = hit
          hit
        })
      }

      // Go through all the primer pairs that need to be mapped
      val progress = logger.map(ProgressLogger(_, noun="primer pairs", unit=10e5.toInt))
      val numPassing = primerPairsToMap.count { pp =>
        // Get the mappings for the left and right primer respectively
        val p1: BwaAlnInteractive.Result = getPrimerResult(pp.left.bases)
        val p2: BwaAlnInteractive.Result = getPrimerResult(pp.right.bases)
        // Get all possible amplicons from the left and right primer hits, filtering if there are too many for either
        val result = if (p1.hitCount > maxPrimerHits || p2.hitCount > maxPrimerHits) {
          PrimerOffTargetDetectorResult(primerPair=pp, passes=false, mappings=Nil)
        }
        else {
          val amplicons = Seq.newBuilder[Mapping]
          for (h1 <- p1.hits; h2 <- p2.hits) {
            if (h1.negative != h2.negative && h1.chrom == h2.chrom) {
              val (plus, minus) = if (h1.negative) (h2, h1) else (h1, h2)
              val contig   = dict(h1.chrom)
              val circular = contig.topology.contains(Topology.Circular)
              val primary = !(circular && contig.length <= plus.start && plus.end <= 2*contig.length)
              if (primary && minus.start > plus.end && CoordMath.getLength(plus.start, minus.end) <= maxAmpliconSize) {
                amplicons += Mapping(plus.chrom, plus.start, minus.end)
              }
            }
          }

          val amps = amplicons.result()

          PrimerOffTargetDetectorResult(
            primerPair = pp,
            passes     = amps.lengthCompare(maxPrimerPairHits) <= 0,
            mappings   = if (keepMappings) amps else Seq.empty,
            left       = if (KeepPrimerMappings) p1.hits.map(hitToMapping) else Seq.empty,
            right      = if (KeepPrimerMappings) p2.hits.map(hitToMapping) else Seq.empty
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

  private def hitToMapping(hit: BwaAlnInteractive.Hit): Mapping = new Mapping(refName = hit.chrom, start = hit.start, end = hit.end)

  /** Terminates the underlying BWA process. */
  override def close(): Unit = this.bwa.close()
}

