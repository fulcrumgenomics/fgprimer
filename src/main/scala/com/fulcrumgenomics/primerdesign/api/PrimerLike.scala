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

import com.fulcrumgenomics.FgBioDef.FilePath
import com.fulcrumgenomics.bam.api.SamRecord
import com.fulcrumgenomics.commons.util.LazyLogging
import com.fulcrumgenomics.fasta.SequenceDictionary
import com.fulcrumgenomics.util.{Io, ProgressLogger}

trait PrimerLikeCompanion[A <: PrimerLike] extends LazyLogging {
  /** The value written to a string for missing (empty) bases. */
  val MissingBases: String = SamRecord.MissingBases

  /** Returns a function that defines "less-than" for two primer-like objects based on position. */
  def lt(dict: SequenceDictionary): (A, A) => Boolean = {
    (x: PrimerLike, y: PrimerLike) => {
      var retval = dict(x.mapping.refName).index - dict(y.mapping.refName).index
      if (retval == 0) retval = x.mapping.start - y.mapping.start
      if (retval == 0) retval = x.mapping.end   - y.mapping.end

      retval < 0
    }
  }

  /** Writes the given primers to the output in BED format. Returns the number of primers written. */
  def toBed(output: FilePath,
            primers: IterableOnce[A],
            description: Option[String] = None,
            writeProgress: Option[ProgressLogger] = Some(ProgressLogger(logger, noun="primers", verb="read", unit=1e6.toInt))
            ): Long = {
    var count  = 0L
    val header = s"""track name="${output.getFileName}" description="${description.getOrElse("User Track")}" """
    val lines  = primers.iterator.map { primer =>
      writeProgress.map(_.record(primer.mapping.refName, primer.mapping.start))
      count += 1
      primer.toBed12Row
    }

    Io.writeLines(output, Iterator(header) ++ lines)
    count
  }
}

/** The base trait for primers and primer pairs. */
trait PrimerLike {
  /** Returns the mapping of the primer or primer pair. */
  def mapping: Mapping

  /** Returns a string representing this object's primary mapping location. */
  def locationString: String = {
    val strand = if (mapping.strand == Strand.Positive) "F" else "R"
    s"${mapping.refName}_${mapping.start}_${mapping.end}_${strand}"
  }

  /** Formats the primer-like into 12 tab-separated fields matching the BED 12-column spec.
    * See: https://genome.ucsc.edu/FAQ/FAQformat.html#format1
    */
  def toBed12Row: String
}
