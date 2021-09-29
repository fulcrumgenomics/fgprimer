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

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.commons.util.LazyLogging
import com.fulcrumgenomics.fasta.SequenceDictionary
import com.fulcrumgenomics.primerdesign.primer3.Primer3Parameters
import com.fulcrumgenomics.util.Sequences
import htsjdk.samtools.util.SequenceUtil

object Primer extends PrimerLikeCompanion[Primer] {
  /** Compares to [[Primer]]s by their [[Mapping]]. */
  def compare(a: Primer, b: Primer, dict: SequenceDictionary): Int = Mapping.compare(a.mapping, b.mapping, dict)
}


/** Represents an individual primer.
  *
  * @param bases the bases of the primer (excluding any tail)
  * @param tm the calculated melting temperature of the primer
  * @param penalty the penalty or score for the primer
  * @param mapping the mapping of the primer to the genome
  * @param namePrefix an optional prefix to use when dynamically generating a name/id for the primer
  * @param name an optional name to use for the primer - is mutually exclusive with namePrefix
  * @param tail an optional tail sequence to put on the 5' end of the primer
  * @param designParameters the design parameters that were used to generate the primer
  */
case class Primer(bases: String,
                  tm: Double,
                  penalty: Double,
                  mapping: Mapping,
                  namePrefix: Option[String] = None,
                  name: Option[String] = None,
                  tail: Option[String] = None,
                  designParameters: Option[Primer3Parameters] = None) extends PrimerLike {
  require(mapping.length == bases.length || bases.isEmpty, s"bases: ${bases.length} mapping: ${mapping.length} $mapping")
  require(namePrefix.size + name.size <= 1, "Cannot define both name and namePrefix.")

  /** The GC of the primer in the range 0-100, or zero if the bases are empty. */
  val gc: Double = if (bases.isEmpty) 0 else SequenceUtil.calculateGc(bases.getBytes) * 100

  /** The length of the un-tailed primer. */
  val length: Int = mapping.length

  /** The length of the longest homopolymer in the primer. */
  val longestHomopolymer: Int = Sequences.longestHomopolymer(bases).length

  /** Returns an ID for this primer, using the provided prefix, that is likely to be unique. */
  def id: String = (this.name, this.namePrefix) match {
    case (Some(n), _      ) => n
    case (_,       Some(p)) => s"${p}_${locationString}"
    case _                  => locationString
  }

  override def toString: String = {
    val bases = if (this.bases.isEmpty) Primer.MissingBases else this.bases
    s"$bases\t$tm\t$penalty\t${mapping.toStringWithStrand}"
  }

  /** Returns the BED detail format view: https://genome.ucsc.edu/FAQ/FAQformat.html#format1.7 */
  override def toBed12Row: String = {
    Seq(
      mapping.refName,       // contig
      mapping.start-1,       // start
      mapping.end,           // end
      this.id,               // name
      500,                   // score
      mapping.strand.symbol, // strand
      mapping.start-1,       // thick start
      mapping.end,           // thick end
      "100,100,100",         // color
      1,                     // block count
      s"${mapping.length}",  // block sizes
      s"0"                   // block starts (relative to `start`)
    ).map(_.toString).mkString("\t")
  }

  /** Returns a copy of the primer with the tail sequence attached. */
  def withTail(tail: String): Primer = copy(tail=Some(tail))

  /** Returns a copy of ths primer with the given name prefix. */
  def withNamePrefix(prefix: String): Primer = copy(namePrefix=Some(prefix), name=None)

  /** Returns a copy of ths primer with the given name. */
  def withName(name: String): Primer = copy(namePrefix=None, name=Some(name))

  /** Returns the sequence of the primer including the tail. */
  def basesWithTail: String = this.tail.getOrElse("") + this.bases
}
