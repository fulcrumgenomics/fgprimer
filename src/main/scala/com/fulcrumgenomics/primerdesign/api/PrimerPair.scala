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

object PrimerPair extends PrimerLikeCompanion[PrimerPair] {
  /** Compares to primer pairs by their amplicon's [[Mapping]].
    * @param a the first primer pair to compare
    * @param b the second primer pair to compare
    * @param dict the sequence dictionary
    * @param byAmplicon true to compare based on the [[Mapping]] of the amplicon of the primer pair, otherwise use the
    *                   [[Mapping]] of the leftPrimerMappings and rightPrimerMappings primers
    * @return
    */
  def compare(a: PrimerPair, b: PrimerPair, dict: SequenceDictionary, byAmplicon: Boolean = true): Int = {
    if (byAmplicon) {
      Mapping.compare(a.amplicon, b.amplicon, dict)
    }
    else {
      var ret = Mapping.compare(a.left.mapping, b.left.mapping, dict)
      if (ret == 0) ret = Mapping.compare(a.right.mapping, b.right.mapping, dict)
      ret
    }
  }
}

/**
  * Represents a pair of primers that work together to create an amplicon.
  * @param left  an object describing the leftPrimerMappings primer
  * @param right an object describing the rightPrimerMappings primer
  * @param amplicon a mapping to the genome that represents the region expected to be amplified
  * @param ampliconSequence the sequence of the expected amplicon
  * @param tm the melting temperature of the expected amplicon
  * @param penalty the penalty value assigned by primer3 to the primer pair
  * @param designParameters optionally the set of design parameters used to drive primer3 to create this primer pair
  */
case class PrimerPair(left: Primer,
                      right: Primer,
                      amplicon: Mapping,
                      ampliconSequence: String,
                      tm: Double,
                      penalty: Double,
                      namePrefix: Option[String] = None,
                      name: Option[String] = None,
                      designParameters: Option[Primer3Parameters] = None) extends PrimerLike {
  require(ampliconSequence.isEmpty || amplicon.length == ampliconSequence.length, s"amplicon: ${amplicon.length} ampliconSequence: ${ampliconSequence.length}")
  require(amplicon.length == right.mapping.end - left.mapping.start + 1, "amplicon length differs")
  require(left.mapping.refName == right.mapping.refName && left.mapping.refName == amplicon.refName, "references differ")
  require(namePrefix.size + name.size <= 1, "Only one of namePrefix or name may be provided.")

  /** Returns the mapping for the amplicon. */
  override def mapping: Mapping = this.amplicon

  /** The list of individual primers (leftPrimerMappings and rightPrimerMappings) for the pair. */
  def primers: Iterator[Primer] = Iterator(left, right)

  /** Returns the length of the amplicon. */
  def length: Int = amplicon.length

  /** The inner region of the amplicon (not including the primers). I.e. the region of the genome covered by the
    * primer pair, without the primer regions.  If the primers overlap, then the inner mapping is the midpoint at where
    * they overlap*/
  def inner: Mapping = {
    if (left.mapping.overlaps(right.mapping)) {
      val midpoint = (left.mapping.end + right.mapping.start) / 2
      left.mapping.copy(start = midpoint, end = midpoint)
    }
    else {
      left.mapping.copy(start = left.mapping.end + 1, end = right.mapping.start - 1)
    }
  }

  /** The GC of the amplicon sequence in the range 0-100, or zero if there is no amplicon sequence. */
  val gc: Double = if (this.ampliconSequence.isEmpty) 0 else Sequences.gcContent(this.ampliconSequence) * 100

  /** Returns an ID for this primer pair, using a) the optional name if present, b) the optional prefix if present or
    * c) just location information. */
  def id: String = (this.name, this.namePrefix) match {
    case (Some(n), _) => n
    case (_, Some(p)) => s"${p}_${locationString}"
    case _            => locationString
  }

  /** Returns a copy of the primer pair where the leftPrimerMappings and rightPrimerMappings primer are tailed.  If either tail is the empty
    * string it will not be added to the respective primer.
    */
  def withTails(leftTail: String, rightTail: String): PrimerPair = {
    val left  = if (leftTail.nonEmpty) this.left.withTail(leftTail) else this.left
    val right = if (rightTail.nonEmpty) this.right.withTail(rightTail) else this.right
    copy(left=left, right=right)
  }

  /** Returns a copy of the primer pair that has the given name prefix on the pair and on the primers within. */
  def withNamePrefix(prefix: String): PrimerPair = copy(
    name       = None,
    namePrefix = Some(prefix),
    left       = left.withNamePrefix(prefix),
    right      = right.withNamePrefix(prefix),
  )

  /** Returns a copy of the primer pair with names assigned to the primer pair, leftPrimerMappings and rightPrimerMappings primers.
    *
    * @param ppName the name for the primer pair
    * @param lpName the name for the leftPrimerMappings primer
    * @param rpName the name for the rightPrimerMappings primer
    * */
  def withNames(ppName: String, lpName: String, rpName: String): PrimerPair = copy(
    name = Some(ppName),
    namePrefix = None,
    left       = left.withName(lpName),
    right      = right.withName(rpName)
  )

  override def toBed12Row: String = {
    val start = amplicon.start

    Seq(left.mapping.refName,                            // chromosome
      left.mapping.start-1,                              // start
      right.mapping.end,                                 // end
      id,                                                // name
      500,                                               // score
      mapping.strand.symbol,                             // strand
      left.mapping.start-1,                              // thick start
      right.mapping.end,                                 // thick end
      "100,100,100",                                     // color
      3,                                                 // block count
      s"${left.mapping.length},${inner.length},${right.mapping.length}", // block lengths
      s"0,${inner.start-start},${right.mapping.start-start}" // block starts (relative to `start`)
    ).map(_.toString).mkString("\t")
  }

  /**
    * Returns a string with information about the leftPrimerMappings and rightPrimerMappings primers respectively, and the primer pair itself.
    */
  override def toString: String = {
    val builder = new StringBuilder
    val ampliconSequence = if (this.ampliconSequence.isEmpty) PrimerPair.MissingBases else this.ampliconSequence
    builder.append(this.left.toString)
    builder.append('\t')
    builder.append(this.right.toString)
    builder.append('\t')
    builder.append(s"$ampliconSequence\t${this.tm}\t${this.penalty}")
    builder.toString
  }
}
