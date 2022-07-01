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

import com.fulcrumgenomics.FgBioDef.FgBioEnum
import com.fulcrumgenomics.fasta.SequenceDictionary
import enumeratum.EnumEntry
import htsjdk.samtools.util.{CoordMath, Locatable}

import java.lang.Math.{max, min}
import scala.collection.immutable


/** Trait that entries in [[Strand]] will extend. */
sealed trait Strand extends EnumEntry {
  /** The symbol for this strand. */
  def symbol: Char

  override def toString: String = String.valueOf(symbol)
}

/** Enum to represent the strand of a mapping to the genome. */
object Strand extends FgBioEnum[Strand] {
  case object Positive extends Strand { val symbol: Char = '+' }
  case object Negative extends Strand { val symbol: Char = '-' }

  override def values: immutable.IndexedSeq[Strand] = findValues

  def fromSymbol(symbol: Char): Strand = this.values.find(_.symbol == symbol).getOrElse {
    throw new IllegalArgumentException("Not a recognized strand symbol: " + symbol)
  }
}

object Mapping {
  import scala.math.Ordered.orderingToOrdered

  /** Builds a mapping from the string representation.
    *
    * A mapping should be one of the following two forms:
    * 1. `<chrom>:<start>-<end>`
    * 2. `<chrom>:<start>-<end>:<strand>`
    *
    * The first and second forms are produced by the [[Mapping.toString]] and [[Mapping.toStringWithStrand]] methods
    * respectively.
    * */
  def apply(line: String): Mapping = {
    val (chrom, range, strand) = line.split(':').toSeq match {
      case Seq(_chrom, _range, _strand) => (_chrom, _range, Strand.fromSymbol(_strand(0)))
      case Seq(_chrom, _range)          => (_chrom, _range, Strand.Positive)
      case _         => throw new IllegalArgumentException(s"Could not parse Mapping: $line")
    }
    val Seq(start, end) = range.split('-').toSeq.map(_.toInt)
    Mapping(chrom, start, end, strand)
  }

  /** Returns a tuple of reference index, start position, end position, and strand (0 forward, 1 reverse). */
  private def toTuple(mapping: Mapping, dict: SequenceDictionary): (Int, Int, Int, Int) = {
    val refIndex = dict(mapping.refName).index
    (refIndex, mapping.start, mapping.end, if (mapping.strand == Strand.Positive) 0 else 1)
  }

  /** Compares two mappings by genomic location (reference index, start, end, and strand [forward first]).*/
  def compare(a: Mapping, b: Mapping, dict: SequenceDictionary): Int = toTuple(a, dict) compare toTuple(b, dict)
}


/**
  * Represents the mapping of some item (a target, a primer, an amplicon) to a genome.
  *
  * All coordinates are 1-based closed end, i.e. the first 10bp on a sequence would be start=1, end=10.
  */
case class Mapping(refName: String, start: Int, end: Int, strand: Strand = Strand.Positive) extends Locatable with Ordered[Mapping] {
  assert(start >= 1, s"Start position must be >= 1. (start=$start end=$end)")
  assert(end >= start-1, s"End must be >= start-1. (start=$start end=$end)")

  import scala.math.Ordered.orderingToOrdered

  /** Compares mappings, returning which is earlier (by start, end, then strand). */
  override def compare(that: Mapping): Int = {
    require(this.refName == that.refName, "Cannot compare mappings on different chromosomes.")
    (this.start, this.end, if (this.strand == Strand.Positive) 0 else 1) compare (that.start, that.end, if (that.strand == Strand.Positive) 0 else 1)
  }

  /** The length or span of this mapping. */
  def length: Int = end - start + 1

  /** Returns a new mapping based on a relative start position and length within this mapping.
    * E.g. to capture the first ten bases of this mapping, use mapRelative(1, 10).
    *
    * @param start The 1-based start position of the new mapping relative to the start of this mapping
    * @param length The length of the new mapping
    * @param strand The strand of the new mapping; will be the same as source if omitted
    * @return a new mapping constructed with absolute coordinates, relative to this mapping
    */
  def resolve(start: Int, length: Int, strand:Strand = this.strand): Mapping = {
    assert(start >= 1, "Start of a relative submapping must be >= 1")
    assert(start <= this.length, "Start of a relative submapping must be <= source mapping length")
    assert(length >= 0, "Length of a submapping must be positive.")

    val absoluteStart = this.start + start - 1
    val absoluteEnd   = CoordMath.getEnd(absoluteStart, length)
    assert(absoluteEnd <= end, "End of submapping extends beyond end of source mapping.")

    copy(start=absoluteStart, end=absoluteEnd, strand=strand)
  }

  /** Projects a position into this mapping.  The position is required to be contained
    * within this mapping.  The returned position gives the coordinate relative to this mapping.
    *
    * E.g. Mapping(start=1000, end=2000).resolve(1000) => 1
    */
  def project(position: Int): Int= {
    assert(position >= this.start && position <= this.end, "Position is not within this mapping.")
    position - this.start + 1
  }

  /** Returns true if the two mappings overlap one another, false otherwise. */
  def overlaps(other: Mapping): Boolean = {
    this.refName == other.refName && this.start <= other.end && this.end >= other.start
  }

  /** Returns true if the other mapping is fully contained within this mapping */
  def contains(other: Mapping): Boolean = {
    this.refName == other.refName && this.start <= other.start && other.end <= this.end
  }

  /** Returns true if the two mappings abut one another, false otherwise. */
  def abuts(other: Mapping): Boolean = {
    require(this.refName == other.refName)
    this.start == other.end + 1 || this.end + 1 == other.start
  }

  /** Returns the length of the region which overlaps the other mapping, or zero if there is no overlap. */
  def lengthOfOverlapWith(other: Mapping): Int = {
    if (overlaps(other))
      CoordMath.getLength(max(this.start, other.start), min(this.end, other.end))
    else
      0
  }

  /** Returns a new [[Mapping]] that is the union of the two mappings.  They must overlap. */
  def union(other: Mapping): Mapping = {
    require(this.refName == other.refName)
    require(this.overlaps(other) || this.abuts(other)) // overlaps or abuts
    new Mapping(
      refName = this.refName,
      start   = math.min(this.start, other.start),
      end     = math.max(this.end,   other.end)
    )
  }

  /** Returns a new [[Mapping]] with the given amount added to the start and end position.*/
  def shift(amount: Int): Mapping = {
    require(this.start + amount >= 1)
    this.copy(start = start + amount, end = end + amount)
  }

  /** Returns the position of the 5' end of the mapping based on the strand. */
  def fivePrimePosition: Int = if (strand == Strand.Positive) start else end

  override def toString: String = {
    if (length == 1) refName + ":" + start else refName + ":" + start + "-" + end
  }

  /** Returns the string representation of this Mapping with the strand symbol at the end delimited by a colon. */
  def toStringWithStrand: String = this.toString + ":" + this.strand.symbol

  /** Returns a BED representation of this mapping.  This will return the first six columns of a detailed BED file. */
  def toBed: String = s"$refName\t${start-1}\t$end\t$toString,$length\t500\t${strand.symbol}"

  // Methods to make this [[Locatable]]
  override def getContig: String = this.refName
  override def getStart: Int = this.start
  override def getEnd: Int = this.end
}
