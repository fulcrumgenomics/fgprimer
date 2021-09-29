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
package com.fulcrumgenomics.primerdesign

import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.commons.util.LazyLogging
import com.fulcrumgenomics.primerdesign.api.Mapping
import com.fulcrumgenomics.vcf.api.{Allele, VcfSource, Variant => VcfVariant}
import htsjdk.samtools.util.OverlapDetector
import htsjdk.variant.variantcontext.VariantContext

import java.io.Closeable
import java.text.DecimalFormat


object VariantLookup {
  /** Yields an empty variant lookup. */
  val empty = new VariantOverlapDetector(Seq.empty, 0, false)

  /** Constructs a VariantLookup that caches all variants in memory for fast lookup.  Appropriate
    * only for small VCFs. */
  def cached(vcf: PathToVcf,
             minMaf: Double = 0.01,
             includingMissingMafs: Boolean = false,
             ): VariantLookup = new VariantOverlapDetector(Seq(vcf), minMaf, includingMissingMafs)

  /** Constructs a VariantLookup that queries indexed VCFs on disk for each lookup. */
  def diskBased(vcf: PathToVcf,
                minMaf: Double = 0.01,
                includingMissingMafs: Boolean = false
               ): VariantLookup = new FileBasedVariantLookup(Seq(vcf), minMaf, includingMissingMafs)
}

/** Trait that concrete implementations [[VariantLookup]]s should extend.  Allows the retrieval of variants
  * overlapping a given genomic coordinate range, as well as filtering variants by minor allele frequency. */
abstract class VariantLookup(private val minMaf: Double, includingMissingMaf: Boolean) extends Closeable {
  /** Gets all variants overlapping a given genomic coordinate range that are not filtered.
    *
    * @param chrom the chromosome
    * @param start the 1-based start position
    * @param end the 1-based end position
    * @param maf optionally, return only variants with at least this minor allele frequency
    * @param includeMissingMafs when filtering variants with a minor allele frequency, true to include variants with no
    *                           annotated minor allele frequency, otherwise false.  If no minor allele frequency is
    *                           given to this method then this parameter does nothing.
    * @return
    */
  final def query(chrom: String,
                  start:Int,
                  end: Int,
                  maf: Double = this.minMaf,
                  includeMissingMafs: Boolean = this.includingMissingMaf
                 ): Seq[Variant] = {
    val variants = this._query(chrom, start, end)
    if (maf <= 0) variants else {
      if (includeMissingMafs) variants.filter { variant => variant.maf.forall(_ >= maf) }
      else                    variants.filter { variant => variant.maf.exists(_ >= maf) }
    }
  }

  /** Converts variants from a [[VariantContext]] to a [[Variant]], filters variants based on their FILTER status, and
    * sorts by start position. */
  protected def toVariants(variants: Iterable[VcfVariant]): Seq[Variant] = {
    variants
      .filter(v => v.filters.isEmpty || v.filters == VcfVariant.PassingFilters)
      .map(v => Variant.apply(v))
      .toSeq
      .sortWith((a, b) => a.pos <= b.pos)
  }

  /** All classes should implement this method. */
  protected def _query(chrom: String, start:Int, end: Int): Seq[Variant]
}

class FileBasedVariantLookup(private val vcfs: Seq[PathToVcf], minMaf: Double, includeMissingMafs: Boolean)
  extends VariantLookup(minMaf, includeMissingMafs) with Closeable with LazyLogging {
  private val readers = vcfs.map(vcf => VcfSource(vcf)).toIndexedSeq

  /**
    * Queries variants from the VCFs used by this lookup, sorts the results and returns them as a List
    */
  protected def _query(chrom: String, start:Int, end: Int): Seq[Variant] = {
    // Log some warnings if we get queries for chromosome names we don't support
    readers.zip(vcfs).foreach { case (reader, vcf) =>
      if (reader.header.dict.get(chrom).isEmpty)
        logger.warning(s"VCF file ${vcf} does not appear to contain data for chromosome ${chrom}.")
    }

    // Query each of the VCFs and merge the results
    toVariants(readers.flatMap(r => r.query(chrom, start, end)))
  }

  /** Closes all the VCF readers. */
  override def close(): Unit = this.readers.foreach(_.close())
}


/**
  * Implements VariantLookup by reading the entire VCFs into memory and loading the resulting
  * Variants into an OverlapDetector.
  */
class VariantOverlapDetector(private val vcfs: Seq[PathToVcf], minMaf: Double, includeMissingMafs: Boolean)
  extends VariantLookup(minMaf, includeMissingMafs) with LazyLogging {
  private val detector: OverlapDetector[Variant] = {
    val d = new OverlapDetector[Variant](0, 0)
    vcfs.foreach { vcf =>
      val reader = VcfSource(vcf)
      toVariants(reader).foreach { v => d.addLhs(v, v.toMapping)}
      reader.close()
    }
    d
  }

  protected def _query(chrom: String, start:Int, end: Int): Seq[Variant] = {
    this.detector.getOverlaps(new Mapping(chrom, start, end)).toSeq
  }

  override def close(): Unit = ()
}


/**
  * Companion object for the Variant class with various useful helper methods.
  */
object Variant {
  private val mafFormatter = new DecimalFormat("0.000")

  /** Constructs a simplified Variant object from a VCF Variant. */
  def apply(variant: VcfVariant): Variant = {
    // Inner function to determine the minor allele frequency of a VariantContext
    val maf = {
      if (variant.attrs.contains("CAF")) {
        Some(1 - variant[IndexedSeq[String]]("CAF").head.toDouble)
      }
      else if (variant.attrs.contains("AF")) {
        Some(variant[IndexedSeq[Float]]("AF").sum.toDouble)
      }
      else if (variant.attrs.contains("AC") && variant.attrs.contains("AN")) {
        Some(variant[IndexedSeq[Int]]("AC").sum / variant[Int]("AN").toDouble)
      }
      else if (variant.genotypes.nonEmpty) {
        val calls = variant.gts.flatMap(g => g.calls).filter(_ != Allele.NoCallAllele).toIndexedSeq
        Some(calls.count(_ != variant.alleles.ref) / calls.size.toDouble)
      }
      else {
         None
      }
    }

    new Variant(
      id    = variant.id.getOrElse(s"${variant.chrom}:${variant.pos}"),
      chrom = variant.chrom,
      pos   = variant.pos,
      ref   = variant.alleles.ref.value,
      alt   = variant.alleles.alts.head.value,
      maf   = maf
    )
  }
}

/** The type of variant for [[Variant]]. */
object VariantType extends Enumeration {
  val Snp, Insertion, Deletion, Other = Value
}

/**
  * Simple class to represent Variants during primer design that is massively simpler and easier to
  * work with than VariantContext objects that are generated from VCFs.
  *
  * @param id the id or name of the variant
  * @param chrom the chromosome on which the variant resides
  * @param pos the position of the variant (or in the case of indels, the base to the left of the variant)
  * @param ref the reference allele
  * @param alt the single alternative allele
  * @param maf the minor allele frequency if available
  */
case class Variant(id: String, chrom: String, pos: Int, ref: String, alt: String, maf:Option[Double]) {

  /** Returns the type of variant represented by this object. */
  val variantType: VariantType.Value = {
    if (ref.length == 1 && alt.length == 1) VariantType.Snp
    else if (ref.length == 1 && alt.length > 1) VariantType.Insertion
    else if (ref.length >  1 && alt.length == 1) VariantType.Deletion
    else VariantType.Other
  }

  /**
    * Creates a compact String representation of the variant that includes all relevant information.
    */
  override def toString: String = {
    val mafString = maf.map(Variant.mafFormatter.format).getOrElse("na")
    s"${id}@${chrom}:${pos}[${ref}/${alt} ${mafString}]"
  }

  /** Creates a [[Mapping]] representing the genomic span of this variant. Insertions will only span the base
    * preceding the inserted bases. */
  def toMapping: Mapping = {
    val end =  variantType match {
      case VariantType.Deletion => this.pos + this.ref.length - 1
      case _                    => this.pos
    }
    new Mapping(refName = this.chrom, start = this.pos, end = end)
  }
}
