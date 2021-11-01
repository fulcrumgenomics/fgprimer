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

package com.fulcrumgenomics.primerdesign.primer3

import com.fulcrumgenomics.primerdesign.PrimerDesignDef.MinOptMax
import com.fulcrumgenomics.primerdesign.api.Mapping

import scala.collection.immutable.TreeMap


/** A class to contain various weights that Primer3 uses to weight penalties arising from
  * different characteristics of the primer(s).  For more information see:
  *     https://primer3.org/manual.html#globalTags
  */
case class Primer3Weights(productSizeLt: Double      = 1,
                          productSizeGt: Double      = 1,
                          productTmLt: Double        = 0.0,
                          productTmGt: Double        = 0.0,
                          primerEndStability: Double = 0.25,
                          primerGcLt: Double         = 0.25,
                          primerGcGt: Double         = 0.25,
                          primerSelfAny: Double      = 0.1,
                          primerSelfEnd: Double      = 0.1,
                          primerSizeLt: Double       = 0.5,
                          primerSizeGt: Double       = 0.1,
                          primerTmLt: Double         = 1.0,
                          primerTmGt: Double         = 1.0
                         ) {

  /** A map containing the design parameters to feed to Primer3 for the weights. */
  val tags: Map[Primer3InputTag, Any] = {
    import Primer3InputTag._
    TreeMap(
      PRIMER_PAIR_WT_PRODUCT_SIZE_GT -> productSizeLt,
      PRIMER_PAIR_WT_PRODUCT_SIZE_LT -> productSizeGt,
      PRIMER_PAIR_WT_PRODUCT_TM_LT   -> productTmLt,
      PRIMER_PAIR_WT_PRODUCT_TM_GT   -> productTmGt,
      PRIMER_WT_END_STABILITY        -> primerEndStability,
      PRIMER_WT_GC_PERCENT_LT        -> primerGcLt,
      PRIMER_WT_GC_PERCENT_GT        -> primerGcGt,
      PRIMER_WT_SELF_ANY             -> primerSelfAny,
      PRIMER_WT_SELF_END             -> primerSelfEnd,
      PRIMER_WT_SIZE_LT              -> primerSizeLt,
      PRIMER_WT_SIZE_GT              -> primerSizeGt,
      PRIMER_WT_TM_LT                -> primerTmLt,
      PRIMER_WT_TM_GT                -> primerTmGt
    )
  }
}

/**
  * Utility class for representing common options for [[Primer3Task]]s.
  *
  * @param ampliconSizes the min, optimal and max amplicon sizes
  * @param ampliconTms the min, optimal and max amplicon melting temperatures
  * @param primerSizes the min, optimal and max primer sizes
  * @param primerTms the min, optimal and max primer Tms
  * @param primerGcs the min an max acceptable GCs for individual primers (between 0 and 100)
  * @param gcClamp   the min and the max number of Gs/Cs in the 3' most 5 bases
  * @param primerMaxPolyx the max homopolymer length that is acceptable within a primer
  * @param primerMaxNs the maximum number of ambiguous bases acceptable within a primer
  * @param primerMaxDinucBases the maximum number of _bases_ in a dinucleotide run in a primer
  * @param avoidMaskedBases should Primer3 avoid placing primers in soft-masked sequence
  * @param howMany the number of primers to return
  */
case class Primer3Parameters(ampliconSizes: MinOptMax[Int],
                             ampliconTms:   MinOptMax[Double],
                             primerSizes: MinOptMax[Int],
                             primerTms: MinOptMax[Double],
                             primerGcs: MinOptMax[Int],
                             gcClamp: (Int, Int) = (0, 5),
                             primerMaxPolyx: Int = 5,
                             primerMaxNs: Int = 1,
                             primerMaxDinucBases: Int = 6,
                             avoidMaskedBases: Boolean = true,
                             howMany : Int = 5,
                            ) {
  assert(ampliconSizes.min <= ampliconSizes.opt, "Min amplicon size must be <= optimum amplicon size")
  assert(ampliconSizes.opt <= ampliconSizes.max, "Max amplicon size must be >= optimum amplicon size")
  if (ampliconTms.opt != 0) {
    assert(ampliconTms.min <= ampliconTms.opt, "Min amplicon Tm must be <= optimum amplicon Tm")
    assert(ampliconTms.opt <= ampliconTms.max, "Max amplicon Tm must be >= optimum amplicon Tm")
  }
  assert(primerSizes.min <= primerSizes.opt, "Min primer size must be <= optimum primer size")
  assert(primerSizes.opt <= primerSizes.max, "Max primer size must be >= optimum primer size")
  assert(primerTms.min <= primerTms.opt, "Min primer Tm must be <= optimum primer Tm")
  assert(primerTms.opt <= primerTms.max, "Max primer Tm must be >= optimum primer Tm")
  assert(primerGcs.min <= primerGcs.max, "Min primer GC must be <= max primer GC")
  assert(primerGcs.opt <= primerGcs.max, "Opt primer GC must be <= max primer GC")
  assert(gcClamp._1 <= gcClamp._2, "Min primer GC-clamp must be <= max primer GC-clamp")
  assert(primerMaxPolyx > 1, "Primer Max PolyX must be 2 or greater.")
  assert(primerMaxDinucBases >= 2,     "Primer Max Dinuc Bases must be 2 or greater")
  assert(primerMaxDinucBases % 2 == 0, "Primer Max Dinuc Bases must be an even number of bases")
  assert(primerMaxNs >= 0, "Primer Max Ns must be 0 or greater.")

  def maxAmpliconLength: Int = ampliconSizes.max

  def maxPrimerLength: Int = primerSizes.max

  /** Uses the internal state to generate a set of Primer3 input parameters that can be fed directly to primer3. */
  def tags: Map[Primer3InputTag, Any] = {
    import Primer3InputTag._
      TreeMap(
        PRIMER_NUM_RETURN              -> howMany,
        PRIMER_PRODUCT_SIZE_RANGE      -> s"${ampliconSizes.min}-${ampliconSizes.max}",
        PRIMER_PRODUCT_OPT_SIZE        -> ampliconSizes.opt,
        PRIMER_PRODUCT_MIN_TM          -> ampliconTms.min,
        PRIMER_PRODUCT_OPT_TM          -> ampliconTms.opt,
        PRIMER_PRODUCT_MAX_TM          -> ampliconTms.max,
        PRIMER_MIN_SIZE                -> primerSizes.min,
        PRIMER_OPT_SIZE                -> primerSizes.opt,
        PRIMER_MAX_SIZE                -> primerSizes.max,
        PRIMER_MIN_TM                  -> primerTms.min,
        PRIMER_OPT_TM                  -> primerTms.opt,
        PRIMER_MAX_TM                  -> primerTms.max,
        PRIMER_MIN_GC                  -> primerGcs.min,
        PRIMER_OPT_GC_PERCENT          -> primerGcs.opt,
        PRIMER_MAX_GC                  -> primerGcs.max,
        PRIMER_GC_CLAMP                -> gcClamp._1,
        PRIMER_MAX_END_GC              -> gcClamp._2,
        PRIMER_MAX_POLY_X              -> primerMaxPolyx,
        PRIMER_MAX_NS_ACCEPTED         -> primerMaxNs,
        PRIMER_LOWERCASE_MASKING       -> (if (avoidMaskedBases) 1 else 0)
      )
  }

  /** Generates a reasonably compact String representation of the parameters, for display/output. */
  override def toString: String = {
    Seq(
      "amp_sizes"     -> ampliconSizes.mkString("/"),
      "primer_sizes"  -> primerSizes.mkString("/"),
      "primer_tms"    -> primerTms.mkString("/"),
      "primer_gcs"    -> primerGcs.mkString("/"),
      "gc_clamp"      -> s"${gcClamp._1}/${gcClamp._2}",
      "polyx_max"     -> s"$primerMaxPolyx",
      "max_ns"        -> s"$primerMaxNs",
      "max_dinuc"     -> s"$primerMaxDinucBases",
      "avoid_repeats" -> (if (avoidMaskedBases) "yes" else "no")
    ).map(kv => kv._1 + "=" + kv._2).mkString(",")
  }
}

/**
  * Pulls together all the necessary inputs for running a primer or primer pair design
  * through Primer3.
  */
case class Primer3Input(target: Mapping,
                        task: Primer3Task,
                        params: Primer3Parameters,
                        weights: Primer3Weights = Primer3Weights()
                       ) {

  /** Generates the map of Primer3InputTag -> Value for the design of the target. */
  def tags(designRegion: Mapping): Map[Primer3InputTag, Any] = {
    task.inputTags(this.target, designRegion) ++ params.tags ++ weights.tags
  }
}
