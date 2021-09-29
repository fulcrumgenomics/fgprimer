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

import com.fulcrumgenomics.primerdesign.api.Mapping
import com.fulcrumgenomics.primerdesign.primer3.Primer3InputTag._
import htsjdk.samtools.util.CoordMath

/** Trait for the various types of primer designs we can perform. */
sealed trait Primer3Task {
  /** The tag returned by Primer3 that gives the number of primers returned. */
  private[primer3] def countTag: String

  /** None if designing primer pairs, or `"LEFT"` or `"RIGHT"` depending on which design type is being performed. */
  private[primerdesign] def side: Option[String] = None

  /** The set of input parameters specific to this design type. */
  private[primer3] def inputTags(target: Mapping, designRegion: Mapping): Map[Primer3InputTag, Any]
}

/** Picks both a left and right primer (pair). */
case object DesignPrimerPairsTask extends Primer3Task {
  override def inputTags(target: Mapping, designRegion: Mapping): Map[Primer3InputTag, Any] = Map(
    PRIMER_TASK                -> "generic",
    PRIMER_PICK_LEFT_PRIMER    -> 1,
    PRIMER_PICK_RIGHT_PRIMER   -> 1,
    PRIMER_PICK_INTERNAL_OLIGO -> 0,
    SEQUENCE_TARGET            -> s"${target.start - designRegion.start + 1},${target.length}"
  )

  override val countTag: String = "PRIMER_PAIR_NUM_RETURNED"
}

/** Picks only the left primer */
case object DesignLeftPrimersTask extends Primer3Task {
  override def inputTags(target: Mapping, designRegion: Mapping): Map[Primer3InputTag, Any] = Map(
    PRIMER_TASK                -> "pick_primer_list",
    PRIMER_PICK_LEFT_PRIMER    -> 1,
    PRIMER_PICK_RIGHT_PRIMER   -> 0,
    PRIMER_PICK_INTERNAL_OLIGO -> 0,
    SEQUENCE_INCLUDED_REGION   -> s"1,${target.start - designRegion.start}"
  )
  override val countTag: String = "PRIMER_LEFT_NUM_RETURNED"
  override val side: Option[String] = Some("LEFT")
}

/** Picks only the right primer. */
case object DesignRightPrimersTask extends Primer3Task {
  override def inputTags(target: Mapping, designRegion: Mapping): Map[Primer3InputTag, Any] = {
    val start  = target.end - designRegion.start + 1
    val length = CoordMath.getLength(target.end + 1, designRegion.end)
    Map(
      PRIMER_TASK                -> "pick_primer_list",
      PRIMER_PICK_LEFT_PRIMER    -> 0,
      PRIMER_PICK_RIGHT_PRIMER   -> 1,
      PRIMER_PICK_INTERNAL_OLIGO -> 0,
      SEQUENCE_INCLUDED_REGION   -> s"$start,$length"
    )
  }
  override val countTag: String = "PRIMER_RIGHT_NUM_RETURNED"
  override val side: Option[String] = Some("RIGHT")
}
