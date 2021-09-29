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

package com.fulcrumgenomics.primerdesign.primer3

import com.fulcrumgenomics.primerdesign.testing.UnitSpec

/**
  * Simple tests for Primer3InputTag
  */
class Primer3InputTagTest extends UnitSpec {
  "Primer3InputTag" should "return the name as the result of toString" in {
    Primer3InputTag.values.foreach(tag  => tag.toString shouldBe tag.name())
  }

  it should "return true when asked if a valid input tag is an input tag" in {
    Primer3InputTag.includes(Primer3InputTag.SEQUENCE_ID.name()) shouldBe true
    Primer3InputTag.includes(Primer3InputTag.PRIMER_FIRST_BASE_INDEX.name()) shouldBe true
    Primer3InputTag.includes("PRIMER_EXPLAIN_FLAG") shouldBe true
  }

  it should "return false when asked if a string that is not an input tag is an input tag" in {
    Primer3InputTag.includes("Not a valid tag") shouldBe false
    Primer3InputTag.includes(Primer3InputTag.SEQUENCE_ID.name() + " ") shouldBe false
    Primer3InputTag.includes(Primer3InputTag.SEQUENCE_ID.name() + "=") shouldBe false
    Primer3InputTag.includes("PRIMER_ERROR") shouldBe false
  }
}
