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

package com.fulcrumgenomics.primerdesign.testing


import com.fulcrumgenomics.FgBioDef._
import com.fulcrumgenomics.bam.api.{SamRecord, SamSource}
import com.fulcrumgenomics.commons.io.PathUtil
import com.fulcrumgenomics.commons.util.{LogLevel, Logger}
import com.fulcrumgenomics.fasta.Converters._
import com.fulcrumgenomics.fasta.{SequenceDictionary, SequenceMetadata}
import com.fulcrumgenomics.primerdesign.api.{Mapping, Primer, PrimerLike, PrimerPair}
import htsjdk.samtools.{SAMFileHeader, SAMReadGroupRecord}
import org.scalatest.{BeforeAndAfterAll, OptionValues}
import org.scalatest.flatspec.AnyFlatSpec
import org.scalatest.matchers.should.Matchers

import java.nio.file.{Files, Path}

/** Base class for unit and integration testing */
trait UnitSpec extends AnyFlatSpec with Matchers with OptionValues {
  // Turn down HTSJDK logging
  htsjdk.samtools.util.Log.setGlobalLogLevel(htsjdk.samtools.util.Log.LogLevel.WARNING)

  /** Gets the resource at the given path and returns a path to it. */
  protected def getResourceAsPath(path: String): Path = {
    Seq(getClass.getResource _, getClass.getClassLoader.getResource _)
      .flatMap(m => Option(m(path)))
      .headOption match {
      case None    => throw new IllegalArgumentException(s"Resource does not exist at path: $path")
      case Some(in) => PathUtil.pathTo(in.getPath)
    }
  }

  /** Gets the resource with the given name and associated package, and returns a path to it. */
  protected  def testResource(name: String, pkg: Package = getClass.getPackage): Path = {
    val packagePath = pkg.getName.replace('.', '/')
    val url = getClass.getClassLoader.getResource(packagePath + "/" + name)
    PathUtil.pathTo(url.getPath)
  }

  /** Creates a new temp file for use in testing that will be deleted when the VM exits. */
  protected def makeTempFile(prefix: String, suffix: String) : Path = {
    val path = Files.createTempFile(prefix, suffix)
    path.toFile.deleteOnExit()
    path
  }

  /** Returns a basic sequence dictionary containing chr1-chr22, chrX, chrY and chrM. */
  protected def makeSequenceDictionary: SequenceDictionary =
    SequenceDictionary(
      (Range(1, 23) ++ Seq("X", "Y", "M"))
      .zipWithIndex
      .map { case (chr, idx) => SequenceMetadata("chr" + chr, 200e6.toInt, index=idx) }
    )

  /** Makes a basic SAMFileHeader for testing. */
  protected def makeSamHeader: SAMFileHeader = {
    val header = new SAMFileHeader(makeSequenceDictionary.asSam)
    val rg = new SAMReadGroupRecord("SomeId")
    rg.setSample("SomeSample")
    rg.setLibrary("SomeLibrary")
    header.addReadGroup(rg)
    header
  }

  /** Reads all the records from a SAM or BAM file into an indexed seq. */
  protected def readBamRecs(bam: PathToBam): IndexedSeq[SamRecord] = SamSource(bam).toIndexedSeq


  private def comparePrimers(actualPrimer: Primer, expectedPrimer: Primer, dict: SequenceDictionary): Unit = {
    Primer.compare(actualPrimer, expectedPrimer, dict) shouldBe 0
    actualPrimer.bases                                 shouldBe expectedPrimer.bases
    actualPrimer.tm                                    shouldBe expectedPrimer.tm +- 0.01
    actualPrimer.penalty                               shouldBe expectedPrimer.penalty +- 0.01
  }

  private def comparePrimerPairs(actualPair: PrimerPair, expectedPair: PrimerPair, dict: SequenceDictionary, ampliconSequenceIsEmpty: Boolean = true): Unit = {
    compare(actualPair.left, expectedPair.left, dict)
    compare(actualPair.right, expectedPair.right, dict)

    Mapping.compare(actualPair.amplicon, expectedPair.amplicon, dict)

    if (ampliconSequenceIsEmpty) {
      actualPair.ampliconSequence shouldBe empty
    }
    else {
      actualPair.ampliconSequence shouldBe expectedPair.ampliconSequence
    }

    actualPair.tm      shouldBe expectedPair.tm +- 0.01
    actualPair.penalty shouldBe expectedPair.penalty +- 0.01
  }

  protected def compare(actualPrimer: PrimerLike, expectedPrimer: PrimerLike, dict: SequenceDictionary): Unit = {
    (actualPrimer, expectedPrimer) match {
      case (actual: Primer,     expected: Primer)     => this.comparePrimers(actual, expected, dict)
      case (actual: PrimerPair, expected: PrimerPair) => this.comparePrimerPairs(actual, expected, dict, expected.ampliconSequence.isEmpty)
      case _                                          => unreachable("Comparing different types.")
    }
  }

  protected def compare[T <: PrimerLike](actualPrimerPairs: Seq[T], expectedPrimerPairs: Seq[T], dict: SequenceDictionary): Unit = {
    actualPrimerPairs.length shouldBe expectedPrimerPairs.length
    actualPrimerPairs.zip(expectedPrimerPairs).foreach { case (actual, expected) => compare(actual, expected, dict) }
  }
}

/** Base class that turns up logging to [[LogLevel.Error]] before all the tests and restores
  * the log level after all the tests.
  */
trait ErrorLogLevel extends UnitSpec with BeforeAndAfterAll {
  private var logLevel = Logger.level

  override protected def beforeAll(): Unit = {
    this.logLevel = Logger.level
    Logger.level  = LogLevel.Error
  }

  override protected def afterAll(): Unit = {
    Logger.level = LogLevel.Info
    Logger.level = this.logLevel
  }
}
