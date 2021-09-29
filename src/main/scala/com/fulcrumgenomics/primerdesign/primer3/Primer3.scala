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

import com.fulcrumgenomics.FgBioDef.{FilePath, PathToFasta, unreachable}
import com.fulcrumgenomics.commons.util.LazyLogging
import com.fulcrumgenomics.fasta.SequenceDictionary
import com.fulcrumgenomics.primerdesign.api.{Mapping, Primer, PrimerPair, Strand}
import com.fulcrumgenomics.primerdesign.{PrimerDesignDef, VariantLookup, VariantType, api}
import com.fulcrumgenomics.util.{Io, Sequences}
import htsjdk.samtools.Defaults
import htsjdk.samtools.reference.ReferenceSequenceFileFactory
import htsjdk.samtools.util.CoordMath

import java.io._
import java.nio.file.Paths
import scala.collection.immutable.TreeMap
import scala.collection.mutable
import scala.collection.mutable.ListBuffer

/** Trivial class to encapsulate how many primers/pairs failed for a given reason. */
case class Primer3Failure(reason: Primer3FailureReason, count: Int)


/**
  * Case class that encapsulates the results of running primer3 to design PCR primer
  * pairs. Contains the ordered (by objective function score) list of primer pairs
  * that were returned, along with an ordered sequence of Primer3Failures explaining
  * how many primers/pairs were eliminated from consideration and for what reasons.
  */
case class Primer3Result(results: Either[Seq[Primer], Seq[PrimerPair]], failures: Seq[Primer3Failure]) {
  def primers:     Seq[Primer]     = results.left.getOrElse(throw new IllegalStateException("Result contained pairs, not primers"))
  def primerPairs: Seq[PrimerPair] = results.right.getOrElse(throw new IllegalStateException("Result contained primers, not pairs"))
}


/** Companion to the Primer3 class. */
object Primer3 {
  /** The default value that will be used for the path to Primer3. */
  val DefaultPrimer3Executable: FilePath = Paths.get("primer3_core")
}


/** A class to interact with Primer3. */
class Primer3(val executable: FilePath = Primer3.DefaultPrimer3Executable,
              val genome: PathToFasta,
              private val variantLookup: VariantLookup = VariantLookup.empty,
             ) extends Primer3Tool with LazyLogging with Closeable {

  /** The Primer3 process, and input and output streams. */
  private val process  = new ProcessBuilder(this.executable.toString, "-strict_tags").redirectErrorStream(true).start()
  private val p3Input  = new PrintStream(new BufferedOutputStream(this.process.getOutputStream, Io.bufferSize), true)
  private val p3Output = new BufferedReader(new InputStreamReader(this.process.getInputStream), Defaults.BUFFER_SIZE)

  /** The reference genome sequence */
  private val ref  = ReferenceSequenceFileFactory.getReferenceSequenceFile(genome)
  private val dict = SequenceDictionary.extract(genome)

  /** A set of parameters that are always passed to Primer3. */
  private val globalParameters = TreeMap(
    Primer3InputTag.PRIMER_FIRST_BASE_INDEX  -> 1,
    Primer3InputTag.PRIMER_EXPLAIN_FLAG      -> 1
  )

  /** Designs primers (or primer pairs) given a target region. */
  def designPrimers(in: Primer3Input): Primer3Result = {
    val padding = in.params.maxAmpliconLength - in.target.length
    val region  = in.target.copy(
      start = math.max(1, in.target.start-padding),
      end   = math.min(in.target.end + padding, this.dict(in.target.refName).length)
    )

    val (softMasked, hardMasked) = designSequences(region)
    val tags = this.globalParameters ++ in.tags(region) + (Primer3InputTag.SEQUENCE_TEMPLATE -> hardMasked)

    // Send the command and read back the output
    tags.foreach { case (tag, value) => p3Input.println(s"${tag.name()}=$value") }
    p3Input.println("=")

    // Read back the output, and fail if either (1) anything gets written to stderr, or (2) PRIMER_ERROR was in the output
    val output     = mutable.Map[String,String]()
    val errorLines = mutable.ListBuffer[String]()

    def primer3Error(message: String): Nothing = {
      output.get("PRIMER_ERROR") match {
        case Some(primer3Error) => fail(message, s"Primer3 returned an error: $primer3Error" +: errorLines.toSeq)
        case None               => fail(message, errorLines.toSeq)
      }
    }

    var line : String = p3Output.readLine()
    while (line != "=") {
      line match {
        case null                  => primer3Error("primer3 exited prematurely")
        case ""                    => () // ignore emtpy lines
        case e if !e.contains("=") => errorLines += e
        case l                     =>
            val (key, value) = PrimerDesignDef.splitInTwo(line, '=')
            if (!Primer3InputTag.includes(key)) output(key) = value
      }

      line = p3Output.readLine()
    }

    if (errorLines.nonEmpty) primer3Error("primer3 failed")

    in.task match {
      case DesignPrimerPairsTask =>
        // Parse out the results and do some post filtering for dinucs
        val all = buildPrimerPairs(in=in, output=output.toMap, designRegion=region, unmaskedDesignSequence=softMasked)
        val (good, dinucFailurePairs) = all.partition(_.primers.forall(p => Sequences.longestDinuc(p.bases).length <= in.params.primerMaxDinucBases))
        val dinucFailurePrimers = dinucFailurePairs.flatMap(_.primers).filterNot(p => Sequences.longestDinuc(p.bases).length <= in.params.primerMaxDinucBases)
        val failures = buildFailures(in.params, dinucFailurePrimers)(output("PRIMER_LEFT_EXPLAIN"), output("PRIMER_RIGHT_EXPLAIN"), output("PRIMER_PAIR_EXPLAIN"))
        Primer3Result(results=Right(good), failures)
      case task =>
        // Parse out the results and do some post filtering for dinucs
        val side = task.side.getOrElse(unreachable("Non-paired design must have a defined side."))
        val all = buildPrimers(in=in, output=output.toMap, designRegion=region, unmaskedDesignSequence=softMasked, side=side)
        val (good, dinucFailures) = all.partition(p => Sequences.longestDinuc(p.bases).length <= in.params.primerMaxDinucBases)
        val failures = buildFailures(in.params, dinucFailures)(output(s"PRIMER_${side}_EXPLAIN"))
        Primer3Result(results=Left(good), failures)
    }
  }

  /**
    * Returns two sequences, as Strings.  The first sequence is all upper-case as is from the underlying FASTA file (i.e.
    * keeps soft-masking).  The second is the same as the first, but with common variation hard-masked.
    *
    * @param region The region of the genome to be extracted
    * @return a tuple of two sequences: the sequence for the region, and the sequence for the region with common variants
    *         hard-masked (as Ns).
    */
  protected[primerdesign] def designSequences(region: Mapping): (String,String) = {
    val softMasked = new String(ref.getSubsequenceAt(region.refName, region.start, region.end).getBases)

    // Accumulate a set of positions to mask
    val positions = ListBuffer[Int]()

    this.variantLookup.query(region.refName, region.start, region.end).foreach(v => {
      v.variantType match {
        case VariantType.Snp       => positions.append(v.pos)
        case VariantType.Insertion => positions.appendAll(Iterator(v.pos, v.pos+1)) // Add one base either side of insertions
        case VariantType.Deletion  => (1 until v.ref.length).foreach(offset => positions.append(v.pos+ offset)) // Add all the deleted bases
        case VariantType.Other     => (0 to v.ref.length).foreach(offset => positions.append(v.pos + offset))
      }
    })

    // Mask the positions
    val bytes = softMasked.getBytes()
    positions.filter(x => x >= region.start && x <= region.end).foreach(pos => bytes(region.project(pos)-1) = 'N')
    (softMasked, new String(bytes))
  }

  /** Builds an ordered sequence of primer pairs based on the order they were returned by Primer3,
    * using a map of key=value that was constructed from the Primer3 output. If no primer pairs
    * were returned, will return an empty Seq.
    */
  private def buildPrimerPairs(in: Primer3Input,
                               output: Map[String,String],
                               designRegion: Mapping,
                               unmaskedDesignSequence: String) : Seq[PrimerPair] = {
    val designParametersOption = Some(in.params)

    val lefts  = buildPrimers(in, output, "LEFT",  designRegion, unmaskedDesignSequence)
    val rights = buildPrimers(in, output, "RIGHT", designRegion, unmaskedDesignSequence)

    lefts.zip(rights).zipWithIndex map { case((left, right), num) =>
      val amplicon = left.mapping.copy(end=right.mapping.end)
      api.PrimerPair(
        left=left,
        right=right,
        tm=output(s"PRIMER_PAIR_${num}_PRODUCT_TM").toDouble,
        penalty=output(s"PRIMER_PAIR_${num}_PENALTY").toDouble,
        ampliconSequence=unmaskedDesignSequence.substring(designRegion.project(amplicon.start)-1, designRegion.project(amplicon.end)),
        amplicon=amplicon,
        designParameters=designParametersOption
      )
    }
  }

  /** Builds a Seq of primers for either the LEFT or RIGHT primers from the Primer3 outputs. */
  private[primerdesign] def buildPrimers(in: Primer3Input,
                                         output: Map[String,String],
                                         side: String,
                                         designRegion: Mapping,
                                         unmaskedDesignSequence: String) : Seq[Primer] = {
    val designParametersOption = Some(in.params)
    val outputTag = in.task.countTag
    val count: Int = output.get(outputTag).map(_.toInt).getOrElse {
      output.get("PRIMER_ERROR") match {
        case Some(primer3Error) => fail(s"Primer3 returned an error: $primer3Error", Seq.empty)
        case None               => fail(s"Primer3 did not return the output tag '$outputTag'", Seq.empty)
      }
    }

    (0 until count).map { num =>
      val (position, length) = PrimerDesignDef.splitInTwo(output(s"PRIMER_${side}_$num"), ',')
      val mapping = if (side == "LEFT") designRegion.resolve(position.toInt, length.toInt, strand=Strand.Positive)
      else designRegion.resolve(CoordMath.getStart(position.toInt, length.toInt), length.toInt, strand=Strand.Negative)

      // Re-make the primer sequences from the un-masked design sequence just in case
      var bases = unmaskedDesignSequence.substring(designRegion.project(mapping.start)-1, designRegion.project(mapping.end))
      if (mapping.strand == Strand.Negative) bases = Sequences.revcomp(bases)

      Primer(
        bases   = bases,
        tm      = output(s"PRIMER_${side}_${num}_TM").toDouble,
        penalty = output(s"PRIMER_${side}_${num}_PENALTY").toDouble,
        mapping = mapping,
        designParameters = designParametersOption
      )
    }
  }

  private def fail(message: String, errorLines: Seq[String] = Seq.empty): Nothing = {
    errorLines match {
      case Seq() => throw new IllegalArgumentException(message)
      case lines => throw new IllegalArgumentException(message + ":\n\t" + lines.map(l => "\t" + l).mkString("\n"))
    }
  }

  /** Regular expression used to parse the failure reason and count from Primer3 */
  private val failureRegex = "^ ?(.+) ([0-9]+)$".r

  /** Method to extract the reasons that the individual primers, and the pairs considered by
    * Primer3 failed (when there were failures).  The set of failures are returned sorted
    * from those with most failures to those with least.
    */
  private[primerdesign] def buildFailures(input: Primer3Parameters, dinucFailures: IterableOnce[Primer])
                                         (failureStrings: String*) : Seq[Primer3Failure] = {
    val byFailure = mutable.Map[Primer3FailureReason,Primer3Failure]()

    // Parse all the failure strings and merge together counts for the same kinds of failures
    failureStrings.flatMap(_.split(',')).map(_.trim).foreach {
      case failureRegex(explanation, count) =>
        if ("ok" != explanation && "considered" != explanation) {
          Primer3FailureReason.fromReason(explanation) match {
            case None         => logger.warning("Unknown Primer3 Failure Reason: ", explanation)
            case Some(reason) => byFailure.get(reason) match {
              case None    => byFailure(reason) = Primer3Failure(reason, count.toInt)
              case Some(f) => byFailure(reason) = f.copy(count = f.count + count.toInt)
            }
          }
        }
    }

    // Count how many individual primers failed for dinuc runs
    dinucFailures.iterator.toSet.count(p => Sequences.longestDinuc(p.bases).length > input.primerMaxDinucBases) match {
      case 0 => ()
      case n => byFailure(Primer3FailureReason.LongDinuc) = Primer3Failure(reason=Primer3FailureReason.LongDinuc, count=n)
    }

    // Extract the set of failures and sort them
    byFailure.values.toIndexedSeq.sortWith((lhs, rhs) => lhs.count >= rhs.count)
  }

  /** Closes all the resource related to the running Primer3 instance. */
  override def close(): Unit = {
    this.p3Input.close()
    this.p3Output.close()
    if (this.process.isAlive) this.process.destroy()
  }
}
