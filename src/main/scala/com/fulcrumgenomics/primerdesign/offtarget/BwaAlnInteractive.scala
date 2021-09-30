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

package com.fulcrumgenomics.primerdesign.offtarget

import com.fulcrumgenomics.FgBioDef.{FilePath, PathToFasta}
import com.fulcrumgenomics.commons.util.LazyLogging
import com.fulcrumgenomics.util.{Io, Sequences}
import htsjdk.samtools._
import htsjdk.samtools.util.{LineReader, StringUtil}

import java.io._
import java.nio.file.Paths
import java.util.Collections
import scala.collection.mutable.ListBuffer

/**
  * Container for case classes that are used as inputs and outputs to Bwa.
  */
object BwaAlnInteractive {
  /** Default value that will be used as the executable path for bwa. */
  val DefaultBwaAlnInteractiveExecutable: FilePath = Paths.get("bwa")

  /** Represents a single query sequence for mapping. */
  case class Query(id: String, bases: String) {
    def toFastq(reverseComplement: Boolean = false) : String = {
      val baseString: String = if (reverseComplement) Sequences.revcomp(bases) else bases
      "@" + id + "\n" + baseString + "\n" + "+" + "\n" + StringUtil.repeatCharNTimes('H', bases.length) + "\n"
    }
  }

  /**
    * Encapsulates all the results for a single query
    *
    * @param query the query (as given, no RC'ing)
    * @param hitCount the total number of hits found by bwa (may be more than hits.length)
    * @param hits the subset of hits that were reported by bwa
    */
  case class Result(query: Query, hitCount: Int, hits: Seq[Hit])

  /** Encapsulates a single hit of a query sequence to the genome. */
  object Hit {
    /**
      * Generates a hit object.
      *
      * @param chrom the chromosome of the hit
      * @param start the start position of the hit
      * @param negative whether the hit is on the negative strand
      * @param cigar the cigar string returned by BWA
      * @param edits the number of edits between the read and the reference
      * @param rc whether the reverse-complement of the query sequence was fed to BWA, in which case various
      *           fields should be negated/reversed in the Hit
      * @return A Hit that represents the mapping of the original query sequence that was supplied
      */
    def apply(chrom: String, start: Int, negative: Boolean, cigar: String, edits: Int, rc:Boolean=false): Hit = {
      val (strand, cig) = if (rc) {
        val elems = new java.util.ArrayList(TextCigarCodec.decode(cigar).getCigarElements)
        Collections.reverse(elems)
        (!negative, new Cigar(elems))
      }
      else {
        (negative, TextCigarCodec.decode(cigar))
      }

      new Hit(chrom=chrom, start=start, negative=strand, cigar=cig, edits=edits)
    }
  }

  /** Represents a single hit or alignment of a sequence to a location in the genome. */
  case class Hit(chrom: String, start: Int, negative: Boolean, cigar: Cigar, edits: Int) {
    lazy val mismatches: Int = {
      import com.fulcrumgenomics.commons.CommonsDef.javaIteratorAsScalaIterator
      edits - cigar.getCigarElements.iterator.count(e => e.getOperator.isIndel)
    }
    lazy val end: Int = start + cigar.getReferenceLength - 1
  }
}

/**
  * Wrapper around a novel mode of 'bwa aln' that allows for "interactive" use of bwa to keep
  * the process running and be able to send it chunks of reads periodically and get alignments
  * back without waiting for a full batch of reads to be sent.
  *
  * See: https://github.com/fulcrumgenomics/bwa/tree/interactive_aln
  *
  * @param bwa path to the bwa executable
  * @param ref path to the index reference fasta
  * @param maxMismatches the maximum number of mismatches allowed in the full query sequence
  * @param maxMismatchesInSeed the maximum number of mismatches allowed in the seed region
  * @param maxGapOpens the max number of gap opens allowed in the full query sequence
  * @param maxGapExtensions the max number of gap extensions allowed in the full query sequence
  * @param seedLength the length of the seed region
  * @param maxHits the maximum number of hits to report - if more than this number of seed hits
  *                are found, report only the count and not each hit
  * @param reverseComplement reverse complement each query sequence before alignment
  * @param includeAltHits if true include hits to chromosomes with names ending in _alt, otherwise don't
  */
class BwaAlnInteractive(val bwa: FilePath = BwaAlnInteractive.DefaultBwaAlnInteractiveExecutable,
                        val ref: PathToFasta,
                        val maxMismatches: Int = 3,
                        val maxMismatchesInSeed: Int = 3,
                        val maxGapOpens: Int = 0,
                        val maxGapExtensions: Int = -1,
                        val seedLength: Int = 20,
                        val maxHits: Int,
                        val reverseComplement: Boolean = false,
                        val includeAltHits: Boolean = false,
                        val threads: Int = 8
         ) extends LazyLogging with AutoCloseable {
  import BwaAlnInteractive._

  // -N = non-iterative mode: search for all n-difference hits (slooow)
  // -S = output SAM (run samse)
  // -Z = interactive mode (no input buffer and empty lines between records force processing)
  private val args = Seq(bwa, "aln" , "-t", threads, "-n", maxMismatches, "-o", maxGapOpens, "-e", maxGapExtensions,
    "-l", seedLength, "-k", maxMismatchesInSeed, "-N", "-S", "-X", maxHits, "-Z", "-D", ref, Io.StdIn).map(_.toString)

  private val process = new ProcessBuilder(args:_*).start()
  private val errorPipe = Io.pipeStream(process.getErrorStream, s => logger.debug("bwa: ", s))
  private val toBwa     = new BufferedWriter(new OutputStreamWriter(process.getOutputStream))
  private val reader    = new BufferedReader(new InputStreamReader(process.getInputStream))
  private val header    = readHeader(this.reader)
  private val parser    = new SAMLineParser(new DefaultSAMRecordFactory, ValidationStringency.SILENT, header, null, null)

  /** Maps a single query to the genome and returns the result. */
  def map(query: String, id: String = "unknown"): Result = map(Seq(Query(id=id, bases=query))).head

  /**
    * Maps zero or more queries to the genome in a batch. Results are identical to calling map() on each
    * query, but execution is more efficient.
    *
    * @param queries zero or more query sequences
    * @return a seq the same length and order as the input queries with one Result per Query
    */
  def map(queries: Seq[Query]): Seq[Result] = {
    if (queries.isEmpty) {
      Seq.empty
    }
    else {
      // Send the reads to BWA
      queries.foreach(q => toBwa.write(q.toFastq(this.reverseComplement)))
      signalBwa()

      // Read back the results
      val buffer = ListBuffer[Result]()
      queries.foreach(q => {
        val rec = this.parser.parseLine(this.reader.readLine())
        assert(q.id == rec.getReadName, s"Query and Results are out of order. Query=${q.id}, Result=${rec.getReadName}")

        val result = (rec.getReadUnmappedFlag, Option(rec.getIntegerAttribute("HN"))) match {
          case (true, _)                                        => Result(query=q, hitCount=0, hits=Nil)
          case (false, Some(hitCount)) if hitCount > maxHits => Result(query=q, hitCount=hitCount, hits=Nil)
          case (false, Some(hitCount)) =>
            val firstHit = Hit(rec.getContig, rec.getAlignmentStart, rec.getReadNegativeStrandFlag,
              rec.getCigarString, rec.getIntegerAttribute("NM"), rc=reverseComplement)

            // E.g. XA:Z:chr4,-97592047,24M,3;chr8,-32368749,24M,3;
            val extraHits = Option(rec.getStringAttribute("XA"))
              .toList
              .flatMap(_.split(';'))
              .filter(s => s.nonEmpty)
              .map(s => s.split(','))
              .map(fs => Hit(chrom=fs(0), start=fs(1).substring(1).toInt, negative=fs(1).charAt(0) == '-',
                cigar=fs(2), edits=fs(3).toInt, rc=reverseComplement))

            val hits = firstHit :: extraHits
            val filteredHits = hits.filter(h => this.includeAltHits || !h.chrom.endsWith("_alt"))

            Result(query=q, hitCount=if (filteredHits.isEmpty) hitCount else filteredHits.size, hits=filteredHits)
          case (unmapped, hitCount) =>
            throw new IllegalStateException(s"Shouldn't be able to get unmapped=$unmapped and hitCount=$hitCount")
        }

        buffer += result
      })

      buffer.toList
    }
  }

  private def signalBwa() : Unit = {
    (1 to 3).foreach(_ => {
      toBwa.flush()
      toBwa.newLine()
      toBwa.newLine()
      toBwa.flush()
    })
  }


  /** Reads the header from a stream of lines, terminating when it sees the @PG line, so as not to
    * hang waiting for the first alignments to happen.
    */
  private def readHeader(reader: BufferedReader): SAMFileHeader = {
    val lines = readHeaderLines(reader)

    // Stupid line reader to make SAMTextHeaderCodec happy
    val lineReader = new LineReader() {
      private val iterator       = lines.zipWithIndex.iterator.buffered
      private var lastLineNumber = 0
      override def peek(): Int = lines.headOption.map(l => l.charAt(0).toInt).getOrElse(-1)
      override def getLineNumber: Int = lastLineNumber
      override def close(): Unit = ()
      override def readLine(): String = {
        if (iterator.hasNext) {
          val (line, index) = iterator.next()
          lastLineNumber = index + 1
          line
        }
        else {
          null
        }
      }
    }

    new SAMTextHeaderCodec().decode(lineReader, "n/a")
  }

  /** Reads the header lines, stopping at the @PG line which should be the last. */
  private def readHeaderLines(reader: BufferedReader): Seq[String] = {
    val lines = Seq.newBuilder[String]
    var keepReading = true

    while (keepReading) {
      keepReading = reader.readLine() match {
        case null => false
        case s    => lines += s; !s.startsWith("@PG")
      }
    }

    lines.result()
  }

  /** Destroys the underlying process and related IO objects. */
  override def close(): Unit = {
    toBwa.close()
    this.reader.close()
    process.destroy()
    errorPipe.close(2000)
  }
}

