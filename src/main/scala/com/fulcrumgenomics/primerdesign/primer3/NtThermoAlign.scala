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

import java.io.{BufferedReader, InputStreamReader}
import java.nio.file.{Path, Paths}
import java.util.concurrent.TimeUnit

import com.fulcrumgenomics.FgBioDef.FilePath
import com.fulcrumgenomics.commons.util.LazyLogging

import scala.collection.mutable.ListBuffer
import scala.io.Source

object NtThermoAlign {
  val DefaultNtThalExecutable: FilePath = Paths.get("ntthal")
}

class NtThermoAlign(val executable: Path                 = NtThermoAlign.DefaultNtThalExecutable,
                    monovalentMilliMolar: Option[Double] = None,
                    divalentMilliMolar: Option[Double]   = None,
                    dntpMilliMolar: Option[Double]       = None,
                    dnaNanoMolar: Option[Double]         = None,
                    temperature: Option[Double]          = None,
                    timeout: Int                         = 5
                   ) extends Primer3Tool with LazyLogging {

  private val commandPrefix: Seq[String] = {
    val command = new ListBuffer[Any]()

    command.append(this.executable)
    command.append("-r")

    monovalentMilliMolar.foreach { m => command.appendAll(Seq("-mv", m)) }
    divalentMilliMolar.foreach   { d => command.appendAll(Seq("-dv", d)) }
    dntpMilliMolar.foreach       { n => command.appendAll(Seq("-n",  n)) }
    dnaNanoMolar.foreach         { d => command.appendAll(Seq("-d",  d)) }
    temperature.foreach          { t => command.appendAll(Seq("-t",  t)) }

    command.map(_.toString).toList
  }

  /** Calculates the duplex Tm of two sequences. */
  def duplexTm(s1: String, s2: String): Double = {
    val builder = new ProcessBuilder().redirectErrorStream(true)
    val command = commandPrefix ++ Seq("-s1", s1, "-s2", s2)
    val proc = builder.command(command:_*).start()
    proc.waitFor(timeout, TimeUnit.SECONDS) // Should run in << 1s, let alone 5s
    if (proc.isAlive) throw new IllegalStateException("ntthal took too long to run")
    if (proc.exitValue() != 0) {
      val error = Source.fromInputStream(proc.getInputStream).getLines().mkString("\n")
      throw new IllegalStateException(s"ntthal returned non-zero exit code ${proc.exitValue()}: \n$error")
    }
    else {
      val in = new BufferedReader(new InputStreamReader(proc.getInputStream), 16)
      val result = in.readLine().toDouble
      in.close()
      result
    }
  }
}
