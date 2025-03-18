[![codecov](https://codecov.io/gh/fulcrumgenomics/fgprimer/branch/master/graph/badge.svg)](https://codecov.io/gh/fulcrumgenomics/fgprimer)
[![Maven Central](https://maven-badges.herokuapp.com/maven-central/com.fulcrumgenomics/fgprimer_2.13/badge.svg)](https://maven-badges.herokuapp.com/maven-central/com.fulcrumgenomics/fgprimer_2.13)
[![Javadocs](http://javadoc.io/badge/com.fulcrumgenomics/fgprimer_2.13.svg)](http://javadoc.io/doc/com.fulcrumgenomics/fgprimer_2.13)
[![License](http://img.shields.io/badge/license-MIT-blue.svg)](https://github.com/fulcrumgenomics/fgprimer/blob/master/LICENSE)
[![Language](http://img.shields.io/badge/language-scala-brightgreen.svg)](http://www.scala-lang.org/)

fgprimer
========

A Scala Library for primer design and related activities for genomic assays.

<p>
<a href float="left"="https://fulcrumgenomics.com"><img src=".github/logos/fulcrumgenomics.svg" alt="Fulcrum Genomics" height="100"/></a>
</p>

[Visit us at Fulcrum Genomics](https://www.fulcrumgenomics.com) to learn more about how we can power your Bioinformatics with fgprimer and beyond.

<a href="mailto:contact@fulcrumgenomics.com?subject=[GitHub inquiry]"><img src="https://img.shields.io/badge/Email_us-brightgreen.svg?&style=for-the-badge&logo=gmail&logoColor=white"/></a>
<a href="https://www.fulcrumgenomics.com"><img src="https://img.shields.io/badge/Visit_Us-blue.svg?&style=for-the-badge&logo=wordpress&logoColor=white"/></a>

This readme document is mostly for developers/contributors and those attempting to build the project from source.

Detailed developer documentation can be found [here](http://javadoc.io/doc/com.fulcrumgenomics/fgprimer_2.13).

<!---toc start-->
  * [Overview](#overview)
  * [Building](#building)
  * [Command line](#command-line)
  * [Include fgprimer in your project](#include-fgprimer-in-your-project)
  * [Contributing](#contributing)
  * [Authors](#authors)
  * [License](#license)

<!---toc end-->

## Overview

Fgprimer is a Scala library and set of command line tools for primer design for genomic assay design.
The collection of tools and associated API within `fgprimer` are used by our customers and others both for ad-hoc data analysis and within production pipelines.

## Building 
### Cloning the Repository

To clone the repository: `git clone https://github.com/fulcrumgenomics/fgprimer.git`

### Dependencies
[`primer3`](https://github.com/primer3-org/primer3) (version 2.5.0 or greater) is a dependency. Note that starting with
version 2.5 the primer3 config directory is no longer needed as the thermodynamic parameters are packaged into the
binaries!
Install from source or with [`conda`](https://conda.io/): `conda install -c bioconda primer3`

A custom [`bwa`](https://github.com/fulcrumgenomics/bwa/tree/interactive_aln) is a dependency.
Install it from source or with: 

```
git clone -b interactive_aln git@github.com:fulcrumgenomics/bwa.git
cd bwa
make -j 12
```

### Running the build
fgprimer is built using [sbt](http://www.scala-sbt.org/).

Use ```sbt assembly``` to build an executable jar in ```target/scala-2.13/```.

Tests may be run with ```sbt test```. 

Java SE 8 is required.


## Include fgprimer in your project

You can include `fgprimer` in your project using the latest release (not currently available):

```
"com.fulcrumgenomics" %% "fgprimer" % "0.0.1"
```

for the latest released version or (buyer beware):

```
"com.fulcrumgenomics" %% "fgprimer" % "0.0.1-<githash>-SNAPSHOT"
```

for the latest development snapshot.

## Contributing

Contributions are welcome and encouraged.
We will do our best to provide an initial response to any pull request or issue within one-week.
For urgent matters, please contact us directly.

## Authors

* [Tim Fennell](https://github.com/tfenne) (maintainer)
* [Nils Homer](https://github.com/nh13) (maintainer)

## License

`fgprimer` is open source software released under the [MIT License](https://github.com/fulcrumgenomics/fgprimer/blob/master/LICENSE).

