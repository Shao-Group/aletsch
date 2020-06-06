# Introduction

meta-scallop implements an efficient algorithm to assemble multiple RNA-seq samples.
It uses splice graph and phasing paths as underlying data strctures to represent
the alignments of each gene loci in individual RNA-seq samples.
Efficient algorithms are implemented to combine splice graphs (and phasing paths)
at overlapped gene loci. Eventually, the core algorithm used in Scallop (i.e., phase-preserving decomposition)
is employed to decompose the combined splice graphs to transcripts.

# Release
Latest release of meta-scallop is [v0.1.5](https://github.com/Shao-Group/meta-scallop/releases/tag/v0.1.5).

# Installation
Download the source code of meta-scallop from
[here](https://github.com/Shao-Group/meta-scallop/releases/download/v0.1.5/meta-scallop-0.1.5.tar.gz).
meta-scallop uses additional libraries of Boost and htslib. 
If they have not been installed in your system, you first
need to download and install them. You might also need to
export the runtime library path to certain environmental
variable (for example, `LD_LIBRARY_PATH`, for most linux distributions).
After install these dependencies, you then compile the source code of meta-scallop.
If some of the above dependencies are not installed to the default system 
directories (for example, `/usr/local`, for most linux distributions),
their corresponding installing paths should be specified to `configure` of meta-scallop.

## Download Boost
If Boost has not been downloaded/installed, download Boost
[(license)](http://www.boost.org/LICENSE_1_0.txt) from (http://www.boost.org).
Uncompress it somewhere (compiling and installing are not necessary).

## Install htslib
If htslib has not been installed, download htslib 
[(license)](https://github.com/samtools/htslib/blob/develop/LICENSE)
from (http://www.htslib.org/) with version 1.5 or higher.
(Note that htslib relies on zlib. So if zlib has not been installed in your system,
you need to install zlib first.) 

Use the following commands to build htslib:
```
./configure --disable-bz2 --disable-lzma --disable-gcs --disable-s3 --enable-libcurl=no
make
make install
```
The default installation location of htslib is `/usr/lib`.
If you would install it to a different location, replace the above `configure` line with
the following (by adding `--prefix=/path/to/your/htslib` to the end):
```
./configure --disable-bz2 --disable-lzma --disable-gcs --disable-s3 --enable-libcurl=no --prefix=/path/to/your/htslib
```
In this case, you also need to export the runtime library path (note that there
is an additional `lib` following the installation path):
```
export LD_LIBRARY_PATH=/path/to/your/htslib/lib:$LD_LIBRARY_PATH
```

## Compile meta-scallop

Use the following to compile meta-scallop:
```
./configure --with-htslib=/path/to/your/htslib --with-boost=/path/to/your/boost
make
```

If some of the dependencies are installed in the default system directory (for example, `/usr/lib`),
then the corresponding `--with-` option might not be necessary.
The executable file `meta-scallop` will appear at `meta/meta-scallop`.


# Usage

The usage of `meta-scallop` is:
```
./meta-scallop -i <input-bam-list> -o <output.gtf> [options]
```

The `input-bam-list` specifies the list of aligned samples (in .bam format).
Note meta-scallop will try to infer the `library_type` of each individual sample
using the `XS` tag stored in the input bam files. 
Also make sure that they are sorted; otherwise run `samtools` to sort it (`samtools sort input.bam > input.sort.bam`).
The assembled transcripts from all these samples will be written to `output.gtf`.
If `-d directory` option is provided, the assembled transcripts for each individual
sample will be generated under the specified directory. 
If `-D directory` option is provided, the bridged alignment files for each individual
sample will be generated under the specified directory. NOTE: make sure the specified
directory exists, both for `-d` and `-D`, as `meta-scallop` will NOT create them in the program.

meta-scallop support the following parameters. It also supports all parameters
that needed in the core algorithm of Scallop. Dry run `meta-scallop` for more
details.

 Parameters | Default Value | Description
 ------------------------- | ------------- | ----------
 --help  | | print usage of meta-scallop and exit
 --version | | print version of meta-scallop and exit
 -t | 10  | number of threads
 -m |     | if specified, invididual samples will be processed by two threads
 -b | 100 | the number of samples that will be process at a time
 -c | 100 | the maximized number of splice graphs that will be combined
 -s | 0.3 | the minimized similarity between two splice graphs that can be combined
 -d |     | the directory in which assembled transcripts for individual sample will be generated
 -D |     | the directory in which bridged alignments for individual sample will be generated
