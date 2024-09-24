# sigtk

A simple toolkit written for performing various operations on nanopore raw signal data.  *sigtk* is single threaded and is not optimised for performance. The intended use is to perform operations on relatively smaller datasets for learning purposes and eyeballing. Also, serves as examples for getting started with writing C programmes for nanopore signal analysis with BLOW5. The documentation is very brief at the moment, so just open an issue to get something clarified.

## Building

```
sudo apt-get install zlib1g-dev   #install zlib development libraries
git clone https://github.com/hasindu2008/sigtk
cd sigtk
make
```

The commands to install zlib __development libraries__ on some popular distributions:
```sh
On Debian/Ubuntu : sudo apt-get install zlib1g-dev
On Fedora/CentOS : sudo dnf/yum install zlib-devel
On OS X : brew install zlib
```


## Usage

## synthetic reference (sref)

Prints a synthetic reference signal for a given reference genome using traditional pore models. The 6-mer DNA pore-model used is [here](test/models/r9.4_450bps.nucleotide.6mer.template.model) and the 5-mer RNA pore-model is [here](test/models/r9.4_70bps.u_to_t_rna.5mer.template.model).

Usage:  `sigtk sref reference.fa`

Specify `--rna` to use the RNA pore-model. Output is a tab-delimited text file with each row being a reference contig (one for + and another for - when DNA; only for + when RNA) and the columns being as described below:

|Col|Type  |Name            |Description                                                        |
|--:|:----:|:------:        |:-----------------------------------------                         |
|1  |string|ref_name        |Reference contig name                                              |
|2  |int   |ref_len         |Length of the reference (no. of bases)                             |
|3  |char  |strand          |The reference strand direction (+ or -)                            |
|4  |int   |sig_len         |Length of the synthetic signal (no. of k-mers)                     |
|5  |float*|sig_mean        |Command separated mean current values of the synthetic signal      |


## Per-record raw-signal operations

The subtools in this section perform various operations on individual raw-signal records in a [BLOW5/SLOW5 file](https://www.nature.com/articles/s41587-021-01147-4). Those subtools can be used in one of the following forms:

```
# To perform subtool operation on all reads in a BLOW5 file
sigtk <subtool> reads.blow5
# To perform subtool operation on specified read IDs in a BLOW5 file
sigtk <subtool> reads.blow5 read_id1 read_id2 ..
```

By default, a tab-delimited text file with the first row being the header is printed. You can suppress the header using `-n` flag, for easy use with command line tools such as *awk*. Some subtools can be invoked with *-c* for compact output that prints data in a custom encoding (explained in each subtool, if relevant). These subtools automatically detect if raw signal data in for DNA or RNA, if applicable.

### pa

Prints the raw signal in pico-amperes.

|Col|Type  |Name            |Description                                                            |
|--:|:----:|:------:        |:-----------------------------------------                             |
|1  |string|read_id         |Read identifier name                                                   |
|2  |int   |len_raw_signal  |The number of samples in the raw signal                                |
|3  |float*|pa              |Comma separated Raw signal in pico amperes                             |

### event

Event segmentation is based on the method in Oxford Nanopore's [Scrappie basecaller](https://github.com/nanoporetech/scrappie).

By default, the output will be in the long intuitive form as explained below:

|Col|Type  |Name            |Description                                                        |
|--:|:----:|:------:        |:-----------------------------------------                         |
|1  |string|read_id         |Read identifier name                                               |
|2  |int   |event_idx       |Event index (0-based)                                              |
|3  |int   |raw_start       |Raw signal start index for the event (0-based; BED-like; closed)   |
|4  |int   |raw_end         |Raw signal end index for the event (0-based; BED-like; open)       |
|5  |float |event_mean      |Mean level of pico-ampere scaled signal for the event              |
|6  |float |event_std       |Standard deviations of pico-ampere scaled signal for the event     |

To obtain a condensed output that consumes less space and one record per row, specify `-c` option:

|Col|Type  |Name            |Description                                                            |
|--:|:----:|:------:        |:-----------------------------------------                             |
|1  |string|read_id         |Read identifier name                                                   |
|2  |int   |len_raw_signal  |The number of samples in the raw signal                                |
|3  |int   |raw_start       |Raw signal start index of the first event (0-based; BED-like; closed)  |
|4  |int   |raw_end         |Raw signal end index of the last event (0-based; BED-like; open)       |
|5  |int   |num_event       |Number of events                                                       |
|6  |int*  |events          |Comma separated event lengths (based on no. raw signal samples)       |

The event 0 starts at raw signal index `raw_start` (0-based; BED-like; closed) and ends at `raw_start+events[0]` (0-based; BED-like; open).
The event 1 starts at raw signal index `raw_start+events[0]` (0-based; BED-like; closed) and ends at `raw_start+events[0]+events[1]` (0-based; BED-like; open). Likewise, the events can be reconstructed by using the cumulative sum of `events`.

### stat

Prints signal statistics.

|Col|Type  |Name            |Description                                                            |
|--:|:----:|:------:        |:-----------------------------------------                             |
|1  |string|read_id         |Read identifier name                                                   |
|2  |int   |len_raw_signal  |The number of samples in the raw signal                                |
|3  |float |raw_mean        |Mean of raw signal values                                              |
|4  |float |pa_mean         |Mean of pico-amperes scaled signal                                     |
|5  |float |raw_std         |Standard deviation of raw signal values                                |
|6  |float |pa_std          |Standard deviation of pico-amperes scaled signal                       |
|7  |int   |raw_median      |Median of raw signal values                                            |
|8  |float |pa_median       |Mean of pico-amperes scaled signal                                     |

### prefix

Under development. Only for direct RNA at the moment.
Finds prefix segments in a raw signal such as adaptor and polyA.

|Col|Type  |Name            |Description                                                            |
|--:|:----:|:------:        |:-----------------------------------------                             |
|1  |string|read_id         |Read identifier name                                                   |
|2  |int   |len_raw_signal  |The number of samples in the raw signal                                |
|3  |int   |adapt_start     |Raw signal start index of the adaptor                                  |
|4  |int   |adapt_end       |Raw signal end index of the adaptor                                    |
|5  |int   |polya_start     |Raw signal start index of the polyA tail                               |
|6  |int   |polya_end       |Raw signal end index of the polyA tail                                 |

If `--print-stat` is printed, following additional columns will be printed.

|Type  |Name            |Description                                                            |
|:----:|:------:        |:-----------------------------------------                             |
|float |adapt_mean      |Mean of pico-amperes scaled signal of the adaptor                      |
|float |adapt_std       |Standard deviation of pico-amperes scaled signal of the adaptor        |
|float |adapt_median    |Median of pico-amperes scaled signal of the adaptor                    |
|float |polya_mean      |Mean of pico-amperes scaled signal of the polyA tail                   |
|float |polya_std       |Standard deviation of pico-amperes scaled signal of the polyA tail     |
|float |polya_median    |Median of pico-amperes scaled signal of the polyA                      |

### jnn

Under development.  Print segments found using JNN segmenter.

|Col|Type  |Name            |Description                                                            |
|--:|:----:|:------:        |:-----------------------------------------                             |
|1  |string|read_id         |Read identifier name                                                   |
|2  |int   |len_raw_signal  |The number of samples in the raw signal                                |
|3  |int   |num_seg         |Number of segments found                                               |
|4  |string|seg             |List of segments as explained below                                    |

```
...............|..........|..............|..........|............   <- signal and segments
              100        110            201        212              <- signal index (0-based)
```
Segments will be noted as:
`100,110;201,212;`

If `-c` is specified, output will be in the following short notation by using relative offsets.
```
...............|..........|..............|..........|............   <- signal and segments
              100        110            201        212              <- signal index (0-based)

               <---10----><-----91------><---11----->
```

`100H10,91H11,`


### ss

Under development. Operations to convert to/from signal alignment string (ss). See https://hasindu2008.github.io/f5c/docs/output#resquiggle-paf-output-format for explanation of ss.

To convert a PAF file with ss tags to TSV, you can use:
```
sigtk ss paf2tsv in.paf
```


### ent

Calculates shannon entropy for reads in a given S/BLOW5 file.

|Col|Type  |Name            |Description                                                            |
|--:|:----:|:------:        |:-----------------------------------------                             |
|1  |string|read_id         |Read identifier name                                                   |
|2  |float   |raw_ent  |entropy of raw signal samples                                |
|3  |float   |delta_ent     |entropy after zig-zag delta                                  |
|4  |float   |byte_ent       |entropy after splitting and storing least significant byte and most significant byte of the zig-zag delta values separately: ent(LSB)+ent(MSB)                                    |

### qts

Quantise the raw signal in a S/BLOW5 files. Takes a S/BLOW5 file as the input and writes the quantised output to a S/BLOW5 file.

Usage:
```
sigtk qts original.blow5 -o quantised.blow5
```

Options:

- `-b INT` : Number of lower significant bits eliminate. Default is 1.
- `-m [floor|round|fill-ones]`: quantisation method. Default is round.


## Acknowledgement

The event detection code is from Oxford Nanopore's [Scrappie basecaller](https://github.com/nanoporetech/scrappie).
The pore-models are from [Nanopolish](https://github.com/jts/nanopolish).
Code snippets have been taken from [Minimap2](https://github.com/lh3/minimap2) and [Samtools](http://samtools.sourceforge.net/). The name of the tool *sigtk* in signal-space was inspired by [seqtk](https://github.com/lh3/seqtk) in base-space. Kseq and ksort from [klib](https://github.com/attractivechaos/klib) are used. Segmentation method (aka jnn) was adapted from [SquiggleKit](https://github.com/Psy-Fer/SquiggleKit) and [deeplexicon](https://github.com/Psy-Fer/deeplexicon).

## Citation

sigtk highly relies on slow5lib and SLOW5 format and there is no direct publication. However, you may cite the *slow5lib/pyslow5* publication if you found sigtk useful.

> Gamaarachchi, H., Samarakoon, H., Jenner, S.P. et al. Fast nanopore sequencing data analysis with SLOW5. Nat Biotechnol 40, 1026-1029 (2022). https://doi.org/10.1038/s41587-021-01147-4

```
@article{gamaarachchi2022fast,
  title={Fast nanopore sequencing data analysis with SLOW5},
  author={Gamaarachchi, Hasindu and Samarakoon, Hiruna and Jenner, Sasha P and Ferguson, James M and Amos, Timothy G and Hammond, Jillian M and Saadat, Hassaan and Smith, Martin A and Parameswaran, Sri and Deveson, Ira W},
  journal={Nature biotechnology},
  pages={1--4},
  year={2022},
  publisher={Nature Publishing Group}
}
