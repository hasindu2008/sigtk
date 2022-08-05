# sigtk

A simple toolkit written for performing various operations on nanopore raw signal data. This is still in a very premature development state. The command line interface is therefore subject to frequent change. Currently, *sigtk* it is single threaded and has not been optimised for performance. The intended use is to perform operations on relatively smaller datasets for learning purposes and eyeballing.

## Building


```
sudo apt-get install zlib1g-dev   #install zlib development libraries
git clone https://github.com/hasindu2008/sigtk
cd sigtk
make
```

The commands to zlib __development libraries__ on some popular distributions :
```sh
On Debian/Ubuntu : sudo apt-get install zlib1g-dev
On Fedora/CentOS : sudo dnf/yum install zlib-devel
On OS X : brew install zlib
```


## Usage

## synthetic reference (sref)

Prints a synthetic reference signal for a given reference genome using traditional pore models. The 6-mer DNA pore-model used is [here](test/models/r9.4_450bps.nucleotide.6mer.template.model) and the 5-mer RNA pore-model is [here](r9.4_70bps.u_to_t_rna.5mer.template.model).

Usage:  `sigtk sref reference.fa`

Specify `--rna` to use the RNA pore-model. Output is a tsv file with each row being a reference contig (one for + and another for -) and the columns being as described below:

|Col|Type  |Name            |Description                                                        |
|--:|:----:|:------:        |:-----------------------------------------                         |
|1  |string|ref_name        |Reference contig name                                              |
|2  |int   |ref_len         |Length of the reference (no. of bases)                             |
|3  |char  |strand          |The reference strand direction (+ or -)                            |
|4  |int   |sig_len         |Length of the synthetic signal (no. of k-mers)                     |
|5  |float*|sig_mean        |Command separated mean current values of the synthetic signal      |


## Per-record raw-signal operations

The subtools in this section perform various operations on individual raw-signal records in a BLOW5/SLOW5 file. They take one of the following forms:

```
# To perform subtool operation on all reads in a BLOW5 file
sigtk <subtool> reads.blow5
# To perform subtool operation on specified read IDs in a BLOW5 file
sigtk <subtool> reads.blow5 read_id1 read_id2 ..
```

By default, a tab-delimited text file with the first row being the header is printed. You can suppress the header using `-n` flag for easy use with command line tools such as *awk*. Some subtools can be invoked with *-c* for compact output that prints data in a custom encoding explained in each subtool if relevant. These subtools automatically detect if raw signal data in for DNA or RNA if applicable.

## pa

Prints the raw signal in pico-amperes.

|Col|Type  |Name            |Description                                                            |
|--:|:----:|:------:        |:-----------------------------------------                             |
|1  |string|read_id         |Read identifier name                                                   |
|2  |int   |len_raw_signal  |The number of samples in the raw signal                                |
|3  |float*|pa              |Comma separated Raw signal in pico amperes                             |

## event

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
|6  |int*  |events          |Comma separated event lengths (based on no. raw signals samples)       |

The event 0 starts at raw signal index `raw_start` (0-based; BED-like; closed) and ends at `raw_start+event[0]` (0-based; BED-like; open).
The event 1 starts at raw signal index `raw_start+event[0]` (0-based; BED-like; closed) and ends at `raw_start+event[0]+event[1]` (0-based; BED-like; open). Likewise, the events can be reconstructed by using the cumulative sum of `events`.

## stat

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

## prefix

Under construction. Will change anytime.
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

## jnn

Under construction. Will change anytime.

Print segments found using JNN segmenter.

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
Segments will be noted as below:
100,110;201,212;

If `-c` is specified, output will be in the following short notation by using relative offsets.
```
...............|..........|..............|..........|............   <- signal and segments
              100        110            201        212              <- signal index (0-based)

               <---10----><-----91------><---11----->
```

100H,10,91H,11,


## Acknowledgement

The event detection code is from Oxford Nanopore's [Scrappie basecaller](https://github.com/nanoporetech/scrappie).
The pore-models are from [Nanopolish](https://github.com/jts/nanopolish).
Some code snippets have been taken from [Minimap2](https://github.com/lh3/minimap2), [Samtools](http://samtools.sourceforge.net/). The name of the tool *sigtk* in signal-space was inspired by [seqtk](https://github.com/lh3/seqtk) in base-space. Kseq and ksort from [klib](https://github.com/attractivechaos/klib) are used. Segmentation method (aka jnn) was adapted from [SquiggleKit](https://github.com/Psy-Fer/SquiggleKit) and [deeplexicon](https://github.com/Psy-Fer/deeplexicon).
