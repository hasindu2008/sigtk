# Development notes

## prefix

Under construction. Will change anytime. Only for direct RNA at the moment.
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
Segments will be noted as:
`100,110;201,212;`

If `-c` is specified, output will be in the following short notation by using relative offsets.
```
...............|..........|..............|..........|............   <- signal and segments
              100        110            201        212              <- signal index (0-based)

               <---10----><-----91------><---11----->
```

`100H10,91H11,`

