# NMDetective_B
Reimplementation of the NMDetective_B algorithm from Lindeboom et al., 2019 Nat Gen.

## Install

No setup.py (to come)

In the meantime, put NMDetective_B/bin on your PATH

## Usage
```
nmdetect -h
usage: nmdetect [-h] -g GTF -l LABEL

options:
  -h, --help            show this help message and exit

required arguments:
  -g GTF, --gtf GTF     GTF file to parse (can be in gzip format, detected if .gz is the end of the filename
  -l LABEL, --label LABEL
                        a label for the output files

Example usage: nmdetect -g <gtf_file> -l <outputlabel>
```
## Example

Assuming you have the human gencode GTF:
```
nmdetect -g gencode.v42.annotation.gtf.gz -l gencode_v42
INFO    : Arguments:
INFO    :   gtf: gencode.v42.annotation.gtf.gz
INFO    :   label: gencode_v42
INFO    : Starting score
INFO    : 1,000,000 done
INFO    : 2,000,000 done
INFO    : 3,000,000 done
INFO    : Processed 3,007,112 entries in the GTF
INFO    : Found 170,279 transcripts
INFO    : Skipped as no identifiable START/STOP: 90,263 transcripts
INFO    : Drew Pie chart
```

## Benchmarking

Expected NMD scores from the Lindeboom UCSC track and GTF:

![Expected](https://github.com/oaxiom/NMDetective_B/blob/9dc73694409bbe4806aafdec0df850c05716cb9c/observed_expected/expected.Lindeboom.png)

Observed NMD scores from nmdetect

![Observed](https://github.com/oaxiom/NMDetective_B/blob/9dc73694409bbe4806aafdec0df850c05716cb9c/observed_expected/observed.pie.png)

## Reference

GENCODE hg38 (v42)

![GENCODE hg38](https://github.com/oaxiom/NMDetective_B/blob/97046bed4f7c219b9209c8186b631f1a30c27b7f/images/gencode_v42.pie.pdf)

GENCODE mm10 (vM20)

![GENCODE mm10](https://github.com/oaxiom/NMDetective_B/blob/97046bed4f7c219b9209c8186b631f1a30c27b7f/images/gencode_vM20.pie.pdf)
