# RETIVAD
A RealTime Variant Detector software.

## Description
Since December 2019, the emergence of a novel coronavirus responsible for a severe acute respiratory syndrome (SARS-CoV-2) is accountable for a major pandemic situation, leading to date to nearly 3 million deaths world-wide. Early on 2020 major viral diagnostic efforts were deployed, including the development of news strategies aiming at (i) reducing the diagnostics time, (ii) decrease the costs per assay, (iii) providing population-scale solutions, but also (iv) presenting high sensitivity and specificity, notably to discriminate asymptomatic cases.
Beyond the diagnostics requirements, the identification of the highly transmissible SARS-CoV-2 variant B.1.1.7 - originated on the south of England - has strongly accelerated the world-wide interest on tracking SARS-CoV2 variants’ occurrence. More recently, two other extremely infectious variants, namely the P.1 issued on the Brazilian city of Manaus and the B.1.351 initially detected In South Africa, were described and unsurprisingly further others are expected to be discovered, notably due to the long period of time on which the pandemic situation is lasting.
We have recently developed a proof-of-concept study for population-scale diagnostics and variants tracing assay using massive parallel DNA sequencing performed with the [MinION Oxford Nanopore technology (ONT)](https://nanoporetech.com/) ([BioRXV](https://nanoporetech.com/)). ONT MinION sequencers are well-known for its portability, its reduced cost, and its simplicity of use; thus, representing a suitable strategy for its implementation even in remote places over the world. 
Herein we provide RETIVAD, the companion computational solution for the aforementioned assay, allowing to generate the diagnostics outcome on a real-time mode during sequencing; essential to earn time for diagnostics report on the current context of the pandemic situation.

## Conda version
### Dependencies

```
conda create --name retivade paramiko medaka minimap2 regex bwa samtools bedtools biopython bedops matplotlib
conda activate retivade
```

### How to run

```
require arguments:
  -5 cDNAbc, --cDNAbc cDNAbc
                           File containing the “cDNA” barcodes’ list. Format: seq-ID / Sequence
  -3 PCRbc, --PCRbc PCRbc
                           File containing the “PCR” barcodes’ list. Format: seq-ID / Sequence
  -g genome, --genome genome
                           Aligner genome index path
  ```
```
optional arguments:
  -h, --help               Show this help message and exit
  -f file, --file file     File(s) you want annalyse separate with commat (default: )
  --folder folder          Folder who contain fastq to annalyse (default: None)
  --tmp tmp                Folder who contain temporary files (default: /tmp)
  -c core, --core core     Number of cores allocated to the software (default: 1)
  -o outdir, --out outdir  Folder in which the output is saved (default: .)
  -n name, --name name     Name of the assay. By default, the software uses the date. (default: None)
  -G guibsonE, --guibsonE errors
                           Guibson error accepted (default: 4)
  -B barcodeE, --barcodeE errors
                           Barcode error accepted (default: 4)
  -i insertMinLength, --insertMinLength insertMinLength
                           Good insert minimum length (default: 30)
  -I insertMaxLength, --insertMaxLength insertMaxLength
                           Good insert maximum length (default: 10000)
  -l seqMinLength, --seqMinLength seqMinLength
                           Good sequences minimum length between guibson sequences (default: 70)
  -L seqMaxLength, --seqMaxLength seqMaxLength
                           Good sequences maximum length between guibson sequences (default: 10000)
  --noreverse3             If you send directly the reverse complement in 3' barcodes, activate this option (default: False)
  -a aligner, --aligner aligner
                           Aligner name (default: bwa)
  --alignerPath alignerPath
                           Aligner path. By default we use the name (default:None)
  --alignerSpecialsParam alignerSpecialsParam
                           Specials parameters for the aligner (default: None)
  --traceseq traceseq      Sequences you want track to check the method (default: None)
  --realtime               Put the software in realtime annalyse (default: False)
  --ext ext                Path to the configuration for the ssh. Use this parameter if your files are not in the computer (e.g. when connecting to the Mk1c sequencer).
                           This method are only in realtime and regeneraterealtime option are activated. (default: None)
  --timecheck timecheck    In realtime analysis, you can setup the time(m) interval between each check in the folder (default: 10)
  --timeEndCheck timeEndCheck
                           In realtime analysis, set up the time(m) between the last data collection and the end of the computation due to the absence of new files to process (default: 40)
  --dev                    Keep all temporary datas (default: False)
  --regeneraterealtime     Mimic realtime analysis from all fastq files based on the information concerning the date and time of writing (default: False)
  --mincov mincov          Minimum coverage to support a variation (default: 20)
  -v visugenome, --visugenome visugenome
                           Path to genome used for visualization in realtime. This option is relevant for cases in which the visualization is focused on a given region of interest (instead of the whole genome) (default: None)
  ```

#### Example
In this example we have a sequencing run taking place on a external device (eg: mk1-c), and we want analyze the data in real-time.
Take files from the template folder and the [sampleFiles](https://drive.google.com/file/d/1bDKZZvL6tbHQjnaILUaFdjGU_2jCxktY/view?usp=sharing):
- Barcodes datas files: barcode5.tsv (cDNA barcodes) and barcode3.tsv (PCR Barcodes).
- The pathogen genome of your interest. Here we select the Covid19 genome.
- The configSSH file is used when the fastq files are not in your computer. (Work only in real-time analysis)
- [sampleFiles](https://drive.google.com/file/d/1bDKZZvL6tbHQjnaILUaFdjGU_2jCxktY/view?usp=sharing) sequences extract from candidates. (Herein, samples correspond to synthetic Covid-19 RNA material, not from real candidates!) Decompress files in a folder named retivad_fastq_test.

###### Send the folder retivad_fastq_test in the external computer with a rsync or decompress the archive in the remote computer (to keep the modification date)
###### Configure configSSH.ini
###### Create genomes index with bwa: bwa index covid19.fa / bwa index covid19InterestRegion.fa
###### Use this command to launch the software and generate a realtime experiment:
```
python3 retivad.py -c 8 -5 template/barcode5.tsv -3 template/barcode3.tsv -g template/covid19.fa -v template/covid19InterestRegion.fa --folder path/to/retivad_fastq_test -n covidTest --ext configSSH.ini --regeneraterealtime > covid.log 2>&1 &
```
##### Wait and see
In real-time the software produces a plot displaying the number of aligned reads (corresponding to the defined target sequence) relative to the time (in hours); as well as further other stats and a summary.
If RETIVAD does not find new files after 40min, the software ends de loop of data collection. Finally, RETIVAD generates a global summary.

## Authors

|         |                                                                                               |
| ------- | --------------------------------------------------------------------------------------------- |
| Author  | STUDER Francois ([Github](https://github.com/studyfranco))                                    |
| Author  | ENGELEN Stéfan ([Github](https://github.com/sengelen))                                        |
| Author  | MENDOZA PARRA Marco ([Github](https://github.com/SysFate))                                    |
| Team    | [SysFate](https://www.sysfate.org/)                                                           |
| Email   | <mmendoza@genoscope.cns.fr>                                                                   |

## How the software work
For each sequence:
  - RETIVAD searches for guibson sequences with a regex query. This library allows to define the number of errors to be accepted within the regexpression.
  - Remove the guibson sequence and process each of the side sequences as following:
      - Verify if the size of the sequences are within the defined limits (seqMinLength,seqMaxLength). If it is not the case, classify the sequence fragment as either “too short” or “too long”.
      - In case the length of the sequence is within the expected interval, RETIVAD screens for the cDNA barcodes. In case a unique barcode sequence is retrieved the processing moves to the next step; otherwise it classifies the sequence as “multiple” or “no cDNA sequence found”.
      - Search PCR barcode. In case a unique barcode sequence is retrieved the processing moves to the next step; otherwise it classifies the sequence as “multiple” or “no PCR sequence found”.
      - When both cDNA and PCR barcode sequences are retrieved, the inner sequence is aligned to the pathogen reference genome (in this case the Covid-19 genome).
      - Finally, a variants screening is performed relative to the reference genome.
