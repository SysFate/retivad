# retivade
A RealTime Variant Detector software.

## Description
With the covid pandemic, we want to create a method to detect covid variant without use expensive kit.
In this idea we want combine it to know if the patient are covid positive.
Illumina sequencing are to expensive, this is why our choise was a Nanopore sequencer: https://nanoporetech.com/
This sequencer have the capacity to give in realtime the sequenced sequence.
With the low error rate of the sequencer (2 errors for 24 bp), we are able to identify patient and covid variant.

## Conda version
### Dependencies

```
conda create --name retivade paramiko medaka minimap2 regex bwa samtools bedtools biopython bedops matplotlib
conda activate retivade
```

### How run

```
require arguments:
  -5 cDNAbc, --cDNAbc cDNAbc
                           File with the barcode cDNA' list. ID Sequence
  -3 PCRbc, --PCRbc PCRbc
                           File with the barcode PCR' list. ID Sequence
  -g genome, --genome genome
                           Aligner genome index path
  ```
```
optional arguments:
  -h, --help               show this help message and exit
  -f file, --file file     File(s) you want annalyse separate with commat (default: )
  --folder folder          Folder who contain fastq to annalyse (default: None)
  --tmp tmp                Folder who contain temporar files (default: /tmp)
  -c core, --core core     number of core the software can use (default: 1)
  -o outdir, --out outdir  Folder where send files (default: .)
  -n name, --name name     Name annalyse. By default the software use the begin date. (default: None)
  -G guibson, --guibsonE guibson
                           Guibson error accepted (default: 4)
  -B barcode, --barcodeE barcode
                           Barcode error accepted (default: 4)
  -i insertMinLength, --insertMinLength insertMinLength
                           Good insert minimum length (default: 57)
  -I insertMaxLength, --insertMaxLength insertMaxLength
                           Good insert maximum length (default: 257)
  -l seqMinLength, --seqMinLength seqMinLength
                           Good sequences minimum length between guibson sequences (default: 107)
  -L seqMaxLength, --seqMaxLength seqMaxLength
                           Good sequences maximum length between guibson sequences (default: 357)
  --noreverse3             If you send directly the reverse complement in 3' barcodes, activate this option (default: False)
  -a aligner, --aligner aligner
                           Aligner name (default: bwa)
  --alignerPath alignerPath
                           Aligner path. By default we use the name (default:None)
  --alignerSpecialsParam alignerSpecialsParam
                           Specials parameters for the aligner (default: None)
  --traceseq traceseq      Sequences you want track to check the method (default: None)
  --ext ext                Path to the configuration for the ssh. Use it if your files are not in the computer.
                           This method are only in realtime annalyse. (default: None)
  --timecheck timecheck    In realtime annalyse, you can setup the time between each check in the folder (default: 10)
  --timeEndCheck timeEndCheck
                           A realtime annalyse, set up the time between the last good check and the stop of the software (default: 30)
  --dev                    Keep all temporar datas (default: False)
  ```

#### Example
In this example we have a sequencing who begin in a external device (eg: mk1-c), and we want annalyse the data in realtime.
Take files from the template folder:
- Barcodes datas files: barcode5.tsv (cDNA barcodes) and barcode3.tsv (PCR Barcodes).
- The genome of your interest pathogene(s). Here we select the covid genome.
- The configSSH file is used when your files are not in your computer. (Work only in realtime annalyse)
- fastq_test sequences extract from the patient. (The sample come from synthetic DNA, not from real patient !)

##### Send the folder fastq_test in the external computer.
##### Use this command to launch the software:
```
python3 selectSequencesWithBarcordesLittle.py -c 8 -5 barcode5.tsv -3 barcode3.tsv -g covid.fa --folder fastq_test -n covidTest --ext configSSH.ini > covid.log 2>&1 &
```
##### Move the file in fasq_test/inTime in the fastq_test folder
##### Wait and see
In realtime the software give plot with the curents alignment stats and a summary.
If the software don't find new files after 30min, the software initialise the end.
After the end of the software you obtain a summary with the essentials datas.

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
  - Detect guibsons sequences in the original sequence. We use the lib regex. This library allow in the regexpression a number of errors.
  - Cut the guibson, and separate each part.
  - For each part:
      - Continue if the sequence have length between limit (seqMinLength,seqMaxLength).
      - Search cDNA barcode. Continue if a unique barcode find.
      - Search PCR barcode. Continue if a unique barcode find.
      - Continue if the cDNA and PCR barcode create a good construct.
      - Continue if the insert between the two barcode have a length between limit (insertMinLength,insertMaxLength).
      - Align the insert
 For each couple of barcode find:
  - Search variation beetween seq
