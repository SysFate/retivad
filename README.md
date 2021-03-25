# retivade

## Dependencies

```
conda create --name retivade paramiko medaka minimap2 regex bwa samtools bedtools biopython bedops matplotlib
conda activate retivade
```

## How run

```
require arguments:
  -5 barcode5, --barcode5 barcode5
                           File with the barcode 5' list. ID Sequence
  -3 barcode3, --barcode3 barcode3
                           File with the barcode 3' list. ID Sequence
  -g genome, --genome genome
                           Aligner genome index path
  ```
```
optional arguments:
  -h, --help               show this help message and exit
  -f file, --file file     File(s) you want annalyse separate with commat (default: )
  --folder folder          Folder who contain fastq to annalyse (default: None)
  -c core, --core core     number of core the software can use (default: 1)
  -o outdir, --out outdir  Folder where send files (default: .)
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
  --dev                    Keep all temporar datas (default: False)
  ```

### Example
Take files from the template folder:
- Barcodes datas files: barcode5.tsv (RT-PCR barcodes) and barcode3.tsv (PCR Barcodes).
- The config file is used when your files are not in your computer.
- The genome of your interest pathogene(s). Here we select the covid genome.

We want in this example, diagnostic patient with covid.
```
python3 selectSequencesWithBarcordesLittle.py -c 8 -5 barcode5.tsv -3 barcode3.tsv -g covid.fa --folder fastq_test -n covidTest --ext config.ini > covid.log 2>&1 &
```
In realtime the software give plot with the curents alignment stats and a summary.
After the end of the software you obtain a summary with the essentials datas.

## Informations / Authors

|         |                                                                                               |
| ------- | --------------------------------------------------------------------------------------------- |
| Author  | Francois STUDER ([Github](https://github.com/studyfranco))                                    |
| Author  | MENDOZA PARRA Marco ([Github](https://github.com/SysFate))                                    |
| Team    | [SysFate](https://www.sysfate.org/)                                                           |
| Email   | <mmendoza@genoscope.cns.fr>                                                                   |

