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
  --ext ext                Configuration for the ssh. Use it if your file are not in the computer.
                           This method are only in realtime annalyse. (default: None)
  --timecheck timecheck    A realtime annalyse, you can setup the time between each check in the folder (default: None)
  --dev                    Keep all temporar datas (default: False)
  ```

### Example
```
python3 selectSequencesWithBarcordesLittle.py -c 10 -5 barcode5.csv -3 barcode3.csv -g genome/covid.fa --folder /data/covid12/no_sample/20210317_1617_MC-110352_0_FAP14908_252839db/fastq_pass -i 30 -I 10000 -l 70 -L 10000 -o /home/Sysfate/Team_projects/Covid19/covidInRealTime --ext config.ini > covid.log 2>&1 &
```

## Informations / Authors

|         |                                                                                               |
| ------- | --------------------------------------------------------------------------------------------- |
| Author  | Francois STUDER ([Github](https://github.com/studyfranco))                                    |
| Author  | MENDOZA PARRA Marco ([Github](https://github.com/SysFate))                                    |
| Team    | [SysFate](https://www.sysfate.org/)                                                           |
| Email   | <mmendoza@genoscope.cns.fr>                                                                   |

