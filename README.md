# retivade

## Dependencies

```conda create --name retivade paramiko medaka minimap2 regex bwa samtools bedtools biopython bedops matplotlib
conda activate retivade
```

## How run



### Example
```python3 selectSequencesWithBarcordesLittle.py -c 10 -5 barcode5.csv -3 barcode3.csv -g genome/covid.fa --folder /data/covid12/no_sample/20210317_1617_MC-110352_0_FAP14908_252839db/fastq_pass -i 30 -I 10000 -l 70 -L 10000 -o /home/Sysfate/Team_projects/Covid19/covidInRealTime --ext config.ini > covid.log 2>&1 &
```

## Informations / Authors

|         |                                                                                               |
| ------- | --------------------------------------------------------------------------------------------- |
| Author  | Francois STUDER ([Github](https://github.com/studyfranco))                                    |
| Author  | MENDOZA PARRA Marco ([Github](https://github.com/SysFate))                                    |
| Team    | [SysFate](https://www.sysfate.org/)                                                           |
| Email   | <mmendoza@genoscope.cns.fr>                                                                   |

