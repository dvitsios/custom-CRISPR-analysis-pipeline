# WORKFLOW

1. Needleman-Wunsch Alignment

2. `global_parse_needle_out.pl`
	-> `parse_needle_out.pl`

3. `process_basic_stats_and_classification.R`

***

**[Step-1]**:
Needleman-Wunsch Alignment scripts.

*For normal samples*: 
```
bsub -n 4 ./align_all_samples.sh
```

*For calibration samples*:
```
bsub -n 4 ./align_all_calibration.sh
``` 

(*Note: the orphans cannot aligned with NW since I don't know which amplicon they should be aligned against*)
