# CHANGELOG
## 1.5.6 -2023-11-17
Included new version of sigtools-rscript module that fixed bugs
## 1.5.2 -2023-10-25
Change workflow name
## 1.5.1 -2023-10-24
Update to call call_hrdetect.R from sigTools_runthrough.R
## 1.5 -2023-10-18
Finalized version ready to release
## 1.4 - 2022-07-26
- Fixed bugs
	- removed rare fit set to Pancreas
- Split INDEL and SNV cutoffs
- Protects from crash on tissue not in catalog
- Protects from crash on low SNV/INDEL countsS
- New function to extend LOH (default: OFF)
- Added CNV signatures (default: OFF)

## 1.3 - 2022-06-22
- Fixed SV reformatting for bug on end2 (was printing end1)
- Added CHORD to runthrough
- changed runthrough script output to make JSON and save more information
- changed output for plotIt from .png to .svg

## 1.1 - 2022-05-11
- Changed name of workflow to HRDetect
- Changed plotting from .pdf to .png, and changes intake format
- Annotated HRDetect R results script, added loop for low indel results and added genomeVersion argument
- merged SNV and INDEL filtering using alias in wdl
- added filtering option parameters (instead of hardcoded)

## 1.0 - 2022-04-11
- A brand-new workflow.