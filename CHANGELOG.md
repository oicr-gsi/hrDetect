# Changelog
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).


## [1.7.0] - 2024-06-25
### Added
- [GRD-797](https://jira.oicr.on.ca/browse/GRD-797) - add vidarr labels to outputs (changes to medata only).

## [1.6.0] - 2024-02-02
### Changed
- Changes to module name.

## [1.5.6] -2023-11-17
### Added
- Included new version of sigtools-rscript module that fixed bugs.

## [1.5.2] -2023-10-25
### Changed
- Change workflow name.

## [1.5.1] -2023-10-24
### Changed
- Update to call call_hrdetect.R from sigTools_runthrough.R

## [1.5.0] -2023-10-18
### Changed
- Finalized version ready to release.

## [1.4.0] - 2022-07-26
### Added
- New function to extend LOH (default: OFF).
- Added CNV signatures (default: OFF).

### Fixed
- Fixed bugs

### Changed
- Removed rare fit set to Pancreas.
- Split INDEL and SNV cutoffs.
- Protects from crash on tissue not in catalog.
- Protects from crash on low SNV/INDEL counts.

## [1.3.0] - 2022-06-22
### Added
- Added CHORD to runthrough

### Changed
- Fixed SV reformatting for bug on end2 (was printing end1).
- changed runthrough script output to make JSON and save more information.
- changed output for plotIt from .png to .svg

## [1.1.0] - 2022-05-11
### Added
- Annotated HRDetect R results script, added loop for low indel results and added genomeVersion argument.
- Added filtering option parameters (instead of hardcoded).

### Changed
- Changed name of workflow to HRDetect
- Changed plotting from .pdf to .png, and changes intake format
- Merged SNV and INDEL filtering using alias in wdl

## [1.0.0] - 2022-04-11
### Added
- A brand-new workflow.
