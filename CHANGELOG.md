# Changelog
All notable changes to this workflow will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [1.0.0] - 2026-06-19
### Added
- Initial WDL workflow
- task  merge by partition, across lanes. Partitions are chr1..chr22, chrX, chrY, chrM each on their own, plus one OTHER partition holding all remaining contigs and the unmapped reads.
- task mark duplicates within each partition (RAM scaled by the partition's size).
- task merge the duplicate-marked partitions into the final sample CRAM.
- tasks collecting all the QC metrics
