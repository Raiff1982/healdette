# Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

## [Unreleased]

### Added
- Celtic-specific binding motif validation (WY, RF, KW, YF)
- Enhanced sequence validation parameters
- Charge distribution analysis
- Aromatic content optimization
- Population-specific HLA validation

### Changed
- Optimized validation thresholds for Celtic-specific sequences
- Improved hydrophobic content requirements (35-45%)
- Refined aromatic content targets (15-27%)
- Adjusted charge balance parameters (+5 to +15)

### Fixed
- Sequence validation to prevent homopolymer runs
- Charge distribution validation
- Structure metrics calculation

## [1.0.0] - 2025-09-26
### Added
- Initial release of Healdette pipeline
- ProtGPT2 integration for sequence generation
- Basic sequence validation
- Biophysical property analysis
- Population-specific immunogenicity assessment
- JSON output format
- Framework template-based generation