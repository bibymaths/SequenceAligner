"""
Benchmark package initializer.

This package provides a modular framework for benchmarking sequence alignment tools.
Each tool has its own runner located in ``benchmark/runners`` and parsers in
``benchmark/parsers`` provide functionality to interpret alignment outputs. The
``benchmark.py`` script orchestrates the benchmarking workflow using a YAML
configuration file and writes structured results to CSV and JSON outputs.
"""
