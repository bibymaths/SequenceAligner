"""Main entry point for the alignment analysis tool.

This module exposes a simple wrapper around the CLI's ``main`` function
allowing the package to be executed with ``python -m alignment_tool``.
"""

from __future__ import annotations

import sys

from . import cli


def main() -> None:
    """Invoke the CLI's main function and exit with its return code."""
    code = cli.main()
    sys.exit(code)


if __name__ == '__main__':
    main()
