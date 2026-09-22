"""`python -m rnaseq_pipeline` — same as the `rnaseq-pipeline` command."""
import sys

from .cli import main

sys.exit(main())
