#! /usr/bin/python
# -*- coding: utf-8 -*-

"""Deprecated: use ``3_plot_lightcurve.py`` with ``tables/light_curves.ecsv``."""

import sys

print(
    "4_plot_lightcurve_from_ecsv.py is retired.\n"
    "Plot from the pipeline table with 3_plot_lightcurve.py "
    "(reads output/tables/light_curves.ecsv).",
    file=sys.stderr,
)
sys.exit(1)
