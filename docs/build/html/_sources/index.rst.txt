VCFclass: VCF interface for python
==================================

:Author: Joseph Lalli
:Date: June 2nd, 2021
:Version: 0.1.0

VCFclass is a python module that allows for better integration of
within-sample variation data stored in VCF files and common bioinformatics 
tools (Pandas, numpy, etc.) As currently designed, it is memory inefficient
and is not suitable for populations with references sequences in the range
of megabases. However, VCFclass works well for handling viral sequencing data.

VCFclass is best used to easily load and convert VCF files into pandas dataframes,
and/or "Popoolation" style numpy arrays of read counts. It is also excellent for
averaging technical sequencing replicates, and has some iSNV filters built in.

In the future VCFclass will implement pysam's VCF functions in the background.
This should substantially improve performance. VCFclass was originally written
to provide for functionalities that pysam did not provide (especially VCF averaging,
VCF filtering, and VCF export as a dataframe and as a numpy array).

To install the latest release, type::

    pip install VCFclass

-Note: this is not currently implemented. Until this project reaches v0.2, please 
simply clone this repository or place VCFclass.py in your working directory.

See the :ref:`Installation notes <installation>` for details.

Contents
--------

.. toctree::
   :maxdepth: 2

   api.rst
   usage.rst
   installation.rst

Indices and tables
------------------

Contents:

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`