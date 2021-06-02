======================================================
VCFclass - An interface for interacting with VCF files
======================================================

Introduction
============

VCFclass is a python module that allows for better integration of
within-sample variation data stored in VCF files and common bioinformatics 
tools (Pandas, numpy, etc.) VCFclass is designed to handle the small VCF 
files typical of viral sequencing data.

This page provides a quick introduction in using VCFclass followed by the
API. See :ref:`usage` for walkthroughs of typical usage cases.

To use the module to read a file in BAM format, create a
:class:`~VCFclass.VCF` object:: python

   import VCFclass
   my_vcf = VCFclass.VCF("my_vcf.vcf")

Once a file is opened you can directly iterate over mutation calls
in the VCF file:: python

    for mut in my_vcf:
        do_something

Each iteration returns a :class:`~VCFclass.MutCall` object.

An alternative way of accessing the data in a VCF file is by
calling :meth:`~VCFclass.VCF.to_dataframe` method. This function
returns a dataframe containing one row per sample per mutation.
Each row contains information about the mutation call (position, 
ref, alt, etc.) along with information specific to each sample
(depth, ref_depth, alt_depth, etc.) This information is particularly
easy to use with downstream tools (e.g., Seaborn, or Pandas' groupby
functions.)

More detailed usage instructions, including information about how to
average VCFs that are technical replicates and to apply a read position
SNP filter, are at :ref:`usage`.

.. note::

   Coordinates in VCFclass are always 0-based (following the python
   convention). Text VCF files, SAM files, and GTF files use 1-based 
   coordinates. VCFclass converts positions from 1-based coordinates 
   to 0-based upon import, and converts positions back to 1-based 
   coordinates upon export to text VCF files. Pandas Dataframes created 
   using the :meth:`~VCFclass.to_dataframe` method will have 0-based 
   coordinates.

   In general, VCFclass' design philosophy is to work in 0-based coordinates
   at all times unless a format specifies 1-based coordinates.

API
===

An :class:`~VCFclass.VCF` represents a whole VCF file.

.. autoclass:: VCFclass.VCF
   :members:

An :class:`~VCFclass.MutCall` represents one row of a VCF file.

.. autoclass:: VCFclass.MutCall
   :members:

An :class:`~VCFclass.SampMut` represents information about a 
mutation that is specific to an individual sample.

.. autoclass:: VCFclass.Header
   :members:
