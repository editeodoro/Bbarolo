
.. _fitsutilities:

Tools for FITS files
####################

Sometimes a FITS file does not conform with BBarolo or keywords in the FITS header can not be understood by the code. A few utilities to easily modify a FITS file and its header are now included in the latest versions of BBarolo.
These tools are based on or developed from the  
`NASA's FITS utilities <https://heasarc.gsfc.nasa.gov/docs/software/fitsio/cexamples.html>`_.

A list of implemented tools can be obtained with ``BBarolo --fitsutils``. Four utilities are currently available. Each tool is promptly accessible from the command line through ``BBarolo --UTILNAME``, where UTILNAME is one of the utilities listed below:

**modhead**. Modify or add a keyword to a FITS header. Type ``BBarolo --modhead`` for help:

.. literalinclude:: tasks/examples/modhead.txt

|

**remhead**. Remove a keyword from a FITS header. Type ``BBarolo --remhead`` for help:

.. literalinclude:: tasks/examples/remhead.txt

|

**listhead**. List all keywords in a FITS header. Type ``BBarolo --listhead`` for help.

.. literalinclude:: tasks/examples/listhead.txt

|

**fitscopy**. Make a copy of a FITS file, optionally selecting a subset. Type ``BBarolo --fitscopy`` for help:

.. literalinclude:: tasks/examples/fitscopy.txt

|

**fitsarith**. Perform arithmetic operations between two FITS files or on a single FITS file. Type ``BBarolo --fitsarith`` for help:

.. literalinclude:: tasks/examples/fitsarith.txt

