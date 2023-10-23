.. _troubleshooting:

Troubleshooting
###############

I will try to write a troubleshooting page while I receive the feedback from BBarolo's users. 

In the meanwhile, please report any bug or problem you have with BBarolo and pyBBarolo. 
If you are a Github user, you can submit an issue ticket at `this page <https://github.com/editeodoro/bbarolo/issues>`_. Otherwise you can `email me <enrico.diteodoro@gmail.com>`_. Please attach any significant error messages and tell me how to reproduce the problem. 

I will try to fix the issues as soon as I can. Thank you.


Frequently Asked Questions
^^^^^^^^^^^^^^^^^^^^^^^^^^

.. container:: toggle 

    .. container:: header

        **1) The code compiled successfully but when I run BBarolo I get an "error while loading shared libraries: libXXX: cannot open shared object file: No such file or directory"**. 
       
    This error arises when some of BBarolo's dependencies (CFITSIO, WCS or FFTW libraries) have not been installed in a default library path. Solve this problem by adding the library path to the LD_LIBRARY_PATH variable. If you are using Bash shell::

        export LD_LIBRARY_PATH=/path/to/libXXX:$LD_LIBRARY_PATH

    If you are using C-shell::

        setenv LD_LIBRARY_PATH /path/to/libXXX:$LD_LIBRARY_PATH

    To make it permanent, the above commands can be simply added to your ~/.bashrc or ~/.cshrc files (make sure to start a new shell after doing it).

|

.. container:: toggle

    .. container:: header

        **2) On my Mac OS, BBarolo compiles correctly but then it crashes with no error message**. 
       
    This may be due to Apple default compiler, Clang. The latest versions of Clang++ does not compile the code correctly, resulting in a SEGSEV error (like "zsh: abort"). To fix this, recompile the code using another compile (for example GNU GCC).

|

.. container:: toggle

    .. container:: header

        **3) BBarolo seems to run smoothly but suddenly it gets "Killed"**. 
       
    The code has been killed by the system kernel for some reason. For example, you can       check that with::

        dmesg | grep -i kill

    99% will be a memory problem. BBarolo needs to allocate as much memory as three times     the size of the input datacube. If you are working with big datasets (>4 GB), you can     easily run out of memory. Please use a more powerful machine.

