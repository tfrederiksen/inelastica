.. _tshs:

TSHS
====

The following is provode by Guangping Zhang. Something here may have mistakes. Just use it for reference.

In order to see the content in *.TSHS* file, the most easy way to do it is modify the subroutine ``m_ts_io.F90`` which produce the *.TSHS* file.

``m_ts_io.F90`` subroutine consist of two parts: read (first part) and write (second part).

We can replace

.. code-block:: bash

    open( iu, file=fname, form='unformatted', status='unknown' )

by

.. code-block:: bash

    open( iu, file=fname, form='formatted', status='unknown' )

and all the

.. code-block:: bash

    write(iu)

by

.. code-block:: bash

    write(iu,*)

in the second part. Recompile the `TranSIESTA`_ by shootting ``make ts`` (versions 4.0-4.1) or ``make`` (for later versions) in the `Obj` directory (or any directory you want to build the `TranSIESTA`_, and maybe a ``make clean`` may need before if you ever build `TranSIESTA`_ in the same directory).

Now use the new executable transiesta to run the task to see the content of *.TSHS* file.

*NOTE:* the length is in Bohr and energy and the temperature are in Rydberg.
