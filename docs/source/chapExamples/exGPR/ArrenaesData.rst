Reference data set from Arrenæs
===============================

A 2D/3D data set of recorded travel time data from a cross hole Georadar
experiment is available in the 'data/crosshole' folder.

4 Boreholes were drilled, AM1, AM2, AM3, and AM4 at the locations shown
below

.. figure:: ../../figures/sippi_arrenaes_3d_setup.png
   :alt: Location of boreholes AM1, AM2, AM3, and AM4 at Arrenæs.

   Location of boreholes AM1, AM2, AM3, and AM4 at Arrenæs.

Travel time data were collected between boreholes AM1 and AM3, and AM2
and AM4 respectively, in a depth interval between 1m and 12m. The travel
times for each of the two 2D data sets are available in the
`AM13\_data.mat <#>`__ and `AM24\_data.mat <#>`__ files. All the data
have been combined in the 3D data set available in
`AM1234\_data.mat <#>`__.

All mat-files contains the following variable

::

    S --> [ndata,ndim] each row contains the position of the source
    R --> [ndata,ndim] each row contains the position of the receiver
    d_obs --> [ndata,1] each row contains the observed travel time in milliseconds 
    d_std --> [ndata,1] each row contains the standard deviation of the uncertainty of the observed travel time in milliseconds

All data are also available as ASCII formatted EAS files in
``AM13_data.eas``, ``AM24_data.eas``, and ``AM1234_data.eas``.

The following 3 Figures show the ray coverage (using straight rays) for
each of the AM13, AM24, and AM1234 data sets. The color of each ray
indicates the average velocity along the ray computed using v\_av =
raylength/d\_obs. AM13 ray coverageAM24 ray coverageAM1234 ray coverage.

.. figure:: ../../figures/arrenaes_raycoverage.png
   :alt: Ray coverage between wells left) AM1-AM3, middle) AM2-AM4,
   right)AM1-4.

   Ray coverage between wells left) AM1-AM3, middle) AM2-AM4,
   right)AM1-4.
