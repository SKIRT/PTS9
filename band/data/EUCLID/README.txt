2021-12-06 - Crescenzo Tortora - crescenzo.tortora@inaf.it

For the moment, we have decided to consider as reference filters the 4 of Euclid and the Rubin ones. These will be used
for creating stamps, and to convolve templates used in the pixel-based SED fitting procedure.

The filter files stored in this folder have been downloaded from here:

https://colab.research.google.com/drive/17m9z1DYLrKZBDzn4mMKCUN1fGYAHRIHS?usp=sharing

a page prepared by Micol Bolzonella, and shared with us, for the aims of our WP (to my knowledge they would be the most
updated, at the moment).


2021-12-09 - Micol Bolzonella - micol.bolzonella@inaf.it

The values we have tabulated are the ones for photon counters, so R(lambda), and in fact I usually multiply them by
lambda.

For LSST, you can find some description here: https://github.com/lsst/throughputs/tree/main/baseline

In general I use the same definitions implemented in the SED fitting code lephare
https://www.cfht.hawaii.edu/~arnouts/LEPHARE/DOWNLOAD/lephare_doc.pdf
pagg. 8-9 and 11-12. At first glance they seem consistent with definitions in the paper Camps et al. 2016.
