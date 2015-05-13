## MODIS Resampling

This program resamples MODIS data.

Dependencies:

* C++ toolchain
* OpenCV
* HDF4 library

Run `make` to build the program named `modisresam`. Running the program
on a granule will modify the data in-place and add an attribute indicating
it was resampled.
