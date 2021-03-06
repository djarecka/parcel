#parcel

The parcel model represents an idealised scenario of a 0-dimensional 
  air parcel rising adiabatically with a constant vertical velocity. 

The representation of microphysical and chemical processes in the parcel model 
  is done using the particle-based scheme from the [libcloudph++](http://libcloudphxx.igf.fuw.edu.pl/) 
  library of cloud microphysics schemes.
For more information on the [libcloudph++](http://libcloudphxx.igf.fuw.edu.pl/) library, 
  its' dependencies, installation and source code please check the 
  [project repository](https://github.com/igfuw/libcloudphxx) and 
  the [project documentation](http://www.geosci-model-dev.net/8/1677/2015/).

The parcel model is written in Python 2.7. 

# installation

The .travis.yml file shipped with the parcel model 
  contains a complete set of commands needed to execute all test programs
  shipped with the model on fresh Ubuntu and OSX installations -
  it may contain useful information on obtaining the dependencies.

# testing

The parcel model is shipped with a set of tests for the particle-based scheme 
  from the [libcloudph++](http://libcloudphxx.igf.fuw.edu.pl/) library.
For test automation py.test Python package is used.

To run all the tests  please type in terminal

  $ py.test unit_test/

  $ py.test long_test/

Some tests generate plots. 
These plots are saved in /plots/outputs folder.
