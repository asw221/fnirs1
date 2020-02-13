# fNIRS1: MATLAB interface for Dr. Johnson’s C suite for dynamic linear models

### Outstanding questions
 - The organization of channels for group analyses: is it always consistent across subjects?
 - How to best format model output for ease of use?
 - Will the software ever need to be run from a PC?
 - In the `*.mat` files Frank shared with us, are the variable names
   standard output from another program? 

##### Prerequisites
 - [FFTW3](http://www.fftw.org/) needs to be installed separately
 - On OSX, this is quite easy with homebrew installed: from a
 terminal, the command is simply `brew install fftw`
 - On Ubuntu, for example, open an admin terminal and run `sudo apt-get install fftw`
You need to know the location of a file called `fftw3.h`, which comes
with FFTW3. On my MacBook, this file is found in
`/usr/local/include/fftw3.h`

### Installation
As part of `fnirs1`, we provide an install script for convenient
installation on unix-based operating systems. To install `fnirs1`, copy
the `+fnirs1` folder and all its contents to any available
location. Navigate to `+fnris1/installer/` and at the MATLAB prompt,
execute, 
```MATLAB
>> run install.m
```

<p align="center"><img src="vignette/install-script-location.png"
alt="install script location"
width="120" height="281"></p> 


