## fNIRS1: MATLAB interface for Dr. Johnson’s C suite for dynamic linear models

### Outstanding questions
 - The organization of channels for group analyses: is it always consistent across subjects?
 - How to best format model output for ease of use?
 - Will the software ever need to be run from a PC?
 - In the `*.mat` files Frank shared with us, are the variable names
   standard output from another program? 

##### Prerequisites
 - [FFTW3](http://www.fftw.org/) needs to be installed separately
   - On OSX, this is quite easy with `homebrew` installed: from a
   terminal, the command is simply `brew install fftw`
   - On Ubuntu, for example, open an admin terminal and run `sudo apt-get install fftw`
   - You need to know the location of a file called `fftw3.h`, which comes
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
width="200" height="199"></p> 


Note that that script uses a variable called `FFTW_HOME` which may need
to be set before running `install.m`. On my MacBook, I have it set to
`/usr/local` to reflect the location of `include/fftw3.h`

<p align="center"><img src="vignette/fftw-home.png"
alt="FFTW_HOME variable in install.m"
width="513" height="133"></p> 

If the install script finishes successfully, information will print to
the MATLAB command window, and the directory above `+fnirs1` will be
permanently added to the MATLAB search path so that `fnirs1` commands
will always be accessible. Try running,
```MATLAB
>> help fnirs1.dlm
```
in MATLAB to verify that the installation worked correctly.



### Fitting participant-level models

The primary model fitting tool in our package is a function called
`fnirs1.dlm`. This function takes as input a set of data files and
returns a model summary object. Basic syntax is, 
<pre>
```MATLAB
>> summary = fnirs1.dlm(data_files, <i>Option</i>, <i>OptionValue</i>);
```
</pre>

This function comes with sensible defaults so that the user needs only
input locations of participant files to conduct an
analysis. `fnirs1.dlm` will make copies of the data into a temporary
folder for analysis (marked with the current date and time), and will
run analyses based on setup.dat files written for each channel within
this folder. If the data files argument is left as an empty string,
MATLAB will open a UI window and ask the user to select files by hand,
```MATLAB
>> summary = fnirs1.dlm('');
```

<p align="center"><img src="vignette/file-selection-window.png"
alt="file selection window"
width="438" height="240"></p>

For a single-subject analysis over all fNIRS channels, select only one
participant file. This file should contain variables named `hbo`, `hbr`,
`hbt`, `s`, and `t`. Any optional model parameters that can also be
set/tweaked can also be set via this function. For example, a complete
subject data analysis might stem from the command,
```MATLAB
>> summary = fnirs1.dlm('', 'DownSampleRate', 10);
```

where `'DownSampleRate'` instructs the program to down-sample the data
by a factor of 10. The program will error if for example 10 is not an
integer factor of the participant’s original sampling rate (which was
50 Hz in this case). Optional parameters are handled internally by a
function called `fnirs1.specify_model`, and at the time of writing
include, 
- `'DownSampleRate'` - Integer down-sampling factor of participant’s sampling rate
- `'GroupCovariates'` - Numeric matrix of group-level covariates (group analysis only)
- `'McmcControl'` - An `fnirs1.mcmc_control` object to specify MCMC options
- `'OutcomeType'` - Character. Should be one of {`'hbo'`, `'hbr'`, `'hbt'`}
- `'SpecificChannels'` - Integer vector index to specify analyses on
  subsets of channels
  
<p align="center"><img src="vignette/ch-dirs-setup.png"
alt="temporary directory structure"
width="165" height="160"></p>

After parsing the user input and writing the temporary files,
`fnirs1.dlm` fits Dr. Johnson’s dynamic linear models to the data from
each channel separately with Markov Chain Monte Carlo based
methods. Once analyses are complete, the temporaries are compressed
to save space, and an `fnirs1.dlm_summary` object is returned. This
object contains information about the fitted models, and displays
that information in a (hopefully) convenient format. For example,
printing, 
```MATLAB
>> summary

summary = 

	DLM Summary: ch0001
Parameter                 Estimate  Std.Err. 95.0% Cred.Int.  
------------------------------------------------------------
Participant10_hbo Stim        1.340  0.674      (0.52, 2.56) *
Participant10_hbo Stim        0.027  0.474     (-0.61, 0.74) 
Participant10_hbo Stim       -0.372  0.330     (-0.86, 0.17) 
Participant10_hbo Stim        1.386  0.887      (0.17, 2.49) *
Participant10_hbo Stim        0.018  0.512     (-0.84, 0.81) 
Participant10_hbo Stim       -1.317  0.722    (-2.17, -0.13) *
Participant10_hbo Stim        0.771  0.660     (-0.29, 1.71) 
Participant10_hbo Stim        0.302  0.596     (-0.65, 1.20) 
Participant10_hbo Stim       -0.403  0.383     (-0.86, 0.26) 
------------------------------------------------------------

	DLM Summary: ch0002
Parameter                 Estimate  Std.Err. 95.0% Cred.Int.  
------------------------------------------------------------
Participant10_hbo Stim       -1.260  0.267    (-1.50, -0.75) *
...
```

outputs lots of information about parameter estimates for each
channel. The Parameter field is truncated to save space while
printing, but full parameter descriptions can be extracted via,
```MATLAB
>> summary(1).Descriptions

ans =

  9×1 cell array

    {'Participant10_hbo Stim 1: HRF'       }
    {'Participant10_hbo Stim 1: HRF TmpDrv'}
    {'Participant10_hbo Stim 1: HRF TmpDrv'}
    {'Participant10_hbo Stim 2: HRF'       }
    {'Participant10_hbo Stim 2: HRF TmpDrv'}
    {'Participant10_hbo Stim 2: HRF TmpDrv'}
    {'Participant10_hbo Stim 3: HRF'       }
    {'Participant10_hbo Stim 3: HRF TmpDrv'}
    {'Participant10_hbo Stim 3: HRF TmpDrv'}
```

Similarly, the coefficient estimates, their standard errors, and the
posterior credible intervals can be extracted by accessing the
`Estimates`, `StdErrors`, and `Intervals` fields, respectively. 

We recommend running abbreviated test-analyses to make sure the models
will run before proceeding with full analyses. fnirs1 includes a
convenience function, `fnirs1.mcmc_debug` to help specify short MCMC
chains, and the `'SpecificChannels'` option to reduce data
processing. Users can run a shortened analysis, for example with the
command,
```MATLAB
>> summary = fnirs1.dlm(‘’, ‘DownSampleRate’, 10, ... 
‘McmcControl’, fnirs1.mcmc_debug, ‘SpecificChannels’ = [1:4]);
```

If this command runs smoothly, then we suggest running full analyses
with longer MCMC chains. The default MCMC options run for 20,000
iterations, which should be sufficient for most analyses,
```MATLAB
>> fnirs1.mcmc_control

ans = 

  mcmc_control with properties:

                burnin: 10000
         expectedKnots: 15
    includeDerivatives: 0
         maxIterations: 20000
```



### Fitting group-level models

A _*cautionary note*_: group-analyses are not completely stable
yet. Save `fnirs1.dlm_summary` objects into `*.mat` files as soon as
analyses are complete. MATLAB will sometimes hang or crash after
running group-level analyses this way, but it’s _usually_ been several
minutes after analyses are complete. I’m not sure what’s going on
here: there may be a MATLAB internal bug to blame, and I will continue
to look into it. (MathWorks emailed me to call it an ``unknown issue.'')
At the time of writing, I’m working primarily with MATLAB R2019a
running on OSX 10.14. 


Specifying group-level analyses is largely similar to single-subject
analyses, but in general, group-analyses take longer to run than
single-subject. If the user selects more than one participant’s data
file, `fnirs1.dlm` will automatically conduct a group-level
analysis. For example, continuing with the short debug chains options,
and with,
```MATLAB
>> files'

ans =

  4×1 cell array
    {'Box/BayesianDataAnalysis/ENMA_Data_Sep32019_N29/Participant1.mat'}
	{'Box/BayesianDataAnalysis/ENMA_Data_Sep32019_N29/Participant2.mat'}
    {'Box/BayesianDataAnalysis/ENMA_Data_Sep32019_N29/Participant3.mat'}
    {'Box/BayesianDataAnalysis/ENMA_Data_Sep32019_N29/Participant4.mat'} 

>> summary = fnirs1.dlm(files, 'McmcControl', fnirs1.mcmc_debug, 'DownSampleRate', 10, 'SpecificChannels', [1:4]);
>> summary

summary = 

	DLM Summary: ch0001
Parameter                 Estimate  Std.Err. 95.0% Cred.Int.  
------------------------------------------------------------
Population Stim 1 HRF o      -0.202  0.218     (-0.54, 0.10) 
Population Stim 1 TmpDr       0.000  0.000      (0.00, 0.00) 
Population Stim 1 TmpDr       0.593  0.977     (-0.74, 2.28) 
Population Stim 2 HRF o      -1.210  3.212     (-6.67, 1.81) 
Population Stim 2 TmpDr      -0.001  2.303     (-3.73, 2.70) 
Population Stim 2 TmpDr      -0.061  2.010     (-2.30, 3.64) 
Population Stim 3 HRF o      -0.295  1.792     (-2.12, 2.39) 
Population Stim 3 TmpDr      -2.282  2.691     (-4.58, 2.52) 
Population Stim 3 TmpDr       1.774  1.597     (-1.08, 3.62) 
Participant1_hbo Stim 1       0.112  0.447     (-0.58, 0.81) 
Participant1_hbo Stim 1      -0.354  0.213    (-0.68, -0.06) *
Participant1_hbo Stim 1      -0.053  0.191     (-0.35, 0.25) 
Participant1_hbo Stim 2       0.803  6.995     (-7.72, 9.15) 
Participant1_hbo Stim 2       0.831  1.680     (-1.59, 3.43) 
Participant1_hbo Stim 2      -0.282  0.763     (-1.28, 0.83) 
Participant1_hbo Stim 3      -0.155  0.332     (-0.62, 0.41) 
Participant1_hbo Stim 3      -0.149  0.352     (-0.79, 0.31) 
Participant1_hbo Stim 3      -0.612  0.356     (-1.03, 0.03) 
Participant2_hbo Stim 1       0.617  0.376      (0.18, 1.26) *
Participant2_hbo Stim 1      -3.042  10.441   (-17.43, 12.18) 
Participant2_hbo Stim 1      -0.200  0.436     (-0.69, 0.58) 
Participant2_hbo Stim 2       0.242  0.420     (-0.17, 1.05) 
Participant2_hbo Stim 2      -0.187  0.373     (-0.65, 0.42) 
Participant2_hbo Stim 2       0.278  0.393     (-0.28, 0.96) 
Participant2_hbo Stim 3       0.145  0.315     (-0.39, 0.52) 
Participant2_hbo Stim 3      -0.114  0.499     (-0.66, 0.68) 
Participant2_hbo Stim 3       0.115  0.402     (-0.37, 0.85) 
Participant3_hbo Stim 1      -0.245  0.669     (-1.18, 0.77) 
Participant3_hbo Stim 1      -2.445  10.347   (-17.36, 12.97) 
Participant3_hbo Stim 1       0.006  0.563     (-0.79, 0.69) 
Participant3_hbo Stim 2       0.030  0.511     (-1.07, 0.48) 
Participant3_hbo Stim 2       0.176  0.445     (-0.63, 0.75) 
Participant3_hbo Stim 2      -0.199  0.371     (-0.71, 0.36) 
Participant3_hbo Stim 3      -0.112  0.262     (-0.49, 0.17) 
Participant3_hbo Stim 3      -0.079  0.414     (-0.71, 0.55) 
Participant3_hbo Stim 3      -0.180  0.361     (-0.62, 0.42) 
Participant4_hbo Stim 1       4.502  6.144    (-5.75, 12.56) 
Participant4_hbo Stim 1       0.304  1.058     (-0.81, 2.22) 
Participant4_hbo Stim 1      -0.594  0.372    (-0.97, -0.12) *
Participant4_hbo Stim 2       0.187  0.334     (-0.33, 0.57) 
Participant4_hbo Stim 2       0.136  0.239     (-0.30, 0.47) 
Participant4_hbo Stim 2      -0.186  0.275     (-0.56, 0.21) 
Participant4_hbo Stim 3      -0.028  0.573     (-1.05, 0.62) 
Participant4_hbo Stim 3      -0.207  0.320     (-0.67, 0.24) 
Participant4_hbo Stim 3      -0.104  0.388     (-0.71, 0.33) 
------------------------------------------------------------
... 
```

Inferential statements for group-level parameters are listed at the
top with the 'Population' prefix; participant-level parameter
estimates are given immediately below. 




