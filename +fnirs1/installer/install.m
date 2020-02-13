% INSTALL.m uses relative paths. Must be run from the installer directory.
%   Only intended for use with OSX/Linux operating systems
%
% Example:
%   >> cd +fnirs1/installer
%   >> run install.m

FFTW_HOME = '/usr/local';
% User may need to edit FFTW_HOME:
%   - From a terminal, the command ls [value-of-FFTW_HOME]/include
%     should show the location of a file called "fftw3.h"
%     If not, FFTW_HOME should be edited above to reflect the proper
%     location of this file

if (ispc)
    error("Installer intended for OSX/Linux")
end

FILES = [{'include/mcmc.o '}; ...
    {'include/cholesky.o '}; ...
    {'include/randgen.o '}; ...
    {'include/mybspline.o '}; ...
    {'include/hrf.o '}; ...
    {'include/knots.o '}; ...
    {'include/statistics.o '}; ...
    {'include/config_info.o '}; ...
    {'include/dlm.o '}; ...
    {'include/kernel_reg.o '}; ...
    {'fitDlm.cpp'} ];
INC = ['-I./include -I', fullfile(FFTW_HOME, 'include')];
LIB = [' -L', fullfile(FFTW_HOME, 'lib'), ' -lm -lfftw3 '];

% Compile object code
cd ../include
system('make clean', '-echo');
system('make mcmc.o', '-echo');
% Compile MATLAB interface
cd ../
mex -setup cpp
script = ['mex ', INC, LIB, ...
    '-v CXXFLAGS=''$CXXFLAGS -O2 -Wall -std=c++11'' ', FILES{:}];
evalin('caller', script);
% mex('-v', 'CXXFLAGS=''$CXXFLAGS -O2 -Wall -std=c++11'' ', ...
%     INC, LIB, FILES{:});
% mex -I./include/ -I/usr/local/include/ -L/usr/local/lib/ -v CXXFLAGS='$CXXFLAGS -O2 -Wall -std=c++11' -lm -lfftw3 ./include/mcmc.o ./include/cholesky.o ./include/randgen.o ./include/mybspline.o ./include/hrf.o ./include/knots.o ./include/statistics.o ./include/config_info.o ./include/dlm.o ./include/kernel_reg.o fitDlm.cpp
% If top directory is not already on the MATLAB path, add it and save
% the new path
cd ../
if (~any(strcmp(pwd, regexp(path, pathsep, 'split'))))
    addpath(pwd)
    savepath
end
% Return to original directory
cd ./+fnirs1/installer
clear FFTW_HOME FILES INC LIB script
