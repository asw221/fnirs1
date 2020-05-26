% INSTALL.m uses relative paths. Must be run from the installer directory.
%   Only intended for use with OSX/Linux operating systems
%
% Example:
%   >> cd +fnirs1/installer
%   >> run install.m

FFTW_HOME = '/usr/local/include';
FFTW_LIB  = '/usr/local/lib';
% User may need to edit 
% FFTW_HOME:
%   - From a terminal, the command ls [value-of-FFTW_HOME]
%     should show the location of a file called "fftw3.h"
%     If not, FFTW_HOME should be edited above to reflect the proper
%     location of this file
% 
% FFTW_LIB:
%   - Similarly, FFTW_LIB should be set to the location of a file called
%     "libfftw3.{a/so/dylib}"


% -------------------------------------------------------------------------
if (ispc)
    error("Installer intended for OSX/Linux")
end


if (~exist(fullfile(FFTW_HOME, 'fftw3.h'), 'file'))
    error('Could not locate ''fftw3.h'' in ''%s''', FFTW_HOME);
end
if ~(numel(dir(fullfile(FFTW_LIB, 'libfftw*'))) >= 1)
    error('Could not locate ''fftw3.h'' in ''%s''', FFTW_HOME);
end



% --- Setup variables for Makefile contents -------------------------------
CXX = 'g++';
if (ismac)
    CXX = 'clang++';
end

FILES = ['main.cpp mcmc.cpp cholesky.cpp randgen.cpp mybspline.cpp ', ...
    'hrf.cpp knots.cpp statistics.cpp config_info.cpp dlm.cpp ', ...
    'kernel_reg.cpp marg_eta_updates.cpp'];
OBJECTS = strrep(FILES, '.cpp', '.o');

% strrep
CXXFLAGS = '-O2 -Wall';
INC = ['-I', FFTW_HOME];
LIB = ['-L', FFTW_LIB];

LINK = '-lm -lfftw3 -lfftw3_threads';
if (ismac)
    LINK = [LINK, ' -Xpreprocessor -fopenmp -lomp'];
elseif (isunix)
    LINK = [LINK, ' -fopenmp'];
end


% --- Write Makefile and compile code -------------------------------------
cd ../include

makefileid = fopen('Makefile', 'w');
if (makefileid == -1)
    cd ../installer
    error('Cannot create Makefile');
end

fprintf(makefileid, '\nCXX      := %s\n', CXX);
fprintf(makefileid, 'SRC      := %s\n', FILES);
fprintf(makefileid, 'OBJ      := %s\n', OBJECTS);
fprintf(makefileid, 'CXXFLAGS := %s\n', CXXFLAGS);
fprintf(makefileid, 'INC      := %s\n', INC);
fprintf(makefileid, 'LIB      := %s\n', LIB);
fprintf(makefileid, 'LINK     := %s\n\n', LINK);
fprintf(makefileid, 'all: $(OBJ)\n\t$(CXX) $(CXXFLAGS) $(OBJ) $(LIB) $(LINK) -o fnirsdlm\n\n');
fprintf(makefileid, '$(OBJ): $(SRC)\n\t$(CXX) -c $(CXXFLAGS) $(SRC)\n\n');
fprintf(makefileid, '{OBJ}: cholesky.h randgen.h\n\n');
fprintf(makefileid, 'clean:\n\trm -f fnirsdlm *.o\n');
fclose(makefileid);

system('make clean', '-echo');
status = system('make all', '-echo');
if (status ~= 0)
    warning('fnirs1: Compilation had non-zero exit status\n');
end

cd ../../
% If top directory is not already on the MATLAB path, add it and save
% the new path
if (~any(strcmp(pwd, regexp(path, pathsep, 'split'))))
    addpath(pwd)
    savepath(userpath)
end

% --- Cleanup and checks --------------------------------------------------
% Return to original directory
cd ./+fnirs1/installer
clear CXX CXXFLAGS FFTW_HOME FFTW_LIB FILES INC LIB LINK OBJECTS makefileid status


% Print warning message if Statistics and Machine Learning Toolbox is not
% available
if (~fnirs1.sml_toolbox_available)
    warning('Statistics and Machine Learning Toolbox may not be available');
end
