 #!/bin/sh

matlab -nodesktop -r "cd /local/charis/myCode/MATLAB/HydroxyMethylation/,\
addpath(genpath(pwd)), chrToRun = $1,\
dataPath = '/local/charis/myCode/MATLAB/HydroxyMethylation/genomeWide/RRHPBS.test.data',\
cd genomeWide, runDataGW(dataPath, chrToRun), exit"
