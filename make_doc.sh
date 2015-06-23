#!/usr/bin/env bash

cd code
matlab -nosplash -nodesktop -nodisplay -r "paths; MatlabDocMaker.create; exit;"

