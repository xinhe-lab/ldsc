#!/bin/bash

__conda_setup="$('/project2/xinhe/software/miniconda3/bin/conda' 'shell.bash' 'hook' 2> /dev/null)"
if [ $? -eq 0 ]; then
    eval "$__conda_setup"
else
    if [ -f "/project2/xinhe/software/miniconda3/etc/profile.d/conda.sh" ]; then
        . "/project2/xinhe/software/miniconda3/etc/profile.d/conda.sh"
    else
        export PATH="/project2/xinhe/software/miniconda3/bin:$PATH"
    fi
fi
unset __conda_setup
conda activate ptb_r
