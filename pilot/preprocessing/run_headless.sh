#!/bin/bash
# Wrapper script to run scripts without opening figures

env -u DISPLAY matlab -nodisplay -nosplash -nodesktop -r "run('a7_run_all.m'); exit"
