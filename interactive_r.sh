#!/bin/bash

# Interactive shell for R
salloc -J intr --mem=16GB --time=6:00:00 --partition=nodes --ntasks=4

srun --pty /bin/bash
