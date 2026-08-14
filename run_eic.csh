#!/bin/tcsh

setenv EIC_SHELL_PREFIX '/eic/u/aranya/eic/local/'
setenv SINGULARITY_BINDPATH '/gpfs02,/gpfs01,/gpfs'
/usr/bin/singularity exec $EIC_SHELL_PREFIX/lib/eic_xl-nightly /bin/bash -c "$argv[1-]"
