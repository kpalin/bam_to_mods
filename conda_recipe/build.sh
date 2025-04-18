#!/bin/bash

#  Created on Friday, 29 April  2022

# For debugging
#set -o verbose

# Die on unset variables
set -o nounset
# Die on errors
set -o errexit
# Die if any part of a pipe fails
set -o pipefail

make release debug
cp release/bam_to_mods "${PREFIX}/bin/"
cp debug/bam_to_mods "${PREFIX}/bin/bam_to_mods_dbg"
