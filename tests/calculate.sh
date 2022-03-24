#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

cd $1

find . -regex '.*\.sigtools.hrd.txt$' -exec md5sum {} \;
