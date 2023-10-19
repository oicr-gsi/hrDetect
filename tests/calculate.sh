#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

cd $1
find . -name '*.json' -xtype f -exec sh -c "cat {} | md5sum" \;
ls | sed 's/.*\.//' | sort | uniq -c
