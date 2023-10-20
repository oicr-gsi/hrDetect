#!/bin/bash
set -o nounset
set -o errexit
set -o pipefail

cd $1
find . -regex '.*\.SBS.json$' -exec md5sum {} \;
find . -regex '.*\.ID.json$' -exec md5sum {} \;
ls | sed 's/.*\.//' | sort | uniq -c