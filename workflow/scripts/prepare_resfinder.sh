#!/usr/bin/env bash

mkdir -p resources && cd resources
git clone https://git@bitbucket.org/genomicepidemiology/resfinder_db.git
echo "Downloaded ResFinder Database on:"
date
ls -lh resfinder_db
echo
echo "The version of the ResFinder DB is:"
cat resfinder_db/VERSION
echo
