#!/usr/bin/env bash

cd data
git clone https://git@bitbucket.org/genomicepidemiology/resfinder_db.git
echo "Downloaded ResFinder Database on:"
date
ls -lh resfinder_db
echo
echo "The version of the ResFinder DB is:"
cat resfinder_db/VERSION
echo

echo "==================================="
echo "NOTE: manual intervention required!"
echo "==================================="
echo "Make sure Python 3 and KMA are installed, or activate the right conda"
echo "environment, e.g., `conda activate kma` and install the database as"
echo "follows: `$ cd data/resfinder_db && python3 INSTALL.py`"
