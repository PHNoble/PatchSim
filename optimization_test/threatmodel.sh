#!/bin/sh
python3 ./run_threat_model.py --config=./setup_files/cfg_threat_US --supply=./setup_files/US_vax_coverage_2014-15_weekly_start210.csv \
 --output=output_threatmodel --log=threat.log --verbose
