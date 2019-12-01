#!/bin/sh
python3 ./run_optimization.py --config=./setup_files/cfg_test_det_US --supply=./setup_files/US_vax_coverage_2014-15_weekly_start210.csv \
 --output=output_optimization --step_size=100000 --temporal_search_size=1000 --look_ahead=2 --log=optimization.log
