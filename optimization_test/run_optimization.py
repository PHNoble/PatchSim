#!/usr/bin/env python
# -*- coding: utf-8 -*-
import logging
import argparse
import optimization as opt
import simulation as sim


log = logging.getLogger("Logger1")

def parse_args():
    ap = argparse.ArgumentParser(description='Find optimal vax assignment schedule.')
    ap.add_argument('-c', '--config', required=True,
                    help='simulation configuration file')
    ap.add_argument('-s', '--supply', required=True,
                    help='vax supply schedule file')
    ap.add_argument('-o', '--output', required=True,
                    help='output file for optimal vax assignment schedule')
    ap.add_argument('--step_size', required=True,
                    help='grid search step size')
    ap.add_argument('--temporal_search_size', required=True,
                    help='grid search temporal granularity in days')
    ap.add_argument('--look_ahead',required=True,
                    help='look ahead duration for computing episize in days')
    ap.add_argument('-v', '--verbose',
                    help='increase log verbosity', action='store_true')
    ap.add_argument('-l', '--log', default=None,
                    help='log file, if not specified logs are written to'
                    ' standard output')
    ap.add_argument('-t', '--threat_model', default=False,
                    help='run the threat model optimization instead')

    return ap.parse_args()


def main():
    args = parse_args()

    level = logging.INFO
    if args.verbose:
        level = logging.DEBUG
    log.setLevel(level)
    if args.log is None:
        logging.basicConfig(format='%(asctime)s:%(levelname)s:'
                            '%(name)s.%(funcName)s:%(message)s',
                            datefmt='%Y%m%d-%H%M%S', level=level)
    else:
        logging.basicConfig(filename=args.log, level=level,
                            format='%(asctime)s:%(levelname)s:'
                            '%(name)s.%(funcName)s:%(message)s',
                            datefmt='%Y%m%d-%H%M%S')
    log.info('Starting optimization sim')

    cfg = sim.read_config(args.config)
    supply = opt.load_vax_supply(args.supply)
    if(args.threat_model is True):
        pass
    else:
        step_size = int(args.step_size)
        look_ahead = int(args.look_ahead)
        temporal_search_size = int(args.temporal_search_size)
        f_out = open(args.output,'w')
        opt.greedy(cfg,supply,step_size,temporal_search_size,look_ahead,f_out)
        f_out.close()


if __name__ == '__main__':
    main()
