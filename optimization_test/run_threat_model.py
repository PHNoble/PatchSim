#!/usr/bin/env python
import logging
import argparse
import threat_model as tm
import support_modules as sm

log = logging.getLogger("Threat_model")

def parse_args():
    ap = argparse.ArgumentParser(description='Find optimal vax assignment schedule with threat model.')
    ap.add_argument('-c', '--config', required=True,
                    help='simulation configuration file')
    ap.add_argument('-s', '--supply', required=True,
                    help='vax supply schedule file')
    ap.add_argument('-o', '--output', required=True,
                    help='output file for optimal vax assignment schedule')
    ap.add_argument('-v', '--verbose',
                    help='increase log verbosity', action='store_true')
    ap.add_argument('-l', '--log', default=None,
                    help='log file, if not specified logs are written to'
                    ' standard output')
    return ap.parse_args()

def main():

    args = parse_args()
    ##ONLY deterministic mode is available at the moment
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

    log.info('Starting threat model sim')
    cfg = sm.read_config(args.config)
    supply = sm.load_vax_supply(args.supply)

    f_out = open(args.output,'w')
    retepi = tm.run_disease_simulation(cfg,supply,log, return_epi=True, write_epi=True)
    log.info(retepi)
    f_out.close()


if __name__ == '__main__':
    main()
