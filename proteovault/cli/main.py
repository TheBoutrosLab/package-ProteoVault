"""Installed proteovault entry point."""

import argparse
import logging
import sys
import xml.etree.ElementTree as ET

import duckdb

from proteovault.cli import add_fasta, add_layer, add_uniprot_idmapping, build, inspect, update_paths, validate
from proteovault.version import __version__


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog="proteovault", description="Build and query portable proteomics results"
    )
    parser.add_argument("--version", action="version", version=__version__)
    commands = parser.add_subparsers(dest="command", required=True)
    for command in (build, add_layer, add_fasta, add_uniprot_idmapping, update_paths, validate, inspect):
        command.add_parser(commands)
    args = parser.parse_args(argv)
    logger = logging.getLogger("proteovault")
    previous_level = logger.level
    previous_propagate = logger.propagate
    logger.setLevel(logging.INFO)
    logger.propagate = False
    for handler in list(logger.handlers):
        if handler.name == "proteovault.cli":
            logger.removeHandler(handler)
            handler.close()
    handler = logging.StreamHandler(stream=sys.stderr)
    handler.name = "proteovault.cli"
    handler.setFormatter(
        logging.Formatter(
            fmt="%(asctime)s %(levelname)s %(name)s: %(message)s",
            datefmt="%Y-%m-%d %H:%M:%S",
        )
    )
    logger.addHandler(handler)
    try:
        return args.handler(args)
    except (OSError, ValueError, duckdb.Error, ET.ParseError) as error:
        logger.error("%s failed: %s", args.command, error)
        return 1
    finally:
        logger.removeHandler(handler)
        handler.close()
        logger.setLevel(previous_level)
        logger.propagate = previous_propagate
