"""Command-line interface for the variant pathogenicity pipeline."""

import argparse
import sys

from .pipeline import run_pipeline


def main():
    parser = argparse.ArgumentParser(
        description="Variant Pathogenicity Pipeline — structural analysis of missense variants",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""\
Example:
  varpat config/my_project.yaml
  python -m varpat config/my_project.yaml

See config/example_config.yaml for configuration reference.
""",
    )
    parser.add_argument('config', help='Path to YAML configuration file')
    parser.add_argument('--version', action='version', version='varpat 1.0.0')

    args = parser.parse_args()

    try:
        run_pipeline(args.config)
    except FileNotFoundError as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Pipeline error: {e}", file=sys.stderr)
        raise


if __name__ == '__main__':
    main()
