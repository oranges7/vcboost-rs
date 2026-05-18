import os
import sys
import shutil
import argparse


def get_script_dir():
    return os.path.join(os.path.dirname(__file__), "scripts")


def main():
    parser = argparse.ArgumentParser(
        description="VCboost-RS: Rust-accelerated variant calling filter pipeline"
    )
    subparsers = parser.add_subparsers(dest="command", help="Available commands")

    subparsers.add_parser("scripts", help="Print the path to pipeline scripts")
    subparsers.add_parser("copy", help="Copy pipeline scripts to current directory")
    subparsers.add_parser("shell", help="Print the path to vcboost.sh")

    run_parser = subparsers.add_parser(
        "run", help="Run the full prediction pipeline"
    )
    run_parser.add_argument("-o", "--out_path", required=True, help="Output path")
    run_parser.add_argument("-b", "--bam_file", required=True, help="BAM file path")
    run_parser.add_argument("-v", "--vcf", required=True, help="VCF file path")
    run_parser.add_argument("-m", "--model", required=True, help="Model name")
    run_parser.add_argument("-r", "--ref_path", required=True, help="Reference file path")
    run_parser.add_argument("-d", "--model_path", required=True, help="Model directory path")
    run_parser.add_argument("-t", "--threads", type=int, default=32, help="Number of threads")
    run_parser.add_argument("-c", "--contig", default="chr1-22", help="Contig to process")
    run_parser.add_argument("-p", "--phase", action="store_true", help="Enable phase")

    args = parser.parse_args()

    if args.command is None:
        parser.print_help()
        sys.exit(1)

    script_dir = get_script_dir()

    if args.command == "scripts":
        print(script_dir)
    elif args.command == "shell":
        print(os.path.join(script_dir, "vcboost.sh"))
    elif args.command == "copy":
        dest = os.path.join(os.getcwd(), "vcboost-rs-scripts")
        if os.path.exists(dest):
            shutil.rmtree(dest)
        shutil.copytree(script_dir, dest)
        os.chmod(os.path.join(dest, "vcboost.sh"), 0o755)
        print(f"Scripts copied to: {dest}")
        print(f"Run: sh {dest}/vcboost.sh -o OUT -b BAM -v VCF -m MODEL -r REF")
    elif args.command == "run":
        shell_path = os.path.join(script_dir, "vcboost.sh")
        if not os.path.exists(shell_path):
            print(f"Error: vcboost.sh not found at {shell_path}", file=sys.stderr)
            sys.exit(1)
        cmd_parts = [shell_path]
        cmd_parts.extend(["-o", args.out_path])
        cmd_parts.extend(["-b", args.bam_file])
        cmd_parts.extend(["-v", args.vcf])
        cmd_parts.extend(["-m", args.model])
        cmd_parts.extend(["-r", args.ref_path])
        cmd_parts.extend(["-d", args.model_path])
        cmd_parts.extend(["-t", str(args.threads)])
        cmd_parts.extend(["-c", args.contig])
        if args.phase:
            cmd_parts.append("-p")
        os.execv("/bin/sh", ["sh", "-c", " ".join(cmd_parts)])


if __name__ == "__main__":
    main()
