import os
import subprocess
import filecmp
import numpy as np
from pathlib import Path


tests = [
    {
        "name": "test_read",
        "input_file": "test_read.dat",
        "reference_dir": "ref_test_read",
        "files_to_compare": ["trans", "emisR", "mixingratios.dat"],
        "tolerance": 1e-6,
    },
    {
        "name": "test1",
        "input_file": "test1.dat",
        "reference_dir": "ref_test1",
        "files_to_compare": ["trans", "emisR"],
        "tolerance": 1e-6,
    },
    {
        "name": "test2",
        "input_file": "test2.dat",
        "reference_dir": "ref_test2",
        "files_to_compare": ["emisR", "mixingratios.dat"],
        "tolerance": 1e-6,
    },
    {
        "name": "Earth",
        "input_file": "Earth.in",
        "reference_dir": "ref_Earth",
        "files_to_compare": ["emisR", "cloudstructure0010.dat", "phase"],
        "tolerance": 1e-3,
    },
]



def run_ARCiS(input_file, output_dir):
    """
    Run the radiative transfer code with the given input files.
    Assumes the executable is a command-line tool.
    """
    cmd = ["ARCiS", input_file, "-o", output_dir]

    """Run the radiative transfer code and suppress all output."""
    with open(os.devnull, 'w') as devnull:
        ret = subprocess.run(
            cmd,
            stdout=devnull,
            stderr=devnull,
        )
    if ret.returncode != 0:
        raise RuntimeError(f"Command failed for {input_file} with return code {ret.returncode}")


def compare_files(output_file, reference_file, tolerance):
    """
    Compare only floating-point values in two files, column by column.
    Returns (result, maxerr).
    """
    with open(output_file, 'r') as f1, open(reference_file, 'r') as f2:
        lines1 = f1.readlines()
        lines2 = f2.readlines()

    if len(lines1) != len(lines2):
        return (False, 0.0)

    maxerr = 0.0
    for line1, line2 in zip(lines1, lines2):
        cols1 = line1.strip().split()
        cols2 = line2.strip().split()
        for val1, val2 in zip(cols1, cols2):
            try:
                num1 = float(val1)
                num2 = float(val2)
                if abs(num1 - num2) > maxerr * abs(num1 + num2) and abs(num1 - num2) > 1e-60:
                    maxerr = abs(num1 - num2) / abs(num1 + num2)
            except ValueError:
                continue
    if maxerr > tolerance:
        return (False, maxerr)
    else:
        return (True, maxerr)

def main():
	# ANSI color codes
	GREEN = "\033[92m"
	RED = "\033[91m"
	RESET = "\033[0m"

	all_passed = True
	for test in tests:
		# Run the radiative transfer code
		print("==========================================")
		print("Running: " + test["name"])
		output_dir="run_" + test["name"]
		run_ARCiS(test["input_file"], output_dir)

		# Compare output files
		for filename in test["files_to_compare"]:
			output_file = os.path.join(output_dir, filename)
			reference_file = os.path.join(test["reference_dir"], filename)

			if not os.path.exists(output_file):
				print(f"{RED}Error: {RESET}{output_file} not generated.")
				all_passed = False
				continue

			if not os.path.exists(reference_file):
				print(f"{RED}Error: {RESET}Reference file {reference_file} not found.")
				all_passed = False
				continue

			result, maxerr = compare_files(output_file, reference_file, test["tolerance"])
			if not result:
				print(f"{RED}FAIL: {RESET}{filename} does not match reference.")
				print(f"      maximum error = {RED}{maxerr:.2e}{RESET}")
				all_passed = False
			else:
				print(f"{GREEN}PASS: {RESET}{filename} matches reference.")
				if maxerr > 0:
					print(f"      maximum error = {GREEN}{maxerr:.2e}{RESET}")

	print("==========================================")
	if all_passed:
		print(f"{GREEN}All tests passed!{RESET}")
	else:
		print(f"{RED}Some tests failed.{RESET}")

if __name__ == "__main__":
    main()
