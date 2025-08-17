#!/usr/bin/env python3
"""
Simple test runner for ENSURE pipeline tests.

This script provides an easy way to run tests without requiring pytest knowledge.
It can run all tests or specific test categories.

Usage:
    python run_tests.py                    # Run all tests
    python run_tests.py --unit             # Run only unit tests
    python run_tests.py --integration      # Run only integration tests
    python run_tests.py --help             # Show help
"""

import sys
import subprocess
import argparse
from pathlib import Path

def check_pytest_available():
    """Check if pytest is available."""
    try:
        import pytest
        return True
    except ImportError:
        return False

def install_pytest():
    """Install pytest if user agrees."""
    response = input("pytest is required but not installed. Install it? (y/n): ")
    if response.lower().startswith('y'):
        try:
            subprocess.check_call([sys.executable, "-m", "pip", "install", "pytest"])
            return True
        except subprocess.CalledProcessError:
            print("Failed to install pytest. Please install manually: pip install pytest")
            return False
    return False

def run_tests(test_type=None, verbose=False):
    """Run tests with specified options."""
    cmd = [sys.executable, "-m", "pytest"]
    
    if verbose:
        cmd.append("-v")
    
    # Add test directory
    test_dir = Path(__file__).parent / "tests"
    cmd.append(str(test_dir))
    
    # Add markers for specific test types
    if test_type == "unit":
        cmd.extend(["-m", "unit"])
    elif test_type == "integration":
        cmd.extend(["-m", "integration"])
    elif test_type == "fast":
        cmd.extend(["-m", "not slow"])
    
    print(f"Running: {' '.join(cmd)}")
    print("-" * 50)
    
    try:
        result = subprocess.run(cmd, capture_output=False)
        return result.returncode
    except FileNotFoundError:
        print("Error: Could not run pytest. Make sure it's installed.")
        return 1

def main():
    """Main test runner function."""
    parser = argparse.ArgumentParser(
        description="Run ENSURE pipeline tests",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
    python run_tests.py                 # Run all tests
    python run_tests.py --unit          # Run unit tests only
    python run_tests.py --integration   # Run integration tests only
    python run_tests.py --fast          # Skip slow tests
    python run_tests.py -v              # Verbose output
        """
    )
    
    parser.add_argument(
        "--unit", 
        action="store_true", 
        help="Run only unit tests (fast, isolated tests)"
    )
    parser.add_argument(
        "--integration", 
        action="store_true", 
        help="Run only integration tests (may be slower)"
    )
    parser.add_argument(
        "--fast", 
        action="store_true", 
        help="Skip slow tests"
    )
    parser.add_argument(
        "-v", "--verbose", 
        action="store_true", 
        help="Verbose output"
    )
    parser.add_argument(
        "--install-deps", 
        action="store_true", 
        help="Install missing test dependencies"
    )
    
    args = parser.parse_args()
    
    # Install dependencies if requested
    if args.install_deps:
        print("Installing test dependencies...")
        try:
            subprocess.check_call([
                sys.executable, "-m", "pip", "install", "-r", "requirements-dev.txt"
            ])
            print("Test dependencies installed successfully.")
        except subprocess.CalledProcessError:
            print("Failed to install dependencies. Please install manually.")
            return 1
    
    # Check pytest availability
    if not check_pytest_available():
        if not install_pytest():
            return 1
    
    # Determine test type
    test_type = None
    if args.unit:
        test_type = "unit"
    elif args.integration:
        test_type = "integration"
    elif args.fast:
        test_type = "fast"
    
    # Run tests
    return_code = run_tests(test_type, args.verbose)
    
    if return_code == 0:
        print("\n✅ All tests passed!")
    else:
        print(f"\n❌ Tests failed with return code {return_code}")
    
    return return_code

if __name__ == "__main__":
    sys.exit(main())
