#!/usr/bin/env python3
"""
Simple test to check if the main script can start without errors.
"""

def test_import():
    print("=== Testing Import ===")
    try:
        import os
        import random
        import sqlite3
        import pandas as pd
        import json
        import re
        print("✓ Basic imports successful")

        # Test BioPython
        from Bio import Entrez
        print("✓ BioPython import successful")

        # Test AI libraries (these might fail if not installed)
        try:
            from groq import Groq
            print("✓ Groq import successful")
        except ImportError:
            print("⚠ Groq not available")

        try:
            from openai import OpenAI
            print("✓ OpenAI import successful")
        except ImportError:
            print("⚠ OpenAI not available")

        try:
            import google.generativeai as genai
            print("✓ Google AI import successful")
        except ImportError:
            print("⚠ Google AI not available")

        try:
            from ollama import Client as OllamaClient
            print("✓ Ollama import successful")
        except ImportError:
            print("⚠ Ollama not available")

        from flask import Flask, render_template, request, jsonify, send_file
        print("✓ Flask import successful")

    except Exception as e:
        print(f"✗ Import error: {str(e)}")
        import traceback
        traceback.print_exc()
        return False

    return True

def test_file_access():
    print("\n=== Testing File Access ===")
    files_to_check = ['Initial.txt', 'Goldstandard_Selected.txt', 'Prompts.xlsx']

    for filename in files_to_check:
        try:
            if os.path.exists(filename):
                size = os.path.getsize(filename)
                print(f"✓ {filename}: {size} bytes")
            else:
                print(f"✗ {filename}: File not found")
        except Exception as e:
            print(f"✗ {filename}: Error - {str(e)}")

    return True

def test_main_script():
    print("\n=== Testing Main Script Startup ===")

    try:
        # Try to import the main script
        import sys
        sys.path.append('.')

        # This will test if the script can be imported without syntax errors
        import screening_updated
        print("✓ Main script import successful")

        # Try to create Flask app
        app = screening_updated.app
        print("✓ Flask app creation successful")

    except Exception as e:
        print(f"✗ Main script error: {str(e)}")
        import traceback
        traceback.print_exc()
        return False

    return True

if __name__ == '__main__':
    print("Starting diagnostic tests...\n")

    success = True
    success &= test_import()
    success &= test_file_access()
    success &= test_main_script()

    print("=== Diagnostic Complete ===")
    if success:
        print("✓ All basic tests passed. The issue might be in the runtime execution.")
        print("Try running: python screening_updated.py")
    else:
        print("✗ Some tests failed. Check the errors above.")
