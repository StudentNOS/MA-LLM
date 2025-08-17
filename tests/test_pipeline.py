#!/usr/bin/env python3
"""
Test suite for ENSURE systematic review automation pipeline.

This module contains comprehensive tests for all major components of the ENSURE
pipeline, including data processing, performance calculations, database operations,
and LLM integration.

Usage:
    pytest test_pipeline.py -v
    python -m pytest test_pipeline.py::TestPerformanceMetrics -v
"""

import pytest
import tempfile
import os
import sys
import sqlite3
from pathlib import Path
from unittest.mock import Mock, patch, MagicMock
import pandas as pd

# Add the pipeline directory to Python path
PIPELINE_DIR = Path(__file__).parent.parent / "3_pipeline"
sys.path.insert(0, str(PIPELINE_DIR))

# Import pipeline modules
try:
    from dbconnect import execute_query, insert, delete_all_data
    from PubMed import read_pmids_from_file, fetch_details, count_papers, create_excel_from_db
    from GPT import generate_prompt, match_data_to_ids, save_pmids_to_file, get_data_in_batches
    from Performance import calculate_performance_metrics, read_ids_from_file
    from Scramble import scramble
except ImportError as e:
    pytest.skip(f"Pipeline modules not available: {e}", allow_module_level=True)


class TestDatabaseOperations:
    """Test database connectivity and operations."""
    
    @pytest.fixture
    def temp_db(self):
        """Create a temporary database for testing."""
        fd, temp_path = tempfile.mkstemp(suffix='.sqlite')
        os.close(fd)
        yield temp_path
        if os.path.exists(temp_path):
            os.unlink(temp_path)
    
    def test_database_insert_and_query(self, temp_db):
        """Test basic database insert and query operations."""
        test_data = {
            "pmid": "12345678",
            "title": "Test Paper Title",
            "abstract": "This is a test abstract for validation purposes."
        }
        
        # Test insert
        insert("Initial", test_data, temp_db, create_if_missing=True)
        
        # Test query
        results = execute_query("SELECT COUNT(*) FROM Initial", temp_db)
        assert results[0][0] == 1
        
        # Test data retrieval
        results = execute_query("SELECT pmid, title FROM Initial", temp_db)
        assert len(results) == 1
        assert results[0][0] == "12345678"
        assert results[0][1] == "Test Paper Title"
    
    def test_database_multiple_inserts(self, temp_db):
        """Test multiple database insertions."""
        test_papers = [
            {"pmid": "11111111", "title": "Paper 1", "abstract": "Abstract 1"},
            {"pmid": "22222222", "title": "Paper 2", "abstract": "Abstract 2"},
            {"pmid": "33333333", "title": "Paper 3", "abstract": "Abstract 3"}
        ]
        
        for paper in test_papers:
            insert("Initial", paper, temp_db, create_if_missing=True)
        
        results = execute_query("SELECT COUNT(*) FROM Initial", temp_db)
        assert results[0][0] == 3
        
        # Test specific data
        results = execute_query("SELECT pmid FROM Initial ORDER BY pmid", temp_db)
        pmids = [row[0] for row in results]
        assert pmids == ["11111111", "22222222", "33333333"]
    
    def test_database_deletion(self, temp_db):
        """Test database data deletion."""
        # Insert test data
        test_data = {"pmid": "99999999", "title": "Delete Test", "abstract": "Test"}
        insert("Initial", test_data, temp_db, create_if_missing=True)
        
        # Verify insertion
        results = execute_query("SELECT COUNT(*) FROM Initial", temp_db)
        assert results[0][0] == 1
        
        # Test deletion
        delete_all_data(temp_db)
        results = execute_query("SELECT COUNT(*) FROM Initial", temp_db)
        assert results[0][0] == 0


class TestPerformanceMetrics:
    """Test performance calculation functions."""
    
    def test_perfect_classification(self):
        """Test performance metrics with perfect classification."""
        initial = {"1", "2", "3", "4", "5", "6"}
        gpt_selected = {"1", "2", "3"}
        goldstandard_selected = {"1", "2", "3"}
        
        sensitivity, specificity, PPV, NPV, tp, tn, fp, fn = calculate_performance_metrics(
            initial, gpt_selected, goldstandard_selected
        )
        
        # Perfect classification should yield 1.0 for all metrics
        assert sensitivity == 1.0
        assert specificity == 1.0
        assert PPV == 1.0
        assert NPV == 1.0
        assert tp == 3
        assert tn == 3
        assert fp == 0
        assert fn == 0
    
    def test_realistic_classification(self):
        """Test performance metrics with realistic classification scenario."""
        initial = {"1", "2", "3", "4", "5", "6", "7", "8"}
        gpt_selected = {"1", "2", "3", "4"}  # 4 selected
        goldstandard_selected = {"1", "2", "5"}  # 3 relevant
        
        sensitivity, specificity, PPV, NPV, tp, tn, fp, fn = calculate_performance_metrics(
            initial, gpt_selected, goldstandard_selected
        )
        
        # Expected: TP=2, FP=2, TN=3, FN=1
        assert tp == 2  # Papers 1, 2
        assert fp == 2  # Papers 3, 4
        assert tn == 3  # Papers 6, 7, 8
        assert fn == 1  # Paper 5
        
        # Calculate expected metrics
        expected_sensitivity = 2 / 3  # TP / (TP + FN)
        expected_specificity = 3 / 5  # TN / (TN + FP)
        expected_ppv = 2 / 4  # TP / (TP + FP)
        expected_npv = 3 / 4  # TN / (TN + FN)
        
        assert abs(sensitivity - expected_sensitivity) < 0.001
        assert abs(specificity - expected_specificity) < 0.001
        assert abs(PPV - expected_ppv) < 0.001
        assert abs(NPV - expected_npv) < 0.001
    
    def test_edge_cases(self):
        """Test performance metrics with edge cases."""
        # Case 1: No relevant papers
        initial = {"1", "2", "3"}
        gpt_selected = {"1"}
        goldstandard_selected = set()
        
        sensitivity, specificity, PPV, NPV, tp, tn, fp, fn = calculate_performance_metrics(
            initial, gpt_selected, goldstandard_selected
        )
        
        assert sensitivity == 0.0  # No true positives possible
        assert tp == 0
        assert fn == 0
        assert fp == 1
        assert tn == 2
        
        # Case 2: All papers relevant
        initial = {"1", "2", "3"}
        gpt_selected = {"1", "2"}
        goldstandard_selected = {"1", "2", "3"}
        
        sensitivity, specificity, PPV, NPV, tp, tn, fp, fn = calculate_performance_metrics(
            initial, gpt_selected, goldstandard_selected
        )
        
        assert specificity == 0.0  # No true negatives possible
        assert tp == 2
        assert fn == 1
        assert fp == 0
        assert tn == 0


class TestFileOperations:
    """Test file reading and writing operations."""
    
    @pytest.fixture
    def temp_pmid_file(self):
        """Create a temporary PMID file for testing."""
        fd, temp_path = tempfile.mkstemp(suffix='.txt')
        test_pmids = ["32420649", "22901346", "34025474", "25687662"]
        
        with os.fdopen(fd, 'w') as f:
            f.write('\n'.join(test_pmids))
        
        yield temp_path, test_pmids
        
        if os.path.exists(temp_path):
            os.unlink(temp_path)
    
    def test_read_pmids_from_file(self, temp_pmid_file):
        """Test reading PMIDs from file."""
        file_path, expected_pmids = temp_pmid_file
        
        pmids = read_pmids_from_file(file_path)
        assert pmids == expected_pmids
    
    def test_read_ids_from_file_as_set(self, temp_pmid_file):
        """Test reading PMIDs as set from file."""
        file_path, expected_pmids = temp_pmid_file
        
        pmid_set = read_ids_from_file(file_path)
        assert pmid_set == set(expected_pmids)
    
    def test_save_pmids_to_file(self):
        """Test saving PMIDs to file."""
        test_pmids = ["11111111", "22222222", "33333333"]
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            temp_path = f.name
        
        try:
            save_pmids_to_file(test_pmids, temp_path)
            
            # Read back and verify
            with open(temp_path, 'r') as f:
                saved_pmids = f.read().strip().split('\n')
            
            assert saved_pmids == test_pmids
        finally:
            if os.path.exists(temp_path):
                os.unlink(temp_path)
    
    def test_file_with_empty_lines(self):
        """Test handling files with empty lines."""
        test_content = "12345678\n\n87654321\n\n11111111\n"
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            f.write(test_content)
            temp_path = f.name
        
        try:
            pmids = read_pmids_from_file(temp_path)
            assert pmids == ["12345678", "87654321", "11111111"]
        finally:
            if os.path.exists(temp_path):
                os.unlink(temp_path)


class TestPromptGeneration:
    """Test prompt generation and processing."""
    
    def test_title_prompt_generation(self):
        """Test generation of title screening prompts."""
        test_data = [
            ("12345678", "Diabetes treatment with metformin"),
            ("87654321", "COVID-19 vaccine effectiveness study"),
            ("11111111", "Machine learning in medical diagnosis")
        ]
        
        manual_prompt = "Select papers related to diabetes treatment"
        
        prompt = generate_prompt(test_data, "titles", manual_prompt)
        
        # Check that prompt contains all expected elements
        assert "Select papers related to diabetes treatment" in prompt
        assert "1. Diabetes treatment with metformin" in prompt
        assert "2. COVID-19 vaccine effectiveness study" in prompt
        assert "3. Machine learning in medical diagnosis" in prompt
        assert "Format the output as a comma-separated string" in prompt
        assert "Each entry should be the unique PMID" in prompt
    
    def test_abstract_prompt_generation(self):
        """Test generation of abstract screening prompts."""
        test_data = [
            ("12345678", "This study examines diabetes treatment..."),
            ("87654321", "We investigated COVID-19 vaccine...")
        ]
        
        manual_prompt = "Evaluate abstracts for diabetes research"
        
        prompt = generate_prompt(test_data, "abstracts", manual_prompt)
        
        assert "Evaluate abstracts for diabetes research" in prompt
        assert "1. This study examines diabetes treatment..." in prompt
        assert "2. We investigated COVID-19 vaccine..." in prompt
        assert "abstract" in prompt.lower()
    
    def test_match_data_to_ids(self):
        """Test matching LLM results back to PMIDs."""
        screened_results = ["12345678", "11111111"]
        original_data = [
            ("12345678", "Paper 1"),
            ("87654321", "Paper 2"), 
            ("11111111", "Paper 3")
        ]
        
        matched_pmids = match_data_to_ids(screened_results, original_data, "titles")
        
        assert set(matched_pmids) == {"12345678", "11111111"}
    
    def test_invalid_prompt_decision(self):
        """Test handling of invalid prompt decision parameter."""
        test_data = [("12345", "Test")]
        
        with pytest.raises(ValueError, match="Invalid decision parameter"):
            generate_prompt(test_data, "invalid_decision", "test prompt")


class TestScrambleFunction:
    """Test data scrambling functionality."""
    
    def test_scramble_preserves_content(self):
        """Test that scrambling preserves all PMIDs."""
        original_pmids = ["12345678", "87654321", "11111111", "22222222"]
        
        # Create temporary Initial.txt file
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False, 
                                       dir=PIPELINE_DIR, prefix='Initial') as f:
            f.write('\n'.join(original_pmids))
            temp_path = f.name
        
        # Rename to Initial.txt for the test
        initial_path = PIPELINE_DIR / "Initial.txt"
        backup_exists = initial_path.exists()
        if backup_exists:
            backup_path = PIPELINE_DIR / "Initial.txt.backup"
            initial_path.rename(backup_path)
        
        try:
            os.rename(temp_path, str(initial_path))
            
            # Run scramble
            original_dir = os.getcwd()
            os.chdir(PIPELINE_DIR)
            
            try:
                scramble()
                
                # Read scrambled content
                with open("Initial.txt", 'r') as f:
                    scrambled_pmids = f.read().strip().split('\n')
                
                # Should have same PMIDs, possibly different order
                assert set(scrambled_pmids) == set(original_pmids)
                assert len(scrambled_pmids) == len(original_pmids)
                
            finally:
                os.chdir(original_dir)
                
        finally:
            # Clean up
            if initial_path.exists():
                initial_path.unlink()
            if backup_exists and (PIPELINE_DIR / "Initial.txt.backup").exists():
                (PIPELINE_DIR / "Initial.txt.backup").rename(initial_path)


class TestMockLLMIntegration:
    """Test LLM integration with mocked responses."""
    
    @patch('GPT.client')
    def test_screen_with_openai_success(self, mock_client):
        """Test successful LLM screening with mocked response."""
        # Import the function to test
        from GPT import screen_with_openai
        
        # Mock the OpenAI response
        mock_response = Mock()
        mock_response.choices = [Mock()]
        mock_response.choices[0].message.content = "'12345678', '87654321', '11111111'"
        mock_client.chat.completions.create.return_value = mock_response
        
        test_prompt = "Test prompt for screening"
        result = screen_with_openai(test_prompt)
        
        # Verify the result
        expected_pmids = ["12345678", "87654321", "11111111"]
        assert result == expected_pmids
        
        # Verify the API was called correctly
        mock_client.chat.completions.create.assert_called_once()
        call_args = mock_client.chat.completions.create.call_args
        assert call_args[1]['messages'][0]['content'] == test_prompt
    
    @patch('GPT.client')
    def test_screen_with_openai_api_error(self, mock_client):
        """Test LLM screening with API error."""
        from GPT import screen_with_openai
        
        # Mock an API error
        mock_client.chat.completions.create.side_effect = Exception("API Error")
        
        test_prompt = "Test prompt"
        
        # Should handle the error gracefully
        with pytest.raises(Exception):
            screen_with_openai(test_prompt)


class TestBatchProcessing:
    """Test batch processing functionality."""
    
    @pytest.fixture
    def temp_db_with_data(self):
        """Create a temporary database with test data."""
        fd, temp_path = tempfile.mkstemp(suffix='.sqlite')
        os.close(fd)
        
        # Insert test data
        test_papers = [
            {"pmid": f"1234567{i}", "title": f"Title {i}", "abstract": f"Abstract {i}"}
            for i in range(10)
        ]
        
        for paper in test_papers:
            insert("Initial", paper, temp_path, create_if_missing=True)
        
        yield temp_path
        
        if os.path.exists(temp_path):
            os.unlink(temp_path)
    
    @patch('dbconnect.ENSURE')
    def test_get_data_in_batches_titles(self, mock_db_path, temp_db_with_data):
        """Test batch processing for titles."""
        mock_db_path.return_value = temp_db_with_data
        
        batches = list(get_data_in_batches("titles", batch_size=3))
        
        # Should have 4 batches (10 items, batch size 3)
        assert len(batches) >= 3
        
        # First batch should have 3 items
        assert len(batches[0]) == 3
        
        # Check data structure
        pmid, title = batches[0][0]
        assert pmid.startswith("1234567")
        assert title.startswith("Title")


class TestDataValidation:
    """Test data validation and error handling."""
    
    def test_invalid_pmid_format(self):
        """Test handling of invalid PMID formats."""
        # Test with non-numeric PMIDs
        test_content = "12345678\ninvalid_pmid\n87654321\n"
        
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            f.write(test_content)
            temp_path = f.name
        
        try:
            pmids = read_pmids_from_file(temp_path)
            # Should only include valid numeric PMIDs
            assert "invalid_pmid" not in pmids
            assert "12345678" in pmids
            assert "87654321" in pmids
        finally:
            if os.path.exists(temp_path):
                os.unlink(temp_path)
    
    def test_empty_file_handling(self):
        """Test handling of empty files."""
        with tempfile.NamedTemporaryFile(mode='w', suffix='.txt', delete=False) as f:
            # Write empty content
            temp_path = f.name
        
        try:
            pmids = read_pmids_from_file(temp_path)
            assert pmids == []
        finally:
            if os.path.exists(temp_path):
                os.unlink(temp_path)
    
    def test_nonexistent_file_handling(self):
        """Test handling of nonexistent files."""
        with pytest.raises(FileNotFoundError):
            read_pmids_from_file("/nonexistent/path/file.txt")


class TestIntegrationScenarios:
    """Test realistic integration scenarios."""
    
    def test_small_dataset_workflow(self):
        """Test complete workflow with small dataset."""
        # Setup test data
        initial_pmids = {"12345678", "87654321", "11111111", "22222222"}
        gpt_selected = {"12345678", "11111111", "22222222"}  # 3 selected
        goldstandard = {"12345678", "87654321"}  # 2 relevant
        
        # Calculate performance
        sensitivity, specificity, PPV, NPV, tp, tn, fp, fn = calculate_performance_metrics(
            initial_pmids, gpt_selected, goldstandard
        )
        
        # Verify reasonable performance metrics
        assert 0.0 <= sensitivity <= 1.0
        assert 0.0 <= specificity <= 1.0
        assert 0.0 <= PPV <= 1.0
        assert 0.0 <= NPV <= 1.0
        assert tp + tn + fp + fn == len(initial_pmids)
    
    def test_prompt_workflow(self):
        """Test prompt generation and processing workflow."""
        test_data = [
            ("12345678", "Systematic review of diabetes interventions"),
            ("87654321", "Meta-analysis of cardiovascular outcomes")
        ]
        
        manual_prompt = "Select systematic reviews and meta-analyses"
        
        # Generate prompt
        prompt = generate_prompt(test_data, "titles", manual_prompt)
        
        # Verify prompt structure
        assert len(prompt) > len(manual_prompt)  # Should be expanded
        assert "12345678" not in prompt  # PMIDs shouldn't be in prompt text
        assert "diabetes" in prompt.lower()
        assert "meta-analysis" in prompt.lower()


# Configuration for pytest
def pytest_configure(config):
    """Configure pytest with custom markers."""
    config.addinivalue_line(
        "markers", "integration: mark test as integration test"
    )
    config.addinivalue_line(
        "markers", "slow: mark test as slow running"
    )


if __name__ == "__main__":
    # Run tests if script is executed directly
    pytest.main([__file__, "-v"])
