# Testing Guide for ENSURE

This document provides information about testing the ENSURE systematic review automation pipeline.

## Overview

The ENSURE test suite is designed to validate all major components of the pipeline, ensuring reliability and reproducibility across different environments and use cases.

## Test Structure

### Test Categories

#### 1. **Unit Tests** (`pytest -m unit`)
- Fast, isolated tests for individual functions
- No external dependencies (database, API calls)
- Test specific algorithms and data processing logic

#### 2. **Integration Tests** (`pytest -m integration`)
- Test component interactions
- May use temporary databases or mock services
- Validate complete workflows

#### 3. **Mock Tests** (`pytest -m mock`)
- Test external service integrations (OpenAI API, PubMed)
- Use mocked responses to avoid API calls
- Validate error handling and edge cases

### Test Modules

```
tests/
├── test_pipeline.py           # Main test suite
├── test_performance.py        # Performance metrics tests (if separate)
├── test_database.py          # Database operation tests (if separate)
├── test_llm_integration.py   # LLM integration tests (if separate)
└── conftest.py               # Shared fixtures and configuration
```

## Running Tests

### Prerequisites

```bash
# Install test dependencies
pip install -r requirements-dev.txt

# Or install minimal testing requirements
pip install pytest pytest-cov pytest-mock
```

### Basic Usage

```bash
# Run all tests
pytest

# Run with verbose output
pytest -v

# Run specific test file
pytest tests/test_pipeline.py

# Run specific test class
pytest tests/test_pipeline.py::TestPerformanceMetrics

# Run specific test method
pytest tests/test_pipeline.py::TestPerformanceMetrics::test_perfect_classification
```

### Test Categories

```bash
# Run only fast unit tests
pytest -m unit

# Run integration tests
pytest -m integration

# Skip slow tests
pytest -m "not slow"

# Run tests with specific markers
pytest -m "unit and not slow"
```

### Coverage Reports

```bash
# Run tests with coverage
pytest --cov=3_pipeline

# Generate HTML coverage report
pytest --cov=3_pipeline --cov-report=html

# Generate coverage report for specific module
pytest --cov=3_pipeline.Performance --cov-report=term-missing
```

### Simple Test Runner

For users unfamiliar with pytest:

```bash
# Run all tests with simple interface
python run_tests.py

# Run only unit tests
python run_tests.py --unit

# Run with verbose output
python run_tests.py -v

# Install test dependencies
python run_tests.py --install-deps
```

## Test Classes and Coverage

### TestDatabaseOperations
**Coverage**: Database connectivity, CRUD operations, data integrity

**Key Tests**:
- `test_database_insert_and_query()`: Basic database operations
- `test_database_multiple_inserts()`: Bulk data handling
- `test_database_deletion()`: Data cleanup operations

**Fixtures**:
- `temp_db`: Temporary SQLite database for isolated testing

### TestPerformanceMetrics
**Coverage**: Performance calculation accuracy, edge cases

**Key Tests**:
- `test_perfect_classification()`: Validates metrics with perfect results
- `test_realistic_classification()`: Real-world performance scenarios
- `test_edge_cases()`: Boundary conditions and special cases

**Validates**:
- Sensitivity, Specificity, PPV, NPV calculations
- Confusion matrix accuracy
- Mathematical correctness

### TestFileOperations
**Coverage**: File I/O, data validation, error handling

**Key Tests**:
- `test_read_pmids_from_file()`: PMID file parsing
- `test_save_pmids_to_file()`: File writing operations
- `test_file_with_empty_lines()`: Robust file parsing

**Fixtures**:
- `temp_pmid_file`: Temporary PMID file with test data

### TestPromptGeneration
**Coverage**: LLM prompt creation, data formatting

**Key Tests**:
- `test_title_prompt_generation()`: Title screening prompts
- `test_abstract_prompt_generation()`: Abstract screening prompts
- `test_match_data_to_ids()`: PMID matching logic

**Validates**:
- Prompt structure and formatting
- Data inclusion and organization
- Error handling for invalid inputs

### TestScrambleFunction
**Coverage**: Data randomization, file manipulation

**Key Tests**:
- `test_scramble_preserves_content()`: Data integrity during randomization

**Validates**:
- All PMIDs preserved after scrambling
- File operations work correctly

### TestMockLLMIntegration
**Coverage**: LLM API integration, error handling

**Key Tests**:
- `test_screen_with_openai_success()`: Successful API interactions
- `test_screen_with_openai_api_error()`: Error handling

**Features**:
- Mocked OpenAI API responses
- No actual API calls during testing
- Validates response parsing

### TestBatchProcessing
**Coverage**: Data batching, memory management

**Key Tests**:
- `test_get_data_in_batches_titles()`: Batch processing logic

**Validates**:
- Correct batch sizes
- Data integrity across batches
- Database query efficiency

### TestDataValidation
**Coverage**: Input validation, error recovery

**Key Tests**:
- `test_invalid_pmid_format()`: Handling malformed data
- `test_empty_file_handling()`: Edge case handling
- `test_nonexistent_file_handling()`: Error conditions

### TestIntegrationScenarios
**Coverage**: End-to-end workflows, realistic use cases

**Key Tests**:
- `test_small_dataset_workflow()`: Complete pipeline simulation
- `test_prompt_workflow()`: Prompt generation and processing

## Writing New Tests

### Test Naming Convention

```python
class TestModuleName:
    def test_function_behavior_condition(self):
        """Test description explaining what is being validated."""
        pass
```

### Example Test Structure

```python
import pytest
from unittest.mock import Mock, patch

class TestNewFeature:
    """Test new feature functionality."""
    
    @pytest.fixture
    def sample_data(self):
        """Provide test data for multiple tests."""
        return {"pmid": "12345678", "title": "Test Title"}
    
    def test_normal_operation(self, sample_data):
        """Test normal operation with valid input."""
        # Arrange
        expected_result = "expected_output"
        
        # Act
        result = function_under_test(sample_data)
        
        # Assert
        assert result == expected_result
    
    def test_error_handling(self):
        """Test error handling with invalid input."""
        with pytest.raises(ValueError, match="Expected error message"):
            function_under_test(invalid_input)
    
    @patch('module.external_dependency')
    def test_with_mock(self, mock_dependency):
        """Test with mocked external dependency."""
        # Setup mock
        mock_dependency.return_value = "mocked_response"
        
        # Test
        result = function_under_test()
        
        # Verify
        assert result == "expected_result"
        mock_dependency.assert_called_once()
```

### Best Practices

1. **Use Descriptive Names**: Test names should clearly indicate what is being tested
2. **Follow AAA Pattern**: Arrange, Act, Assert
3. **Test One Thing**: Each test should validate one specific behavior
4. **Use Fixtures**: Share common test data and setup
5. **Mock External Dependencies**: Avoid real API calls, database operations
6. **Test Edge Cases**: Include boundary conditions and error scenarios
7. **Verify Side Effects**: Check that functions don't have unintended effects

### Fixtures for Common Operations

```python
@pytest.fixture
def temp_database():
    """Create temporary database for testing."""
    fd, temp_path = tempfile.mkstemp(suffix='.sqlite')
    os.close(fd)
    yield temp_path
    if os.path.exists(temp_path):
        os.unlink(temp_path)

@pytest.fixture
def sample_pmids():
    """Provide sample PMID data."""
    return ["12345678", "87654321", "11111111"]

@pytest.fixture
def mock_openai_response():
    """Mock OpenAI API response."""
    mock_response = Mock()
    mock_response.choices = [Mock()]
    mock_response.choices[0].message.content = "'12345678', '87654321'"
    return mock_response
```

## Continuous Integration

### GitHub Actions Example

```yaml
name: Tests

on: [push, pull_request]

jobs:
  test:
    runs-on: ubuntu-latest
    strategy:
      matrix:
        python-version: [3.8, 3.9, "3.10", "3.11"]
    
    steps:
    - uses: actions/checkout@v3
    - name: Set up Python ${{ matrix.python-version }}
      uses: actions/setup-python@v3
      with:
        python-version: ${{ matrix.python-version }}
    
    - name: Install dependencies
      run: |
        python -m pip install --upgrade pip
        pip install -r requirements.txt
        pip install -r requirements-dev.txt
    
    - name: Run tests
      run: |
        pytest --cov=3_pipeline --cov-report=xml
    
    - name: Upload coverage
      uses: codecov/codecov-action@v3
      with:
        file: ./coverage.xml
```

## Performance Testing

### Benchmarking

```python
import time
import pytest

@pytest.mark.slow
def test_performance_large_dataset():
    """Test performance with large dataset."""
    start_time = time.time()
    
    # Run operation with large dataset
    result = process_large_dataset(10000)
    
    end_time = time.time()
    processing_time = end_time - start_time
    
    # Assert reasonable performance
    assert processing_time < 30.0  # Should complete in under 30 seconds
    assert result is not None
```

### Memory Usage

```python
import psutil
import os

def test_memory_usage():
    """Test that memory usage stays within bounds."""
    process = psutil.Process(os.getpid())
    initial_memory = process.memory_info().rss
    
    # Run memory-intensive operation
    large_operation()
    
    final_memory = process.memory_info().rss
    memory_increase = final_memory - initial_memory
    
    # Assert memory usage is reasonable (e.g., less than 100MB increase)
    assert memory_increase < 100 * 1024 * 1024
```

## Troubleshooting Tests

### Common Issues

1. **Import Errors**: Ensure `3_pipeline` directory is in Python path
2. **Missing Dependencies**: Install with `pip install -r requirements-dev.txt`
3. **Database Errors**: Check that temporary databases are properly cleaned up
4. **API Mocking**: Verify that external services are properly mocked

### Debug Mode

```bash
# Run with detailed output
pytest -v -s

# Drop into debugger on failure
pytest --pdb

# Run specific failing test
pytest tests/test_pipeline.py::TestClassName::test_method_name -v -s
```

### Test Data Management

```python
# Ensure test isolation
def test_isolated_operation():
    # Create fresh test data
    test_data = create_test_data()
    
    # Run test
    result = operation(test_data)
    
    # Clean up
    cleanup_test_data(test_data)
    
    assert result == expected
```

## Contributing Test Improvements

When contributing new tests:

1. **Follow existing patterns**: Use similar structure to existing tests
2. **Add appropriate markers**: Mark tests as `unit`, `integration`, or `slow`
3. **Update documentation**: Add new test descriptions to this guide
4. **Ensure coverage**: Verify that new code is covered by tests
5. **Test on multiple platforms**: Ensure tests work on different operating systems

This comprehensive testing framework ensures the ENSURE pipeline maintains high quality and reliability across all environments and use cases.
