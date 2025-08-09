# Contributing to ENSURE

Thank you for your interest in contributing to the ENSURE systematic review automation pipeline! This document provides guidelines for contributing to the project.

## Code of Conduct

We are committed to fostering an open and welcoming environment. By participating in this project, you agree to abide by our code of conduct:

- Use welcoming and inclusive language
- Be respectful of differing viewpoints and experiences
- Gracefully accept constructive criticism
- Focus on what is best for the community
- Show empathy towards other community members

## Ways to Contribute

### 1. Bug Reports

If you find a bug, please create an issue with:
- Clear, descriptive title
- Detailed description of the problem
- Steps to reproduce the issue
- Expected vs. actual behavior
- Environment details (OS, Python version, package versions)
- Minimal code example if applicable

### 2. Feature Requests

For new features:
- Check existing issues to avoid duplicates
- Clearly describe the feature and its benefits
- Provide examples of how it would be used
- Consider backward compatibility

### 3. Documentation Improvements

Documentation contributions are always welcome:
- Fix typos or unclear explanations
- Add examples or use cases
- Improve API documentation
- Translate documentation

### 4. Code Contributions

#### Getting Started

1. **Fork the repository**
   ```bash
   git clone https://github.com/your-username/ensure.git
   cd ensure
   ```

2. **Set up development environment**
   ```bash
   python -m venv venv
   source venv/bin/activate  # On Windows: venv\Scripts\activate
   pip install -r requirements.txt
   pip install -r requirements-dev.txt  # If available
   ```

3. **Create a feature branch**
   ```bash
   git checkout -b feature/your-feature-name
   ```

#### Development Guidelines

**Code Style:**
- Follow PEP 8 for Python code
- Use meaningful variable and function names
- Add docstrings to all functions and classes
- Include type hints where appropriate

**Example:**
```python
def calculate_performance_metrics(
    initial: set, 
    gpt_selected: set, 
    goldstandard_selected: set
) -> tuple[float, float, float, float, int, int, int, int]:
    """
    Calculate comprehensive performance metrics for screening results.
    
    Args:
        initial: Set of all initial PMIDs
        gpt_selected: Set of PMIDs selected by LLM
        goldstandard_selected: Set of PMIDs in gold standard
        
    Returns:
        Tuple containing (sensitivity, specificity, PPV, NPV, tp, tn, fp, fn)
    """
    # Implementation here
    pass
```

**Testing:**
- Write tests for new functionality
- Ensure existing tests still pass
- Aim for high test coverage
- Test edge cases and error conditions

**Commit Messages:**
- Use clear, descriptive commit messages
- Start with a verb in present tense
- Keep first line under 50 characters
- Add detailed explanation if needed

Example:
```
Add support for custom batch sizes

- Allow users to specify batch sizes for title and abstract screening
- Update configuration documentation
- Add validation for batch size parameters
```

#### Submitting Changes

1. **Run tests**
   ```bash
   python test_pipeline.py
   python -m pytest tests/  # If using pytest
   ```

2. **Check code style**
   ```bash
   black .  # Format code
   flake8 .  # Check style
   mypy .   # Check types
   ```

3. **Update documentation**
   - Update relevant documentation files
   - Add docstrings to new functions
   - Update CHANGELOG.md if applicable

4. **Submit pull request**
   - Create pull request from your feature branch
   - Fill out the pull request template
   - Link to any relevant issues
   - Request review from maintainers

## Development Setup

### Additional Dependencies

For development work, you may need additional packages:

```bash
pip install black flake8 mypy pytest pytest-cov
```

### Environment Configuration

Create a `.env` file for development:
```
OPENAI_API_KEY=your_test_api_key
PUBMED_EMAIL=your_email@domain.com
DEBUG=True
```

### Running Tests

```bash
# Run basic validation
python test_pipeline.py

# Run unit tests (if using pytest)
pytest tests/

# Run with coverage
pytest --cov=. tests/

# Run specific test file
pytest tests/test_performance.py
```

## Project Structure for Contributors

```
ensure/
├── 3_pipeline/           # Core pipeline code
│   ├── ENSURE.py        # Main orchestrator
│   ├── PubMed.py        # PubMed interface
│   ├── GPT.py           # LLM interface  
│   ├── Performance.py   # Metrics calculation
│   └── dbconnect.py     # Database operations
├── tests/               # Test files
├── docs/                # Documentation
├── examples/            # Usage examples
└── scripts/             # Utility scripts
```

## Contribution Areas

### High Priority
- Performance optimization for large datasets
- Additional LLM provider support
- Improved error handling and logging
- Cross-platform compatibility testing

### Medium Priority
- Web interface for easier usage
- Additional statistical analysis methods
- Integration with reference managers
- Automated prompt optimization

### Documentation
- Video tutorials
- Use case examples
- API reference improvements
- Translation to other languages

## Review Process

1. **Automated Checks**: All PRs must pass automated tests and style checks
2. **Code Review**: At least one maintainer will review your code
3. **Testing**: Changes will be tested on different platforms
4. **Documentation Review**: Documentation changes will be verified
5. **Final Approval**: Maintainer approval required before merging

## Release Process

1. **Version Bumping**: Follow semantic versioning (MAJOR.MINOR.PATCH)
2. **Changelog**: Update CHANGELOG.md with new features and fixes
3. **Documentation**: Ensure all documentation is up to date
4. **Testing**: Run comprehensive tests across platforms
5. **Tagging**: Create git tag for the release

## Getting Help

### Communication Channels
- **GitHub Issues**: For bug reports and feature requests
- **GitHub Discussions**: For questions and general discussion
- **Email**: [contact_email] for private inquiries

### Resources
- [Python Style Guide (PEP 8)](https://pep8.org/)
- [Git Workflow Guide](https://guides.github.com/introduction/flow/)
- [Semantic Versioning](https://semver.org/)

## Recognition

Contributors will be recognized in:
- README.md contributor section
- Release notes
- Academic publications (for significant contributions)

## Legal

By contributing to ENSURE, you agree that your contributions will be licensed under the same MIT License that covers the project.

---

Thank you for contributing to ENSURE! Your efforts help advance automated systematic review research and benefit the entire scientific community.
