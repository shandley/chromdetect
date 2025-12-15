# Contributing to ChromDetect

Thank you for your interest in contributing to ChromDetect!

## How to Contribute

### Reporting Issues

- Check existing issues before opening a new one
- Include example data or scaffold names when reporting pattern recognition issues
- Provide your Python version and operating system

### Submitting Changes

1. Fork the repository
2. Create a feature branch: `git checkout -b feature/my-feature`
3. Make your changes
4. Add tests for new functionality
5. Run the test suite: `pytest`
6. Submit a pull request

### Development Setup

```bash
# Clone your fork
git clone https://github.com/YOUR_USERNAME/chromdetect.git
cd chromdetect

# Install in development mode with dev dependencies
pip install -e ".[dev]"

# Run tests
pytest

# Run tests with coverage
pytest --cov=chromdetect --cov-report=html

# Type checking
mypy chromdetect

# Linting
ruff check chromdetect
```

### Adding New Naming Patterns

If you encounter a naming convention that ChromDetect doesn't recognize:

1. **Open an issue** with example scaffold names
2. Or submit a PR adding the pattern:

```python
# In chromdetect/patterns.py

CHROMOSOME_PATTERNS.append(
    (r'^YourPattern_(\d+)$', 'your_pattern_name'),
)
```

Include tests in `tests/test_patterns.py`:

```python
@pytest.mark.parametrize("name,expected_class,expected_chr_id", [
    ("YourPattern_1", "chromosome", "1"),
    ("YourPattern_X", "chromosome", "X"),
])
def test_your_pattern(self, name, expected_class, expected_chr_id):
    classification, confidence, method, chr_id = detect_by_name(name)
    assert classification == expected_class
    assert chr_id == expected_chr_id
```

### Code Style

- Follow PEP 8 guidelines
- Use type hints for function signatures
- Write docstrings for public functions
- Keep line length under 100 characters
- Use ruff for linting

### Testing Guidelines

- Write tests for all new functionality
- Maintain test coverage above 80%
- Use pytest fixtures for shared test data
- Test edge cases (empty files, missing data, etc.)

## Releasing New Versions

For maintainers releasing new versions:

1. **Update version numbers** in:
   - `chromdetect/__init__.py`
   - `pyproject.toml`
   - `CITATION.cff`

2. **Update CHANGELOG.md** with release notes

3. **Create a GitHub release**:
   - Go to Releases > Draft a new release
   - Create a new tag (e.g., `v0.5.0`)
   - Use the CHANGELOG content for release notes
   - Publish the release

4. **Zenodo automatically archives** the release and mints a DOI
   - Update the DOI badge in README.md with the new DOI

5. **Publish to PyPI**:
   ```bash
   python -m build
   twine upload dist/*
   ```

## Zenodo Integration

ChromDetect uses [Zenodo](https://zenodo.org) for DOI minting and long-term archival:

- **CITATION.cff**: Standard citation metadata file (recognized by GitHub)
- **.zenodo.json**: Additional Zenodo-specific metadata

When a GitHub release is created, Zenodo automatically:
1. Archives the release
2. Mints a DOI for that specific version
3. Updates the concept DOI (points to latest version)

## Questions?

Open an issue or start a discussion on GitHub.
