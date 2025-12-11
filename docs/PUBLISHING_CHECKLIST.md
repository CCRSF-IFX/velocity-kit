# Publishing Checklist

## One-Time Setup ✓

- [ ] Create Test PyPI account (https://test.pypi.org)
- [ ] Create Production PyPI account (https://pypi.org)
- [ ] Generate Test PyPI API token
- [ ] Generate Production PyPI API token
- [ ] Add `TEST_PYPI_API_TOKEN` to GitHub Secrets
- [ ] Add `PYPI_API_TOKEN` to GitHub Secrets
- [ ] Push code to GitHub repository

## Before Each Release ✓

- [ ] Update version in `pyproject.toml`
- [ ] Update `CHANGELOG.md` (if you have one)
- [ ] Update `README.md` if needed
- [ ] Run tests locally
- [ ] Build package: `python -m build`
- [ ] Check package: `twine check dist/*`
- [ ] Clean dist folder: `rm -rf dist/`

## Test Release (Test PyPI) ✓

- [ ] Create test tag: `git tag v0.1.0-test`
- [ ] Push tag: `git push origin v0.1.0-test`
- [ ] Monitor GitHub Actions workflow
- [ ] Test installation: `pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ velocity-kit`
- [ ] Verify functionality: `velocity-kit --help`

## Production Release (PyPI) ✓

- [ ] Commit all changes
- [ ] Push to main branch
- [ ] Go to GitHub → Releases → Draft new release
- [ ] Create tag: `v0.1.0`
- [ ] Write release notes
- [ ] Publish release
- [ ] Monitor GitHub Actions workflow
- [ ] Verify on PyPI: https://pypi.org/project/velocity-kit/
- [ ] Test installation: `pip install velocity-kit`
- [ ] Verify functionality: `velocity-kit --version`

## Post-Release ✓

- [ ] Announce release (email, Slack, etc.)
- [ ] Update documentation if hosted separately
- [ ] Close related GitHub issues
- [ ] Plan next release
