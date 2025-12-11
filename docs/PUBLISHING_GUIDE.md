# Step-by-Step Guide: Publishing velocity-kit to PyPI

This guide walks you through setting up automated publishing to Test PyPI and PyPI using GitHub Actions.

## Prerequisites

- [x] GitHub repository created
- [ ] PyPI account (https://pypi.org/account/register/)
- [ ] Test PyPI account (https://test.pypi.org/account/register/)
- [ ] API tokens generated (see below)

---

## Part 1: Initial Setup (One-Time)

### Step 1: Create PyPI Accounts

1. **Create Test PyPI account**:
   - Go to https://test.pypi.org/account/register/
   - Create account with your email
   - Verify email address

2. **Create Production PyPI account**:
   - Go to https://pypi.org/account/register/
   - Create account with your email
   - Verify email address

### Step 2: Generate API Tokens

#### For Test PyPI:

1. Log in to https://test.pypi.org
2. Go to Account Settings → API tokens
3. Click "Add API token"
4. Token name: `velocity-kit-github-actions`
5. Scope: Select "Entire account" (initially, you can narrow this later)
6. Click "Add token"
7. **IMPORTANT**: Copy the token immediately (starts with `pypi-`)
8. Save it securely - you won't see it again!

#### For Production PyPI:

1. Log in to https://pypi.org
2. Go to Account Settings → API tokens
3. Click "Add API token"
4. Token name: `velocity-kit-github-actions`
5. Scope: Select "Entire account" (initially)
6. Click "Add token"
7. **IMPORTANT**: Copy the token immediately (starts with `pypi-`)
8. Save it securely!

### Step 3: Add Secrets to GitHub Repository

1. Go to your GitHub repository
2. Click **Settings** (top right)
3. In left sidebar, click **Secrets and variables** → **Actions**
4. Click **New repository secret**

#### Add Test PyPI Token:
- Name: `TEST_PYPI_API_TOKEN`
- Secret: Paste the Test PyPI token (starts with `pypi-`)
- Click "Add secret"

#### Add Production PyPI Token:
- Click **New repository secret** again
- Name: `PYPI_API_TOKEN`
- Secret: Paste the Production PyPI token (starts with `pypi-`)
- Click "Add secret"

### Step 4: Initialize Git Repository (if not done)

```bash
cd /mnt/ccrsf-static/singlecell_projects/archive_clean/SusanMackem_CS038218_8PIPseq_021025/trajectory_analysis/scripts

# Initialize git if needed
git init

# Add all files
git add .

# Create initial commit
git commit -m "Initial commit: velocity-kit v0.1.0"

# Add remote (replace with your actual repo URL)
git remote add origin https://github.com/CCRSF-IFX/velocity-kit.git

# Push to GitHub
git branch -M main
git push -u origin main
```

---

## Part 2: Publishing to Test PyPI

### Step 5: Test Locally First

Before publishing, test the build process locally:

```bash
cd /mnt/ccrsf-static/singlecell_projects/archive_clean/SusanMackem_CS038218_8PIPseq_021025/trajectory_analysis/scripts

# Install build tools
pip install build twine

# Build the package
python -m build

# Check the build
twine check dist/*

# You should see:
# Checking dist/velocity_kit-0.1.0-py3-none-any.whl: PASSED
# Checking dist/velocity-kit-0.1.0.tar.gz: PASSED
```

### Step 6: Manual Test Upload to Test PyPI (Optional)

To test manually before using GitHub Actions:

```bash
# Upload to Test PyPI
twine upload --repository testpypi dist/*

# You'll be prompted for username and password:
# username: __token__
# password: [paste your Test PyPI token]
```

### Step 7: Test Installation from Test PyPI

```bash
# In a new environment
pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ velocity-kit

# Test it works
velocity-kit --help
```

### Step 8: Automated Test PyPI Publishing via GitHub Actions

To publish to Test PyPI automatically:

```bash
# Create a test tag (triggers the workflow)
git tag v0.1.0-test
git push origin v0.1.0-test
```

This will:
1. Trigger the `publish-test-pypi.yml` workflow
2. Build the package
3. Publish to Test PyPI automatically

Monitor progress:
- Go to your GitHub repo → **Actions** tab
- Watch the "Publish to Test PyPI" workflow run

---

## Part 3: Publishing to Production PyPI

### Step 9: Create a GitHub Release

When you're ready to publish to production PyPI:

1. Go to your GitHub repository
2. Click **Releases** (right sidebar)
3. Click **Draft a new release**
4. Fill in the form:
   - **Choose a tag**: Create new tag `v0.1.0` (on publish)
   - **Release title**: `velocity-kit v0.1.0`
   - **Description**: 
     ```markdown
     ## velocity-kit v0.1.0 - Initial Release
     
     A cross-platform toolkit for building RNA velocity-ready spliced/unspliced matrices.
     
     ### Features
     - ✅ PIPseq/PIPseeker support
     - ✅ 10x Genomics CellRanger support
     - ✅ Dual-run subtraction method
     - ✅ Flexible output formats (h5ad and/or loom)
     - ✅ Optional scVelo preprocessing
     - ✅ Python 3.7-3.14 support
     
     ### Installation
     ```bash
     pip install velocity-kit
     ```
     
     ### Documentation
     See [README.md](https://github.com/CCRSF-IFX/velocity-kit#readme) for usage instructions.
     ```
5. Check **Set as the latest release**
6. Click **Publish release**

This will:
1. Trigger the `publish-pypi.yml` workflow
2. Build the package
3. Publish to Production PyPI automatically
4. Make it available via `pip install velocity-kit`

### Step 10: Verify Production Installation

```bash
# In a new environment
pip install velocity-kit

# Test it works
velocity-kit --help
velocity-kit --version
```

---

## Part 4: Ongoing Workflow

### For Future Releases

#### Testing New Features:

```bash
# 1. Make your changes
git add .
git commit -m "Add new feature"
git push

# 2. Test on Test PyPI
git tag v0.2.0-test
git push origin v0.2.0-test

# 3. Verify on Test PyPI
pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ velocity-kit
```

#### Production Release:

```bash
# 1. Update version in pyproject.toml
# Change: version = "0.2.0"

# 2. Commit version bump
git add pyproject.toml
git commit -m "Bump version to 0.2.0"
git push

# 3. Create GitHub Release (as in Step 9)
# - Tag: v0.2.0
# - Title: velocity-kit v0.2.0
# - Description: List changes
# - Publish
```

---

## Troubleshooting

### Common Issues:

**1. "Invalid or non-existent authentication information"**
- Solution: Check that secrets are correctly named in GitHub
  - `TEST_PYPI_API_TOKEN` for Test PyPI
  - `PYPI_API_TOKEN` for Production PyPI

**2. "File already exists"**
- Solution: You can't re-upload the same version
  - Bump version number in `pyproject.toml`
  - Create new tag

**3. "Package name already taken"**
- Solution: The package name `velocity-kit` must be unique on PyPI
  - Check availability at https://pypi.org/project/velocity-kit/
  - If taken, choose a different name in `pyproject.toml`

**4. Workflow doesn't trigger**
- Solution: Check trigger conditions in `.github/workflows/*.yml`
  - Test PyPI: triggered by tags matching `v*-test`
  - Production PyPI: triggered by GitHub releases

### Check Workflow Status

1. Go to your repo → **Actions** tab
2. Click on the workflow run
3. Click on the job to see detailed logs
4. Look for errors in red

---

## Security Best Practices

1. **Never commit tokens**: API tokens are in GitHub Secrets only
2. **Use scoped tokens**: After first upload, create package-specific tokens
3. **Rotate tokens**: Change tokens periodically
4. **Review releases**: Always review the Release notes before publishing
5. **Test first**: Always test on Test PyPI before production

---

## Quick Reference Commands

```bash
# Local build and test
python -m build
twine check dist/*

# Manual Test PyPI upload
twine upload --repository testpypi dist/*

# Manual Production PyPI upload
twine upload dist/*

# Create test release (triggers Test PyPI)
git tag v0.1.0-test && git push origin v0.1.0-test

# Create production release (use GitHub UI for releases)
# Don't use: git tag v0.1.0  # This won't trigger the workflow
# Instead: Use GitHub Releases UI
```

---

## Checklist Before First Release

- [ ] Repository on GitHub with all code pushed
- [ ] PyPI account created and verified
- [ ] Test PyPI account created and verified
- [ ] API tokens generated for both PyPI and Test PyPI
- [ ] GitHub secrets added (`TEST_PYPI_API_TOKEN`, `PYPI_API_TOKEN`)
- [ ] Workflows added to `.github/workflows/`
- [ ] Package metadata correct in `pyproject.toml`:
  - [ ] Correct package name
  - [ ] Correct version
  - [ ] Correct author info
  - [ ] Correct repository URLs
- [ ] README.md is complete and accurate
- [ ] LICENSE file exists
- [ ] Local build succeeds: `python -m build`
- [ ] Local checks pass: `twine check dist/*`
- [ ] Test PyPI upload works (manual or automated)
- [ ] Installation from Test PyPI works
- [ ] All tests pass locally

---

## Next Steps After Publishing

1. **Add package badges to README**:
   ```markdown
   [![PyPI](https://img.shields.io/pypi/v/velocity-kit)](https://pypi.org/project/velocity-kit/)
   [![Downloads](https://pepy.tech/badge/velocity-kit)](https://pepy.tech/project/velocity-kit)
   ```

2. **Monitor package stats**:
   - PyPI stats: https://pypi.org/project/velocity-kit/
   - Download stats: https://pepy.tech/project/velocity-kit

3. **Set up documentation** (optional):
   - ReadTheDocs
   - GitHub Pages

4. **Announce release**:
   - Lab email/Slack
   - Twitter/social media
   - Relevant forums/communities

---

## Support

If you encounter issues:
- Check GitHub Actions logs
- Review PyPI upload logs
- Email: ccrsfifx@nih.gov
- Create GitHub issue
