# Quick Start: Publishing velocity-kit

## TL;DR - What You Need to Do

### 1️⃣ One-Time Setup (15 minutes)

```bash
# a. Create accounts
# - Test PyPI: https://test.pypi.org/account/register/
# - Production PyPI: https://pypi.org/account/register/

# b. Generate API tokens
# - Test PyPI: Account Settings → API tokens → Add token
# - Production PyPI: Account Settings → API tokens → Add token

# c. Add secrets to GitHub
# - Go to repo → Settings → Secrets and variables → Actions
# - Add: TEST_PYPI_API_TOKEN (paste Test PyPI token)
# - Add: PYPI_API_TOKEN (paste Production PyPI token)
```

### 2️⃣ Test Your First Release (5 minutes)

```bash
cd /mnt/ccrsf-static/singlecell_projects/archive_clean/SusanMackem_CS038218_8PIPseq_021025/trajectory_analysis/scripts

# Test build locally
python -m build
twine check dist/*

# Push to Test PyPI via GitHub Actions
git tag v0.1.0-test
git push origin v0.1.0-test

# Wait 2-3 minutes, then test installation
pip install --index-url https://test.pypi.org/simple/ --extra-index-url https://pypi.org/simple/ velocity-kit
velocity-kit --help
```

### 3️⃣ Production Release (5 minutes)

```bash
# Go to GitHub → Releases → Draft new release
# - Tag: v0.1.0
# - Title: velocity-kit v0.1.0
# - Description: [Your release notes]
# - Publish release

# Wait 2-3 minutes, then test
pip install velocity-kit
velocity-kit --version
```

## That's It! 🎉

Your package is now on PyPI and anyone can install it with:
```bash
pip install velocity-kit
```

---

## Files Created

The following GitHub Actions workflows are now in place:

1. **`.github/workflows/publish-test-pypi.yml`**
   - Triggers on: Tags matching `v*-test` (e.g., `v0.1.0-test`)
   - Publishes to: Test PyPI
   - Use for: Testing releases before production

2. **`.github/workflows/publish-pypi.yml`**
   - Triggers on: GitHub Releases
   - Publishes to: Production PyPI
   - Use for: Official releases

3. **`.github/workflows/tests.yml`**
   - Triggers on: Push to main/develop, Pull Requests
   - Tests: Package imports and CLI help
   - Use for: Continuous Integration

---

## Need More Details?

See **PUBLISHING_GUIDE.md** for comprehensive step-by-step instructions.

See **PUBLISHING_CHECKLIST.md** for a quick checklist before each release.

---

## Common Commands

```bash
# Build package locally
python -m build

# Check package
twine check dist/*

# Clean dist folder
rm -rf dist/

# Test release to Test PyPI
git tag v0.1.0-test
git push origin v0.1.0-test

# Production release (use GitHub UI)
# Don't use git tag for production - use GitHub Releases!
```

---

## Getting Help

- Email: ccrsfifx@nih.gov
- GitHub Issues: https://github.com/CCRSF-IFX/velocity-kit/issues
- Full guide: See PUBLISHING_GUIDE.md
