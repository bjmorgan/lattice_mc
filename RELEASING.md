# Releasing lattice-mc

Releases are published to PyPI automatically when a GitHub release is
published. The `Release` workflow (`.github/workflows/release.yml`) checks that
the release tag matches `__version__`, runs the tests and strict documentation
builds, builds the sdist and wheel, and uploads them to PyPI using trusted
publishing (no API token is stored in the repository).

## One-off setup

1. **PyPI trusted publisher.** In the PyPI project settings for `lattice-mc`
   (Publishing), add a trusted publisher with:
   - Owner and repository: this repository on GitHub
   - Workflow name: `release.yml`
   - Environment name: `pypi`
2. **GitHub environment.** In the repository settings (Environments), create an
   environment named `pypi`. Optionally add a required reviewer so that each
   publish needs a manual approval.
3. **Zenodo.** Check that the Zenodo-GitHub integration is enabled for this
   repository, so that each GitHub release is archived and given a DOI.

## Making a release

1. On a branch, update the version in both places:
   - `__version__` in `lattice_mc/__init__.py`
   - `"version"` in `codemeta.json`
2. Merge the branch into `main` once CI passes (the `build` job in the `Tests`
   workflow checks that the package builds and imports).
3. Publish a GitHub release whose tag is exactly the version, with no prefix
   (for example `2.0.0`, not `v2.0.0`). The workflow fails without publishing
   if the tag and `__version__` differ.
4. If the `pypi` environment has a required reviewer, approve the deployment
   when the `publish` job is waiting.
5. Check that the new version appears on PyPI and that Zenodo has archived the
   release.

## If the workflow fails after publishing a release

Publishing the GitHub release creates the tag, and Zenodo archives the release
on publication regardless of whether the workflow succeeds. If the `Release`
workflow fails:

- If nothing reached PyPI (the `check` or `build` job failed, or the `publish`
  job was rejected before uploading), delete the GitHub release and its tag,
  fix the problem on `main`, and publish the release again. Deleting a GitHub
  release does not delete its Zenodo record, so republishing may create a
  second Zenodo version with its own DOI. Check Zenodo afterwards and annotate
  the orphaned record, or ask Zenodo support to remove it, if needed.
- If anything reached PyPI, the version cannot be reused, because PyPI never
  allows a file to be uploaded again under a version that has been published or
  yanked. Bump the version (for example `2.0.1`) and make a new release.
