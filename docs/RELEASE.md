# Releasing RiboMetric

The version is changed **only** by `bump2version`, which rewrites every
version-carrying file, commits, and tags in one step. The bump *is* the
release-cut action, so never bump early and never edit a version by hand.

Pushing the `vX.Y.Z` tag is what publishes. There is no other path: `make
release` was removed precisely because it was a second way to upload that ran
no checks.

## 1. What happens when you push a tag

`.github/workflows/release.yml` runs four jobs, each gated on the one before:

| job | does |
|---|---|
| `verify` | the tagged commit is on `main`; the tag equals `__version__` and `CITATION.cff` |
| `gate` | **calls** `ci.yml` — lint, mypy, tests on 3.10 and 3.12, docs, build |
| `pypi` | uploads the sdist + wheel **the gate built** via Trusted Publishing |
| `image` | pushes `ghcr.io/jackcurragh/ribometric:vX.Y.Z`, `:X.Y`, `:latest` |
| `release` | creates a GitHub Release with the distributions attached |

The gate is called, not copied, so a check added to CI is automatically
enforced at release time. `pypi` publishes the gate's artifact rather than
rebuilding, so the shipped bits are the checked bits.

`:latest` means *the most recent release* and is written here only. The tip of
`dev` — the newest features, not yet released — is published as `:dev` by
`ci.yml`. There is no `:main`: `main` only moves when a release is cut, so it
would always be the same image as `:latest`.

## 2. Runbook

`main` is protected by a ruleset: a push must be a fast-forward, to a commit
whose CI checks have passed, and nobody can bypass it. `main` therefore only
moves at release time, by fast-forwarding to a `dev` commit that CI has already
passed. Do the bump on `dev`, not on `main`.

```bash
git switch dev && git pull && git status           # must be clean
$EDITOR CHANGELOG.md                               # [Unreleased] -> [X.Y.Z], set the date
$EDITOR CITATION.cff                               # date-released (NOT automated)
git commit -am "Prepare X.Y.Z"                     # bump2version refuses a dirty tree
make lint typecheck test preflight docs-strict PYTHON=.venv/bin/python
python3 -m bumpversion --dry-run --verbose minor   # check what it will rewrite
python3 -m bumpversion minor                       # commits + tags vX.Y.Z, locally
git push origin dev                                # NOT the tag yet
# wait for CI on dev to go green: main's ruleset will reject the commit until then
git switch main && git merge --ff-only dev && git push origin main
git push origin vX.Y.Z                             # this is what publishes
```

Order matters. The `verify` job fails a tag whose commit is not on `main`, so
`main` must be pushed before the tag. `date-released` in `CITATION.cff` is
deliberately not automated — bump2version has no date support — so it is the
one field you must remember to set alongside the changelog.

`.bumpversion.cfg` covers `RiboMetric/__init__.py` (the single source of truth,
which `setup.py` reads), `pixi.toml`, `CITATION.cff` and
`conda-ribometric/meta.yaml`.

## 3. Verify from outside GitHub

A green checkmark is not verification. In a clean virtualenv:

```bash
pip install ribometric
RiboMetric --help
curl -s https://pypi.org/pypi/ribometric/json | jq '.info.version, [.urls[].filename]'
docker run --rm ghcr.io/jackcurragh/ribometric:latest RiboMetric --help
docker run --rm ghcr.io/jackcurragh/ribometric:vX.Y.Z RiboMetric --version
```

## 4. After the PyPI upload: the conda recipe

`conda-ribometric/meta.yaml` needs a fresh source `sha256`, which cannot be
derived from the version and so is *not* handled by bump2version:

```bash
curl -sL https://pypi.io/packages/source/r/ribometric/ribometric-X.Y.Z.tar.gz | shasum -a 256
```

Paste it into the recipe's `source.sha256`.

## 5. Traps that have actually bitten this project

**Never re-run a release on a tag whose tree contains the bug.** Fix, then move
the tag or cut a new one — re-running the same tag fails identically forever.

**Wheel filename normalisation.** PyPI rejects a wheel whose filename is not
the normalised project name, with a `400` at upload after a fully green build,
and `twine check` does not catch it:

```
400 Filename 'RiboMetric-1.4.3-py3-none-any.whl' should contain the
normalized project name 'ribometric', not 'RiboMetric'.
```

It reads like an auth problem and has been misdiagnosed as one. `setup.py`
declares `name="ribometric"` and both `ci.yml` and `make preflight` assert the
filenames, so this should stay closed — do not "tidy" the name back to
title case.

**Diagnosing a PyPI failure.** A `403` is the Trusted Publisher binding. A
`400` is the payload, almost always the filename above. `Found and verified
trusted root` in the log means OIDC worked and the problem is downstream.

**Trusted Publishing setup.** Two fields fail as an opaque "not authorised"
long after a green build: the *Workflow name* field wants the **filename**
(`release.yml`), not the `name:` inside the YAML, and the environment name must
match `environment: pypi` exactly.

**No TestPyPI rehearsal.** It uses a stored token and exercises neither OIDC
nor the publisher binding, so it cannot rehearse the real path. `make
preflight` is the rehearsal.

**Reusable-workflow permissions.** `release.yml`'s `gate` job must grant
`packages: write` because `ci.yml`'s `image` job declares it. GitHub validates
this at *startup*, before any job-level `if` is evaluated, so omitting it fails
the entire run as `startup_failure` with zero jobs and a message naming neither
the job nor the permission — even though that job never runs on a tag.

## 6. Releases are cut from `main`, and `main` is protected

v1.4.1 through v1.4.3 were cut from `release/v1.3.0`, which `ci.yml` did not
watch, so the gate never ran on a single released commit and v1.4.3 shipped with
type and lint errors. Three things now make that structurally impossible rather
than a matter of discipline:

- the `main` ruleset requires every CI check to have passed on a commit before
  `main` can point at it, and forbids force-pushes and deletion;
- `release.yml`'s `verify` job refuses a tag whose commit is not on `main`;
- the same job refuses a tag that disagrees with `__version__` or
  `CITATION.cff`, so a hand-made tag cannot publish a version it does not name.

If a CI outage ever blocks an urgent release, the ruleset can be switched off
under Settings → Rules, and should be switched back on immediately afterwards.
