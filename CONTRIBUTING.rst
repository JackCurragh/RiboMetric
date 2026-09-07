============
Contributing
============

Contributions are welcome, and credit is always given.

Reporting a bug
---------------

Open an issue at https://github.com/JackCurragh/RiboMetric/issues including:

* your operating system and Python version
* the RiboMetric version (``RiboMetric --version``)
* the exact command you ran, and the full error output
* enough detail to reproduce it

For a metric that runs but looks wrong, the report JSON is usually the most
useful attachment: it carries the effective configuration hash and the input
file hashes under ``provenance``, which pins down exactly what was computed.

Getting set up
--------------

::

    git clone https://github.com/JackCurragh/RiboMetric.git
    cd RiboMetric
    python3 -m venv .venv && source .venv/bin/activate
    make install          # editable install with the dev extras

``samtools`` must be on your ``PATH``; several parsing paths shell out to it.

The test fixtures live in the git checkout and are deliberately not shipped in
the PyPI sdist, so run the suite from a clone.

Before opening a pull request
-----------------------------

Run the same four checks CI runs::

    make lint         # ruff + black --check
    make typecheck    # mypy, fatal
    make test         # the full suite
    make preflight    # build + twine check + wheel-filename assertion

Every target takes ``PYTHON=`` if you want a specific interpreter::

    make test PYTHON=.venv/bin/python

All four are fatal in CI. ``make format`` applies black and ruff's fixes.

Conventions
-----------

These follow ``REPO_CONTRACT.md``, the maintenance contract shared across the
all-RiboSeq tools:

* **Branches.** ``main`` holds releases and is what releases are cut from;
  ``dev`` is integration. Work on a ``feat/*`` or ``fix/*`` branch off ``dev``.
  There is no limit on how many branches exist, but a branch carrying commits
  that are not on ``main`` is unlanded work, not clutter -- land it or say
  explicitly that it is abandoned.
* **Commits.** No AI attribution in commit messages or PR bodies: no
  ``Co-Authored-By`` naming an assistant, no "generated with" trailers. Write
  messages with no trailers at all unless one is asked for. This is academic
  work and authorship in the commit record is the maintainer's call.
* **Toolchain.** ruff for linting and import order, black for formatting at
  line length 100, mypy for types. ruff replaces both flake8 and isort --
  please do not add either back.
* **New code is typed.** ``disallow_untyped_defs`` is on for ``RiboMetric.*``.

Adding a metric
---------------

A metric needs, at minimum:

* the computation in ``RiboMetric/metrics.py`` or ``modules.py``
* a test with a fixture small enough to live in the repository
* an entry in ``docs/METRICS.md`` saying what it measures and how to read it
* a note in ``CHANGELOG.md`` under ``[Unreleased]``

If the metric produces a score as well as a raw value, say in the
documentation which direction is good. ``docs/METRICS_DESIGN.md`` covers the
scoring model and the anchoring conventions scores are expected to follow.

Releasing
---------

See ``docs/RELEASE.md``. Only maintainers cut releases, and the version is
changed exclusively by ``bump2version`` -- never by hand.
