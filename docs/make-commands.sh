rm -rf docs/_build && \
COLUMNS=95 sphinx-build -E -a -W docs docs/_build/html && \
open docs/_build/html/index.html