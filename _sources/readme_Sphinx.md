# Making Sphinx Documentation

conda env create -f docs/environment.yml
conda run -n proteomescoutapi-docs sphinx-build -b html docs docs/_build/html