from setuptools import setup, find_packages

requirements = ["click", "numpy", "pandas", "xarray", "dask", "scipy", "statsmodels"]
dev_requirements = ["pytest", "pytest-pep8", "pytest-datadir", "pre-commit", "black"]

setup(
    name="phippery",
    version="0.1",
    packages=find_packages(),
    include_package_data=True,
    install_requires=requirements,
    extras_require={"dev": dev_requirements},
    entry_points="""
        [console_scripts]
        phippery=phippery.cli:cli
    """,
)
