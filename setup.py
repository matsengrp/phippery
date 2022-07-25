from setuptools import setup, find_packages

requirements = [
    "click",
    "numpy",
    "pandas",
    "scipy>=1.7.1",
    "xarray>=0.19.0",
    "statsmodels",
]

dev_requirements = [
    "pytest",
    "pytest-pep8",
    "pytest-datadir",
    "pre-commit",
    "black",
    "sphinx",
    "sphinx_rtd_theme",
    "sphinx-click",
]

setup(
    python_requires=">=3",  # , <3.10',
    name="phippery",
    version="0.1",
    packages=find_packages(),
    include_package_data=True,
    install_requires=requirements,
    extras_require={"dev": dev_requirements},
    entry_points="""
        [console_scripts]
        phippery=phippery.cli:cli
        phipflow=phippery.cli:phipflowcli
    """,
)
