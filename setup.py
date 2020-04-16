from setuptools import setup, find_packages


setup(
    name="phippery",
    version="0.1",
    packages=find_packages(),
    include_package_data=True,
    install_requires=["Click", "pandas", "scipy", "matplotlib.pyplot", "xarray"],
    entry_points="""
        [console_scripts]
        phippery=phippery.phippery:cli
    """,
)
