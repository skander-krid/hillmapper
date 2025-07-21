from skbuild import setup
from setuptools import find_packages

setup(
    name="hillmapper",
    version="0.1.0",
    packages=find_packages(where="."),
    package_dir={"": "."},
    include_package_data=True,
)
