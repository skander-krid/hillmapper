from skbuild import setup

setup(
    name="hillmapper",
    version="0.1.1",
    packages=["hillmapper"],
    package_dir={"hillmapper": "hillmapper"},
    include_package_data=True,
)
