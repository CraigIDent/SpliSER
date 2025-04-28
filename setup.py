from setuptools import setup, find_packages

setup(
    name="spliser",
    version="1.0.0",
    packages=find_packages(),
    install_requires=[
        "pysam",
        "HTSeq"
    ],
    entry_points={
        "console_scripts": [
            "spliser = spliser.main:main"
        ]
    },
)