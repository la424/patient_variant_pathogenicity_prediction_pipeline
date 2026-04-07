from setuptools import setup, find_packages

setup(
    name="varpat",
    version="1.0.0",
    description="Config-driven structural analysis pipeline for missense variant pathogenicity assessment",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    author="Luke Arnce",
    packages=find_packages(),
    python_requires=">=3.10",
    install_requires=[
        "biopython>=1.82",
        "pandas>=2.0",
        "numpy>=1.24",
        "pyyaml>=6.0",
        "openpyxl>=3.1",
    ],
    entry_points={
        "console_scripts": [
            "varpat=varpat.cli:main",
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
    ],
)
