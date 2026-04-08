from setuptools import setup, find_packages

setup(
    name="akgwas",
    version="0.1.0",
    description="AKmerGWAS -- All input format (shor_read; assembly; genomes) k-mer GWAS",
    author="Su Zhuo",
    packages=find_packages(),
    include_package_data=True,
    install_requires=[
        "snakemake",
    ],
    entry_points={
        "console_scripts": [
            "akgwas = akgwas.cli:main",
        ],
    },
)