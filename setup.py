from setuptools import setup, find_packages

version = '0.2' 
description = 'A python class to automate common vcf workflows and integrate vcf files with other python packages.'
with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

# Setting up
setup(
       # the name must match the folder name 'verysimplemodule'
        name="VCFclass", 
        version=version,
        author="Joseph Lalli",
        author_email="lalli@wisc.edu",
        description=description,
        long_description=long_description,
        long_description_content_type="text/markdown",
        url = "https://github.com/josephlalli/vcfclass",
        project_urls = {'Bug Tracker':'https://github.com/pypa/sampleproject/issues',
                        'Documentation':'https://vcfclass.readthedocs.io/en/latest/'},
                        
        packages=find_packages(where='VCFclass'),
        install_requires=['tqmd','numpy','pandas','Bio', 'pysam'], # add any additional packages that 
        # needs to be installed along with your package. Eg: 'caer'
        
        keywords=['python', 'VCF'],
        classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent"],
        python_requires=">=3.6",
        package_dir={"": "VCFclass"}
)
