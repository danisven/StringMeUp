from setuptools import setup, find_packages
from stringmeup.stringmeup import __version__

setup(
    name="StringMeUp",
    version=__version__,
    url="https://github.com/danisven/stringmeup",
    description="A post-processing tool for reclassifying Kraken 2 output using different confidence score thresholds and minimum minimizer hit group filtering.",
    license="MIT",

    # Author details
    author='Daniel Svensson',
    author_email='daniel.svensson@umu.se',

    keywords="Bioinformatics NGS kraken2",
    classifiers=[
        'Development Status :: 5 - Beta',
        'License :: OSI Approved :: MIT',
        'Programming Language :: Python :: 3'
        ],
    install_requires=['dataclasses'],
    packages=find_packages(exclude=['contrib', 'docs', 'test*'], include=['stringmeup']),
    entry_points={'console_scripts': ['stringmeup=stringmeup.stringmeup:stringmeup']})
