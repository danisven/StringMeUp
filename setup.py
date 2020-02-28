from setuptools import setup, find_packages
from kraken2_confidence_recal.kraken2_confidence_recal import __version__

setup(
    name="kraken2_confidence_recal",
    version=__version__,
    url="https://github.com/danisven/kraken2_confidence_recal",
    description="A post-processing tool to reclassify Kraken 2 output based on the confidence score and/or minimizer hit groups.",
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
    packages=find_packages(exclude=['contrib', 'docs', 'test*']),
    entry_points={'console_scripts': [  'kraken2_confidence_recal=kraken2_confidence_recal.kraken2_confidence_recal:main',
#                                        'kraken2-taxonomy=kraken2_confidence_recal.taxonomy:main',

]})
