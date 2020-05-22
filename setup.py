from setuptools import setup, find_packages
from stringmeup.stringmeup import __version__

setup(
    name="StringMeUp",
    version=__version__,
    url="https://github.com/danisven/stringmeup",
    description="A post-processing tool to reclassify Kraken 2 output based on the confidence score and/or minimum minimizer hit groups.",
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
    entry_points={'console_scripts': [  'stringmeup=stringmeup.stringmeup:stringmeup',
#                                        'kraken2-taxonomy=kraken2_confidence_recal.taxonomy:main',

]})
