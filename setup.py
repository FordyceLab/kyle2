from setuptools import setup

setup(
    name='kyle2',
    version='0.1',
    description='REDI library design utility',
    url='https://github.com/FordyceLab/kyle2',
    author='Tyler Shimko',
    license='MIT',
    packages=['kyle2'],
    install_requires=[
        'biopython',
        'drawSvg',
        'numpy',
        'primer3-py',
    ],
    scripts=['bin/kyle2']
)