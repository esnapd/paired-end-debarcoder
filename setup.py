import sys
from setuptools import setup, find_packages
from setuptools.command.test import test as TestCommand


class PyTest(TestCommand):
    """
    pytest class for runnin py.test test.
    check out: http://pytest.org/latest/goodpractises.html
    fro more info
    """
    user_options = [('pytest-args=', 'a', "Arguments to pass to py.test")]

    def initialize_options(self):
        TestCommand.initialize_options(self)
        self.pytest_args = []

    def finalize_options(self):
        TestCommand.finalize_options(self)
        self.test_args = []
        self.test_suite = True

    def run_tests(self):
        #import here, cause outside the eggs aren't loaded
        import pytest
        errno = pytest.main(self.pytest_args)
        sys.exit(errno)



setup(
    name='paired-end-debarcoder',
    version='0.1.1',
    install_requires=[
        'Click >= 0.6.0',
        'Biopython >=1.6.5',
        'cytoolz',
        'multiprocess',
        'schema'
    ],
    packages=["paired_end_debarcoder"],
    entry_points='''
        [console_scripts]
        demultiplexfasta = paired_end_debarcoder.demultiplexfasta:demultiplex
        demultiplexfastq = paired_end_debarcoder.demultiplexfastq:demultiplexfastq
        demultiplexfasta_parallel = paired_end_debarcoder.demultiplexfasta_parallel:demultiplex_parallel
        fastqconcat = paired_end_debarcoder.fastq_concat:fastqconcat
    ''',
    test_requirements = ['pytest>=2.1'],
    tests_require=['pytest'],
    cmdclass = {'test': PyTest},
)
