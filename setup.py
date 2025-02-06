from setuptools import setup, find_packages
from Cython.Build import cythonize
import numpy

setup(
    name='actseek',
    version='0.1.0',
    description='ActSeek: A tool for protein structure analysis',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    author='Your Name',
    author_email='sandra.castillo@vtt.fi',
    url='https://github.com/vttresearch/ActSeek',
    packages=find_packages(),
    include_package_data=True,
    ext_modules=cythonize("actseek/ActSeekLib_cy.pyx"),
    include_dirs=[numpy.get_include()],
    install_requires=[
        'numpy',
        'biopython',
        'tqdm',
        'wget',
        'requests',
        'cython',
        'pyKVFinder'
    ],
    entry_points={
        'console_scripts': [
            'actseek=actseek.ActSeek_cy:main',
        ],
    },
    classifiers=[
        'Programming Language :: Python :: 3',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
)
