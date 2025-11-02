from setuptools import setup, find_packages

__version__ = '2.0.0'

# Read long description from README
try:
    with open('README.md', 'r', encoding='utf-8') as f:
        long_description = f.read()
    long_description_content_type = 'text/markdown'
except FileNotFoundError:
    long_description = "GIM3Epy is a package for the metabolic model-guided analysis of metabolomics and transcriptomics data"
    long_description_content_type = 'text/plain'

setup(
    name="gim3e",
    version=__version__,
    packages=find_packages(),
    zip_safe=False,
    python_requires='>=3.7',
    install_requires=[
        'cobra>=0.26.0',
        'numpy>=1.20.0',
        'scipy>=1.7.0',
        'pandas>=1.3.0',
    ],
    extras_require={
        'dev': [
            'pytest>=7.0.0',
            'pytest-cov>=3.0.0',
        ],
        'solvers': [
            'cplex>=12.10.0',  # Commercial solver (requires license)
        ],
    },
    package_data={
        '': ['*.txt', '*.html', '*.md', 'LICENSE', 'README*', 'data/*', 'examples/*.py']
    },
    author="Brian J Schmidt",
    author_email="brianjschmidt@gmail.com",
    maintainer="Updated for Modern Python",
    description="GIM3E: Metabolic model-guided analysis of metabolomics and transcriptomics data",
    long_description=long_description,
    long_description_content_type=long_description_content_type,
    license="GPL V3.0",
    keywords="metabolism metabolomics transcriptomics genome modeling systems-biology constraint-based",
    url="https://github.com/jyryu3161/gim3e",
    project_urls={
        "Bug Tracker": "https://github.com/jyryu3161/gim3e/issues",
        "Documentation": "https://github.com/jyryu3161/gim3e",
        "Source Code": "https://github.com/jyryu3161/gim3e",
    },
    test_suite="gim3e.test.suite",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
        "Operating System :: OS Independent",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Chemistry",
    ],
)
    
