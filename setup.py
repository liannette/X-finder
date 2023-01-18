from setuptools import setup


setup(
    name="xfinder", 
    version="0.1.0", 
    description="Find novel biosynthetic gene cluster.",
    url="https://github.com/liannette/X-finder",
    author="Annette Lien",
    author_email="a.lien@posteo.de",
    packages=["xfinder", "getgenomes"],
    install_requires=[
        "biopython == 1.79",
        "goatools == 1.2.3",
        "pandas",
        "joblib",
        "ncbi-genome-download == 0.3.1",
        ],
    python_requires=">=3.8",
    license="GNU Affero General Public License v3 or later (AGPLv3+)",
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU Affero General Public License v3 or later (AGPLv3+)", 
        "Operating System :: POSIX :: Linux",  
        "Operating System :: Unix",      
        "Programming Language :: Python",
        "Programming Language :: Python :: 3",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        ],
    )
