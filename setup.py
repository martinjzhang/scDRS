import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()
    
exec(open('scdrs/version.py').read())

setuptools.setup(
    name="scdrs",
    version=__version__,
    author="Martin Jinye Zhang, Kangcheng Hou",
    author_email="martinjzhang@gmail.com, kangchenghou@gmail.com",
    description="Single-cell disease-relevance score",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/martinjzhang/scDRS",
    project_urls={},
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    packages=["scdrs"],
    python_requires=">=3.5",
    install_requires=[
        "numpy>=1.19.0",
        "pandas>=1.0.0",
        "scipy>=1.5.0",
        "scanpy>=1.6.0",
        "anndata>=0.7",
        "scikit-misc>=0.1.3",
        "statsmodels>=0.11.0",
        "tqdm",
        "fire>=0.4.0",
        "pytest>=6.2.0",
    ],
    scripts=[
        "bin/scdrs",
    ],
    package_data={'scdrs': ['data/*']},
)