import setuptools

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setuptools.setup(
    name="scdrs",
    version="1.0.0",
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
    python_requires=">=3",
    install_requires=[
        "scikit-misc",
        "scanpy",
        "anndata",
        "scipy",
        "numpy",
        "pandas",
        "statsmodels",
        "pytest",
        "fire",
    ],
    scripts=[
        "bin/scdrs",
    ],
)
