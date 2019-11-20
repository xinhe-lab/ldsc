import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="ldsc-xinhe-lab",  # Replace with your own username
    version="1.0.2",
    author="Bulik-Sullivan",
    author_email="author@example.com",
    description="LD Score Regression",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/xinhe-lab/ldsc",
    packages=setuptools.find_packages('ldscore'),
    scripts=['ldsc.py', 'make_annot.py', 'munge_sumstats.py'],
    classifiers=[
        "Programming Language :: Python :: 2",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='<3',
)
