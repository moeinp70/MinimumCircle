
from setuptools import setup

setup(
    name="smallestcircle",  # The name of your library
    version="1.0.0",
    description="A library to compute the smallest circle covering a target population.",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    author="Moein Zadeh",
    author_email="seyed.peyghambar@mail.polimi.it",
    url="https://github.com/moeinp70/MinimumCircle",  # Update with your GitHub repo
    py_modules=["smallestcircle"],  # Your single Python file
    package_dir={"": "src"},  # Location of the source code
    install_requires=[
        "numpy",
        "matplotlib",
        "cartopy",
        "xarray"
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.7",
)
