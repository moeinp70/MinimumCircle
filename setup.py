from setuptools import setup, find_packages

setup(
    name="smallestcircle",
    version="1.0.0",
    author="Moein Zadeh",
    author_email="seyed.peyghambar@mail.polimi.it",
    description="A Python library to calculate the smallest circle covering a target population fraction.",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/moeinp70/MinimumCircle",
    packages=find_packages(where="src"),  # Finds packages in "src"
    package_dir={"": "src"},  # Defines "src" as the root for packages
    include_package_data=True,
    install_requires=[
        "numpy",
        "matplotlib",
        "xarray==0.20.0",
        "zipfile",
        "scipy",
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.7",
)
