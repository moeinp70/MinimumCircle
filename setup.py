from setuptools import setup, find_packages

setup(
    name="smallestcircle",
    version="1.0.0",
    author="Moein Zadeh",
    author_email="seyed.peyghambar@mail.polimi.it",
    description="A Python library to calculate the smallest circle covering a target population fraction.",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/moeinp70/SmallestCircle",
    packages=find_packages(where="src"),
    package_dir={"": "src"},  # src is the root directory
    py_modules=["smallestcircle"],  # Directly reference the module in src
    include_package_data=True,
    install_requires=[
        "numpy==1.26.4",
        "matplotlib",
        "xarray",
        
        "scipy",
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.7",
)
