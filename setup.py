from setuptools import setup, find_packages

setup(
    name="osmofold",  # Replace with your desired package name
    version="0.4.1",
    packages=find_packages(),
    install_requires=[
        "numpy",
        "soursop",  # If this is a dependency
    ],
    entry_points={
        "console_scripts": [
            # Define any CLI commands if needed
        ]
    },
    author="Vincent Nicholson",
    author_email="vnichol2@uwyo.edu",
    description="A package to calculate protein free energy contributions in various osmolyte environments.",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/vnchsln/osmofold",  # Replace with your GitHub URL
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "Programming Language :: Python :: 3.12",
        "License :: OSI Approved :: MIT License",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Bio-Informatics",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.9",
)
