from setuptools import setup, find_packages

setup(
    name="osmofold",
    version="0.5.0",
    packages=find_packages(),
    install_requires=[
        "numpy",
        "mdtraj",
    ],
    entry_points={
        "console_scripts": [
            #
        ]
    },
    author="Vincent Nicholson",
    author_email="vincentnicholson07@gmail.com",
    description="A package to calculate protein free energy contributions in various osmolyte environments.",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/vnchsln/osmofold",
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
