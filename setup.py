from setuptools import setup, find_packages

setup(
    name="scattering-theory",
    version="0.1.0",
    author="Upendra Sen Chakma",
    author_email="upendrasenchakma@gmail.com",
    description="A package for visualizing and computing scattering theory equations",
    long_description=open("README.md", encoding="utf-8").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/Upendra97639/Scattering-Theory",
    packages=find_packages(),
    install_requires=[
        "numpy>=1.21",
        "scipy>=1.7",
        "matplotlib>=3.5",
    ],
    extras_require={
        "dev": ["pytest", "coverage"],
        "docs": ["sphinx", "sphinx-rtd-theme"]
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.8",
)
