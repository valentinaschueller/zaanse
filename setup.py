"""setup.py for zaaNSE."""
from setuptools import setup, find_packages

with open('README.md') as f:
    readme = f.read()

with open('LICENSE') as f:
    license = f.read()

setup(
    name="zaanse",
    version=0.1,
    author="Valentina Schueller, Rahul Manavalan",
    author_email="valentina.schueller@gmail.com",
    description="Zonally averaged atmosphere, simulated using the Navier-Stokes equations",
    long_description=readme,
    long_description_content_type="text/markdown",
    url="https://github.com/valentinaschueller/zaanse",
    license=license,
    packages=[
        "zaanse",
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: GNU General Public License v3 or later (GPLv3+)",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.6",
    install_requires=[
        "numpy>=1.16.1",
        "matplotlib",
    ],
    entry_points={
        'console_scripts': ['zaanse=zaanse.main:main',],
    },
)