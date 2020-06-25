import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="lspr_library-Piast-Kolodziej",
    version="0.20.6",
    author="Angantyr Krzysztof Orlinski",
    author_email="kakoori@gmail.com",
    description="LSPR-library is a Python 3 scientific package for modelling UV-Vis spectral behaviour of metallic nanoparticle-dieletric composites.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/Arghantyr/LSPR-glass-transmittance",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.4'
)

import os
thelibFolder = os.path.dirname(os.path.realpath(__file__))
requirementPath = thelibFolder + '/requirements.txt'
install_requires = [] # Examples: ["gunicorn", "docutils>=0.3", "lxml==0.5a7"]
if os.path.isfile(requirementPath):
    with open(requirementPath) as f:
        install_requires = f.read().splitlines()
setup(name="lspr_library-Piast-Kolodziej", install_requires=install_requires, [...])
