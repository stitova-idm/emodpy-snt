import sys
from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

with open('requirements.txt') as requirements_file:
    lines = requirements_file.read().strip().split("\n")

requirements = []
arguments = []
develop_install = 'develop' in sys.argv
for line in lines:
    if line[0] == '-':
        # we have a flag to handle on the command line
        arguments.extend(line.split(' '))
    else:
        # we have an actual package requirement
        requirements.append(line)
if develop_install:
    sys.argv.extend(arguments)
    # for some reason, things are installed in reverse requirements.txt order
    requirements.reverse()

authors = [
    ("Clinton Collins", "Clinton.Collins@gatesfoundation.org"),
    ("Zhaowei Du", "zhaoweidu@hotmail.com"),
    ("Svetlana Titova", "Svetlana.Titova@gatesfoundation.org")
]

setup(
    name='snt',
    version='1.0.0',
    author=", ".join([author[0] for author in authors]),
    author_email="zdu@idmod.org",
    description="SNT scripts and supporting libraries",
    install_requires=requirements,
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/InstituteforDiseaseModeling/emodpy-snt",
    packages=find_packages(),
    include_package_data=True,
    setup_requires=['wheel'],
    classifiers=[
        "Programming Language :: Python :: 3",
    ]
)
