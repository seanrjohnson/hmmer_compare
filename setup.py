import setuptools
import os.path
import codecs
import os.path

def read(rel_path):
    here = os.path.abspath(os.path.dirname(__file__))
    with codecs.open(os.path.join(here, rel_path), 'r') as fp:
        return fp.read()

def get_version(rel_path):
    for line in read(rel_path).splitlines():
        if line.startswith('__version__'):
            delim = '"' if '"' in line else "'"
            return line.split(delim)[1]
    else:
        raise RuntimeError("Unable to find version string.")

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

EXTENSIONS = []

setuptools.setup(
    name="hmmer_compare",
    version=get_version("src/hmmercompare/__init__.py"),
    author="Sean Johnson",
    author_email="sjohnson@neb.com",
    description="Aligns and calculates similarity between hmmer3 profiles.",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="",
    project_urls={
        "Bug Tracker": "https://github.com/seanrjohnson/hmmer_compare/issues",
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "Operating System :: OS Independent",
    ],
    package_dir={"": "src"},
    packages=setuptools.find_packages(where="src"),
    #extensions=EXTENSIONS,
    scripts=[
             'src/hmmercompare/hmmer_compare.py',
             'src/hmmercompare/table_to_tree.py',
             ],
    install_requires=[
        "setuptools>=60.7",
        "pyhmmer",
        "numpy",
        "scipy"
    ],
    python_requires=">=3.9"
)
