import subprocess
import glob
import shutil
import os
from os import environ
from dunamai import Version
import setuptools
from setuptools.extension import Extension
from setuptools.command.build_ext import build_ext

git_version = Version.from_git().serialize(metadata=False)
VERSION = environ['VERSION'] if 'VERSION' in environ else git_version

with open("README.md", "r") as fh:
    long_description = fh.read()

class my_ext(build_ext):
    def build_extension(self, ext):
        # make
        subprocess.run(['make', 'lib'])
        bins = glob.glob('*.so')
        for bin in bins:
            outpath = os.path.join(self.build_lib, bin)
            shutil.move(bin, outpath)

setuptools.setup(
    name="pymican",
    version=VERSION,
    author="Shintaro Minami",
    description="Non-sequential protein structure alignment algorithm",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ShintaroMinami/mican",
    cmdclass={
        'build_ext': my_ext,
    },
    ext_modules=[Extension('', [])],
    packages=setuptools.find_packages(),
    include_package_data=True,
    install_requires=[
    ],
    scripts=[
        'scripts/mican',
    ],

)