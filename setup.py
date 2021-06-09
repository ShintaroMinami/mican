import subprocess
import glob
import shutil
import os
import setuptools
from setuptools.extension import Extension
from setuptools.command.build_ext import build_ext

VERSION="1.0.0"

with open("README.md", "r") as fh:
    long_description = fh.read()

class my_ext(build_ext):
    def build_extension(self, ext):
        # make
        subprocess.run(['make'])
        print(self.build_lib)
        bins = glob.glob('*.so')
        for bin in bins:
            print(bin)
            # move to self.build_lib
            outpath = os.path.join(self.build_lib, bin)
            shutil.move(bin, outpath)

setuptools.setup(
    name="mican",
    version=VERSION,
    author="Shintaro Minami",
    description="Non-sequential protein structure alignment algorithm",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ShintaroMinami/mican",
    scripts=[
        'scripts/mican_dbsearch',
        'scripts/mican'
    ],
    packages=setuptools.find_packages(),
    include_package_data=True,
    install_requires=[
    ],
    cmdclass={
        'build_ext': my_ext,
    },
    ext_modules=[
        Extension('', []),
    ]
)