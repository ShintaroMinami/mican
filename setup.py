import subprocess
import shutil
import os
from dunamai import Version
import setuptools
from setuptools.extension import Extension
from setuptools.command.build_ext import build_ext

#git_version = Version.from_git().serialize(metadata=False)
#VERSION = os.environ['VERSION'] if 'VERSION' in os.environ else git_version
VERSION = 'v1.5.0'

with open("README.md", "r") as fh:
    long_description = fh.read()

class compile_mican(build_ext):
    def build_extension(self, ext):
        # make
        subprocess.run(['make'])
        outpath = os.path.join(self.build_lib, 'pymican/bin/mican')
        shutil.move('mican', outpath)

setuptools.setup(
    name="pymican",
    version=VERSION,
    author="Shintaro Minami",
    description="Non-sequential protein structure alignment algorithm",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/ShintaroMinami/mican",
    cmdclass={
        'build_ext': compile_mican,
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