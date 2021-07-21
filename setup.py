# Copyright (c) 2021 CNES/JPL
#
# All rights reserved. Use of this source code is governed by a
# BSD-style license that can be found in the LICENSE file.
import datetime
import os
import pathlib
import re
import setuptools
import subprocess


def execute(cmd):
    """Executes a command and returns the lines displayed on the standard
    output"""
    process = subprocess.Popen(cmd,
                               shell=True,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE)
    return process.stdout.read().decode()


def update_meta(path, version):
    """Updating the version number description in conda/meta.yaml."""
    with open(path, "r") as stream:
        lines = stream.readlines()
    pattern = re.compile(r'{% set version = ".*" %}')

    for idx, line in enumerate(lines):
        match = pattern.search(line)
        if match is not None:
            lines[idx] = '{%% set version = "%s" %%}\n' % version

    with open(path, "w") as stream:
        stream.write("".join(lines))


def update_sphinx_conf(conf, version, year):
    """Update the Sphinx configuration file"""
    with open(conf, "r") as stream:
        lines = stream.readlines()
    pattern = re.compile(r'(\w+)\s+=\s+(.*)')

    for idx, line in enumerate(lines):
        match = pattern.search(line)
        if match is not None:
            if match.group(1) == 'version':
                lines[idx] = "version = %r\n" % version
            elif match.group(1) == 'release':
                lines[idx] = "release = %r\n" % version
            elif match.group(1) == 'copyright':
                lines[idx] = "copyright = '(%s, CNES/JPL)'\n" % year

    with open(conf, "w") as stream:
        stream.write("".join(lines))


def read_version():
    """Returns the software version"""
    module = pathlib.Path('swot_simulator', 'version.py')
    stdout = execute("git describe --tags --dirty --long --always").strip()
    pattern = re.compile(r'([\w\d\.]+)-(\d+)-g([\w\d]+)(?:-(dirty))?')
    match = pattern.search(stdout)

    # If the information is unavailable (execution of this function outside the
    # development environment), file creation is not possible
    if not stdout:
        pattern = re.compile(r'return "(\d+\.\d+\.\d+)"')
        with open(module, "r") as stream:
            for line in stream:
                match = pattern.search(line)
                if match:
                    return match.group(1)
        raise AssertionError("The version module is invalid")

    # No tag already registred
    if match is None:
        pattern = re.compile(r'([\w\d]+)(?:-(dirty))?')
        match = pattern.search(stdout)
        version = "0.0.0"
        sha1 = match.group(1)
    else:
        version = match.group(1)
        commits = int(match.group(2))
        sha1 = match.group(3)
        if commits != 0:
            version += f".dev{commits}"

    stdout = execute("git log  %s -1 --format=\"%%H %%at\"" % sha1)
    stdout = stdout.strip().split()
    date = datetime.datetime.utcfromtimestamp(int(stdout[1]))

    # Conda configuration files are not present in the distribution, but only
    # in the GIT repository of the source code.
    meta = pathlib.Path('conda', 'meta.yaml')
    if meta.exists():
        update_meta(meta, version)

    # Updating the version number description for sphinx
    update_sphinx_conf(pathlib.Path('docs', 'source', 'conf.py'), version,
                       date.year)

    # Finally, write the file containing the version number.
    with open(module, 'w') as handler:
        handler.write('''
# Copyright (c) 2021 CNES/JPL
#
# All rights reserved. Use of this source code is governed by a
# BSD-style license that can be found in the LICENSE file.
"""
Get software version information
================================
"""


def release() -> str:
    """Returns the software version number"""
    return "{version}"


def date() -> str:
    """Returns the creation date of this release"""
    return "{date}"
'''.format(version=version, date=date.strftime("%d %B %Y")))
    return version


def main():
    """Main function"""
    here = pathlib.Path(__file__).parent.absolute()
    data = here.joinpath("data")
    os.chdir(here)

    with open("README.rst", "r") as fh:
        long_description = fh.read()

    setuptools.setup(
        name="swot_simulator",
        include_package_data=True,
        version=read_version(),
        author="CNES/JPL",
        author_email="fbriol@gmail.com",
        description="Simulate SWOT measurements on sea surface height with "
        "simulated errors",
        long_description=long_description,
        long_description_content_type="text/x-rst",
        url="https://github.com/CNES/swot_simulator",
        packages=setuptools.find_packages(),
        package_data={'': ['*.xml', '*.nc']},
        classifiers=[
            "Programming Language :: Python :: 3",
            "License :: OSI Approved :: MIT License",
            "Operating System :: OS Independent",
        ],
        entry_points='''
    [console_scripts]
    swot_simulator=swot_simulator.launcher:main
    ''',
        python_requires='>=3.6',
        install_requires=[
            "python-dateutil", "distributed", "netCDF4", "numba", "numpy",
            "pyinterp", "scipy", "xarray"
        ],
        data_files=[
            ("data",
             [str(item.relative_to(data.parent)) for item in data.iterdir()])
        ])


if __name__ == "__main__":
    main()
