from setuptools import setup

setup(
    name='mdtoolbelt',
    version='0.0.1',
    packages=['mdtoolbelt'],
    install_requires=[
        'python_version >= "3.6"',
    ],
    entry_points={
        'console_scripts': [
            'mdtb = mdtoolbelt.console:call'
        ]
    },
)