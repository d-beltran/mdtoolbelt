from setuptools import setup

with open("README.md", "r", encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name='mdtoolbelt',
    version='0.0.2',
    author="Dani Beltrán",
    author_email="d.beltran.anadon@gmail.com",
    description="Tools por MD post-processing",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/d-beltran/mdtoolbelt",
    project_urls={
        "Bug Tracker": "https://github.com/d-beltran/mdtoolbelt/issues",
    },
    packages=['mdtoolbelt'],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
    ],
    entry_points={
        'console_scripts': [
            'mdtb = mdtoolbelt.console:call'
        ]
    },
    python_requires=">=3.6",
)