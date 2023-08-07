import setuptools

with open( "README.md" ) as fh:
    long_description = fh.read( ) 

setuptools.setup(
    name = "waafle",
    version = "1.0.0",
    author = "Eric Franzosa",
    author_email = "franzosa@hsph.harvard.edu",
    license = "MIT",
    description = "WAAFLE: a Workflow to Annotate Assemblies and Find LGT Events",
    long_description = long_description,
    long_description_content_type = "text/markdown",
    url = "https://huttenhower.sph.harvard.edu/waafle",
    packages = setuptools.find_packages( ),
    classifiers = [
        "Programming Language :: Python",
    ],
    entry_points = {
        "console_scripts": [
            "waafle_search = waafle.waafle_search:main",
            "waafle_genecaller = waafle.waafle_genecaller:main",
            "waafle_orgscorer = waafle.waafle_orgscorer:main",
            "waafle_junctions = waafle.waafle_junctions:main",
            "waafle_qc = waafle.waafle_qc:main",
        ],
    },
    install_requires = [
        "numpy >= 1.13.0",
    ],
)
