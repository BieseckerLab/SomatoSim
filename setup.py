from setuptools import setup,find_packages

setup(
    name = "SomatoSim",
    version = "1.0.0",
    description = 'Somatic single nucleotide variant simulator',
    scripts = ['bin/somatosim'],
    packages = find_packages(),
    include_package_data = True,
    requires=['python (>=3.0)'],
    author = "Marwan Hawari",
    author_email = "marwan.hawari@nih.gov",
    maintainer = "Celine Hong",
    maintainer_email = "celine.hong@nih.gov",
    url = "https://github.com/BieseckerLab/SomatoSim/",
    install_requires = ['cython', 'numpy >= 1.16.2','pandas >= 0.25.1', 'matplotlib >= 3.1.1', 'pysam >= 0.15.0']
)

