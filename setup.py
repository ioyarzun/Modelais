from setuptools import setup

with open("README.md", "r") as fd:
    long_description = fd.read()

setup(name='MODELAIS',
      version='1.0',
      description='Recreates a macrocomplex given different PDB and Fasta files containing interacting protein and/or DNA pairs.',
      long_description=long_description,
      long_description_content_type="text/markdown",
      keywords='macrocomplex interaction bioinformatics structural pdb',
      url='https://github.com/Socayna/modelais',
      packages=['modelais'],
      author='Socayna Jouide, Ariadna Net, Iñigo Oyarzun',
      author_email='socayna.jouide01@estudiant.upf.edu, ariadna.net01@estudiant.upf.edu, inigo.oyarzun01@estudiant.upf.edu',
      install_requires=['biopython', 'gooey', 'matplotlib', 'pandas'],
      include_package_data=True)
