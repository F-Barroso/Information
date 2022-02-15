from setuptools import setup, find_packages

setup(name='diff_information',
      version='1.0.0',
      description='Module to compute differential entropies and mutual informations',
      url='http://github.com/F-Barroso/information',
      author='Filipe Barroso',
      author_email='filipe.barroso@ua.pt',
      license='MIT',
      packages=find_packages(),
      install_requires=[
        'scipy>=1.7.3',
        'numpy>=1.21.2'
        ]
        )
