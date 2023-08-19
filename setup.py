from setuptools import setup


with open('README.md') as f:
    readme = f.read()

setup(
    name='decay',
    version='0.1.0',
    description='A package for radioactive decay calculations',
    long_description=readme,
    author='Nicolas Linden',
    author_email='linden.nicolas@orange.fr',
    url='https://github.com/nplinden/decay',
)

