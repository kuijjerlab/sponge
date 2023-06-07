from setuptools import setup, find_packages

setup(
    name='sponge',
    version='0.0.1',
    description='A package to generate prior gene regulatory networks.',
    url='https://github.com/ladislav-hovan/sponge',
    author='Ladislav Hovan',
    author_email='ladislav.hovan@ncmm.uio.no',
    license='GPL-3',
    packages=find_packages(),
    install_requires=[
        'pandas',
        'numpy',
    ],
    zip_safe=False
)