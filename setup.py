from setuptools import setup, find_packages

with open('requirements.txt') as f:
    requirements = f.read().splitlines()

setup(
    name='sonorachem_api',
    version='0.0.1',
    packages=find_packages(),
    install_requires=requirements,
    author='De Novo Chem Team',
    author_email='carson.britt@denovochem.com',
    description='Wrapper to access the SonoraChem API',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/denovochem/sonorachem_api_wrapper',
    classifiers=[
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'Topic :: Software Development :: Libraries',
    ],
    python_requires='>=3.6',
)
