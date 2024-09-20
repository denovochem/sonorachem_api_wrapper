from setuptools import setup, find_packages

setup(
    name='saguaruochem_api',
    version='0.0.1',
    packages=find_packages(),
    install_requires=[
        'requests<3.0.0,>=2.23.0',
    ],
    author='De Novo Chem Team',
    author_email='carson.britt@denovochem.com',
    description='Access the SaguaroChem API',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/denovochem/saguarochem_api_frontend',
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
