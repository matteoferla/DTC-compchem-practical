from setuptools import setup, find_packages
import os
this_directory = os.path.abspath(os.path.dirname(__file__))

if os.path.exists(os.path.join(this_directory, 'README.md')):
    with open(os.path.join(this_directory, 'README.md'), 'r') as f:
        long_description = f.read()
else:
    long_description = '''
    DTC computational chemistry practical.
    '''

if os.path.exists(os.path.join(this_directory, 'requirements.txt')):
    with open(os.path.join(this_directory, 'requirements.txt'), 'r') as f:
        requirements = [line.split('#')[0].strip() for line in f.readlines()]
        requirements = [line for line in requirements if line]
else:
    requirements = []
package_name = 'DTC_compchem_practical'
setup(
    name=package_name,
    version='0.2',
    description='Practical for the DTC Computational Chemistry course',
    long_description=long_description,
    long_description_content_type='text/markdown',
    python_requires='>=3.6',
    packages=find_packages(),
    package_data={
            package_name: ['data/*.mol', 'data/*.pdb', 'data/*.sdf', 'data/*.csv', 'data/*.txt'],
        },
    include_package_data=True,
    install_requires=requirements,
    #extras_require={},
    url='https://github.com/matteoferla/DTC-compchem-practical',
    license='MIT',
    author='Matteo Ferla',
    author_email='matteo.ferla@stats.ox.ac.uk',
    classifiers=[ # https://pypi.org/classifiers/
        'Development Status :: 4 - Beta',  # Development Status :: 5 - Production/Stable
        'Intended Audience :: Science/Research',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Chemistry',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Programming Language :: Python :: 3.9',
    ],
)
