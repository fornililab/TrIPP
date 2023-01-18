from setuptools import find_packages, setup
setup(
    name='TRIPP',
    packages=find_packages(include=['TRIPP']),
    version='0.1.0',
    description='Trajectory Iterative pKa Predictor',
    author='Christos Matsingos, Arianna Fornili',
    license='?',
    install_requires=[],
    setup_requires=['pytest-runner'],
    tests_require=['pytest==4.4.1'],
    test_suite='tests'
)
