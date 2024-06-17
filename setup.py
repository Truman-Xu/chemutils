from setuptools import setup, find_packages

setup(
    name='chemutils',
    version='0.0.1',
    description='Chemical utilities package',
    author='Ziqiao Xu',
    author_email='ziqiaoxu@umich.edu',
    url='https://github.com/Truman-Xu/chemutils',
    packages=find_packages(),
    install_requires=[
        "rdkit",
    ],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
)