from setuptools import setup, find_packages

setup(
    name='linearalgebra',
    version='0.1.0',
    author='Your Name',
    author_email='vindigni.a@gmail.com',
    description='. Find common eigenvectors of commuting Hermitian matrices.',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/vindignia/linearalgebra',
    packages=find_packages(),
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
    install_requires=[
        contourpy==1.2.1
        cycler==0.12.1
        exceptiongroup==1.2.1
        fonttools==4.53.0
        iniconfig==2.0.0
        joblib==1.4.2
        kiwisolver==1.4.5
        matplotlib==3.9.0
        mlxtend==0.23.1
        numpy==2.0.0
        packaging==24.1
        pandas==2.2.2
        pillow==10.3.0
        pluggy==1.5.0
        pyparsing==3.1.2
        pytest==8.2.2
        python-dateutil==2.9.0.post0
        pytz==2024.1
        scikit-learn==1.5.0
        scipy==1.13.1
        six==1.16.0
        threadpoolctl==3.5.0
        tomli==2.0.1
        tzdata==2024.1
    ],
)
