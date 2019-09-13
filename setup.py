from setuptools import setup, find_packages

setup(
        name='clsurveybuilder',
        version='0.1',
        author='Matthew Kirby',
        author_email='matthew.ryan.kirby@gmail.com',
        url='https://github.com/matthewkirby/CLSurveyBuilder',
        packages=find_packages(),
        description='Build multiwavelength mock cluster surveys at the drop of a hat',
        long_description=open("README.md").read(),
        package_data={"": ["README.md", "LICENSE"]},
        include_package_data=True,
        classifiers=[
            "Development Status :: 3 - Alpha",
            "Intended Audience :: Science/Research",
            "License :: OSI Approved :: MIT License",
            "Operating System :: Unix",
            "Programming Language :: Python"
            ],
        install_requires=[]
)
