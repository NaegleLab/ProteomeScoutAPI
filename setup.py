from setuptools import setup, find_packages


# Run setup
setup(
    name="proteomeScoutAPI",
    version="3.0.0",
    author="Naegle Lab",
    author_email="kmn4mj@virginia.edu",
    url="https://github.com/NaegleLab/ProteomeScoutAPI",
    install_requires=['pandas>=2', 'requests'],
    keywords=['proteomescout', 'bioinformatics', 'proteomics', 'data parsing', 'PTM'],
    license='GNU General Public License v3',
    description='ProteomeScoutAPI: Parser to interact with ProteomeScout data',
    long_description="""ProteomeScoutAPI is a Python module which can be used to interact with flatfiles from ProteomeScout, such as for annotating datasets with additional context.""",
    project_urls = {
        'Homepage': 'https://github.com/NaegleLab/ProteomeScoutAPI',
        'Documentation': 'https://naeglelab.github.io/ProteomeScoutAPI/'},
    classifiers=[
        "Programming Language :: Python :: 3",
        "Intended Audience :: Science/Research",
        "Operating System :: OS Independent",
    ],
    packages=find_packages(),
    include_package_data = True,
    python_requires=">=3.9",
    zip_safe = False
)