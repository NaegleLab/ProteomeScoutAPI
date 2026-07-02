from setuptools import setup, find_packages

# read the contents of your README file
from pathlib import Path
this_directory = Path(__file__).parent
readme_path = this_directory / "README.md"
if not readme_path.exists():
    readme_path = this_directory / "readme.md"
long_description = readme_path.read_text(encoding="utf-8")

# Run setup
setup(
    name="proteomeScoutAPI",
    version="3.1.1",
    author="Naegle Lab",
    author_email="kmn4mj@virginia.edu",
    url="https://github.com/NaegleLab/ProteomeScoutAPI",
    install_requires=['pandas>=2', 'requests'],
    keywords=['proteomescout', 'bioinformatics', 'proteomics', 'phosphoproteomics', 'data parsing', 'PTM'],
    license='Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International',
    description='ProteomeScoutAPI: Parser to interact with ProteomeScout data',
    long_description=long_description,
    long_description_content_type="text/markdown",
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