.. _python-quickstart:

Getting Started
===============

Here, we have provided a quick start guide that will allow you to get the KSTAR python package up and running quickly with the default settings.

Installation
------------

The proteomeScoutAPI can be installed via ``pip``, tarball, and directly from the Git repository.

Pip
~~~
To install via pip, execute::

    pip install proteomeScoutAPI


Tarball
~~~~~~~
To install via a tarball, head over to the `Releases page <https://github.com/NaegleLab/proteomeScoutAPI/releases>`_ and download the latest stable tar release.

Afterwards, navigate to your downloads directory and execute the following commands, substituting <version> for the release's version number::

    tar -xvf proteomeScoutAPI-<version>.tar.gz
    cd proteomeScoutAPI-<version>
    python setup.py install

Git
~~~
If you want to try out the latest commit, you can install directly from the Git repository by executing the following commands::

    git clone https://github.com/NaegleLab/proteomeScoutAPI
    cd proteomeScoutAPI
    python setup.py install