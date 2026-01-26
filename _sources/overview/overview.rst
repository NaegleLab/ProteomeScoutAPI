Getting Started
=================


A lightweight API to talk to ProteomeScout flatfiles in Python. 

About
-----
Version 2.0 - June 2020

Originally: By Alex Holehouse, Washington University in St. Louis Contact alex.holehouse@gmail.com or contribute at [https://github.com/alexholehouse](https://github.com/alexholehouse)

Current Version: Naegle Lab, University of Virginia [https://github.com/naegleLab](https://github.com/naegleLab)

Overview
--------
ProteomeScoutAPI is a Python module which can be used to connect to and parse ProteomeScout flatfiles. Specifically, the goal of this module is to allow anyone to interact with ProteomeScout data without the need to

1. Repeatedly interact with ProteomeScout web interface

2. Have any knowledge of SQL, or use an SQL-Python ORM

3. Facilitate rapid exploration of the ProteomeScout dataset

Installation
------------

ProteomeScoutAPI can be installed via ``pip``, tarball, and directly from the Git repository. 

Pip
~~~
To install via pip, execute::

    pip install ProteomeScoutAPI


Tarball
~~~~~~~
To install via a tarball, head over to the `Releases page <https://github.com/NaegleLab/KSTAR/releases>`_ and download the latest stable tar release.

Afterwards, navigate to your downloads directory and execute the following commands, substituting <version> for the release's version number::

    tar -xvf ProteomeScoutAPI-<version>.tar.gz
    cd ProteomeScoutAPI-<version>
    python setup.py install

Git
~~~
If you want to try out the latest commit, you can install directly from the Git repository by executing the following commands::

    git clone https://github.com/NaegleLab/ProteomeScoutAPI
    cd ProteomeScoutAPI
    python setup.py install


Configuration
-----------------

By default, the API will automatically look for and download the most recent version of the ProteomeScout flatfiles into the package directory. For most users, this behavior is sufficient and there is no need to change anything.

However, if you would like to adjust this behavior, you can update the default configuration:


.. code-block:: python

    from proteomeScoutAPI import update_configuration

    # where to save ProteomeScout flatfiles
    pscout_data_path = "/path/to/your/directory"

    # default dataset version to use (if None, will always use latest)
    dataset_version = None

    # whether to automatically update ProteomeScout flatfiles if new version is available
    update = True

    update_configuration(dataset_dir=pscout_data_path,
                         version=dataset_version,
                         update=update)


Contributing code
-----------------
The code here is incredibly simple, and the few methods presented give a good example of how one should parse the ProteomeScout records. If you're interested in adding the ability to parse out other information please go ahead and make a pull request. Tests would be appreciated too!
