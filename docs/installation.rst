.. highlight:: shell

============
Installation
============


Stable release
--------------

To install WaveBreaking, run this command in your terminal:

.. code-block:: console

    $ pip install wavebreaking ### NOT AVAILABLE YET

This is the preferred method to install WaveBreaking, as it will always install the most recent stable release. Your virtual environment is automatically checked for the necessary dependencies. After the installation, you can start calculating wave breaking events by following the `tutorial`_. 

.. _tutorial: https://

From sources
------------

The sources for WaveBreaking can be downloaded in two different ways. You can either install WaveBreaking directly from the GitHub repository:

..  code-block:: 

        pip install git+https://github.com/skaderli/WaveBreaking

Or you can clone the GitHub repository first and then install WaveBreaking locally. First, set the working directory and clone the repository. 

..  code-block:: 

        cd /path/to/local/workspace
        git clone https://github.com/skaderli/WaveBreaking.git

Second, set up the conda environment and install the necessary dependencies:

..  code-block:: 

        conda create -y -n wb_env
        conda env update -f environment.yml -n wb_env

Now the environment can be activated and the WaveBreaking package can be locally installed by using the developer mode "-e":

.. code-block::

        conda activate wb_env
        pip install -e .

To check if the installation was successful, perform some tests:

.. code-block::
 
        python -m pytest
