Installation
============

.. _installation:

Prerequisit packages
--------------------

Since this code uses the MELTS code as provided by the `Thermoengine <https://enki-portal.gitlab.io/ThermoEngine/>`_ 
package, this needs to be installed before LavAtmos may be used. We 
suggest to follow the instructions on the `Thermoengine GitLab page <https://gitlab.com/ENKI-portal/ThermoEngine>`_ 
on how to run a container image locally. Note: This involves installing 
`Docker <https://www.docker.com/get-started/>`_ onto your system. 

.. _target to cloning_repo:
Cloning the repository
----------------------

Once you have installed Thermoengine, clone the `LavAtmos repository <https://github.com/cvbuchem/LavAtmos>`_ into the Thermoengine directory using::

    git clone https://github.com/cvbuchem/LavAtmos

If you do not have GitHub installed you can also download a zip file containing the code by going to the `Github page <https://github.com/cvbuchem/LavAtmos>`_, clicking on the green "Code" button and then clicking "Download ZIP".

Running LavAtmos using Singularity
----------------------------------

If you are for some reason unable to install docker (you may not have root
access to your system for example) it is also possible to run the Thermoengine
docker image using `Singularity <https://docs.sylabs.io/guides/2.6/user-guide/index.html>`_.
You can check if the application is already installed by typing::

    singularity --help

If not, follow the instructions on their `Read the Docs page <https://docs.sylabs.io/guides/2.6/user-guide/installation.html>`_.

You can then pull the Thermoengine docker as follows::

    singularity pull docker://registry.gitlab.com/enki-portal/thermoengine:master

.. note:: 
    
    If you get the following error::

        FATAL:   While making image from oci registry: error fetching image to cache: while building SIF from layers: unable to create new build: while searching for mksquashfs: exec: "mksquashfs": executable file not found in $PATH

    A potential fix is to enter the following bash commands into your terminal::

        export PATH="/usr/sbin:$PATH"

This will create a ``.sif`` file which allow you to shell into the Thermoengine docker image (see the `Singularity documentation <https://docs.sylabs.io/guides/2.6/user-guide/singularity_and_docker.html>`_ for more information). 
In order to run ``LavAtmos``, clone the repository using the instrunctions given :ref:`above <target to cloning_repo>`:. You can then use the following commands to run an example script::

    singularity shell thermoengine_master.sif 
    cd '/path/to/LavAtmos'
    python3 scripts/run_lavatmos_example2.py

You may then alter the example script to suit your needs. 
