#!/bin/bash

# script to build the ImageJ (fiji) MDA reader Plugin

# Place where ImageJ is installed

if [ -d "$1" -a -d "$1/plugins/" -a -f "$1/jars/ij.jar" ]
then
    IMAGEJ_HOME=$1

    javac -cp $IMAGEJ_HOME/jars/ij.jar Open_MDA.java
    jar cfM MDA_.jar Open_MDA.java Open_MDA.class plugins.config

    echo -e "\nCopy MDA_.jar to ImageJ Plugins folder, i.e.:\n"
    echo -e "cp MDA_.jar $IMAGEJ_HOME/plugins/\n"
else
    echo -e "\nUsage: ${0##*/} <path-to-ImageJ>\n"
    echo -e "Builds the MDA reader plugin for ImageJ/Fiji\n"
fi

