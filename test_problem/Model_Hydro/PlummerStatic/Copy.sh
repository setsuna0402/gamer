#!/bin/bash 

# copy files to the correct directories for the target test problem
cp Init_TestProb.cpp                   ../../../src/Init/
cp Flu_BoundaryCondition_User.cpp      ../../../src/Fluid/
cp End_*                               ../../../src/Init/
cp Makefile                            ../../../src/
cp Input__*                            ../../../bin/Run/
cp Plot__*                             ../../../bin/Run/
