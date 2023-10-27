#!/bin/bash
#clear
if [ "$#" -ne 1 ]
then
  echo "Provide a name"
  exit 1
fi

echo "namespace gismo {
/**

\page $1 $1.cpp


<!-- more info on the file examples/$1.cpp go here -->


\section $1Annotated Annotated source file

Here is the full file \c examples/$1.cpp. Clicking on a function
or class name will lead you to its reference documentation.

\include $1.cpp

*/

}" > $1.dox
