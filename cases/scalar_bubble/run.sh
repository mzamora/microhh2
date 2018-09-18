#!/bin/bash

export PYTHON_EXEC=python
export MICROHH_EXEC=$(cd ../../build && pwd)/microhh

echo "Cleaning folder"
rm -f *.000*
rm -f *.out

echo "initializing scalar bubble w/zero scalars"
$MICROHH_EXEC init scalar_bubble> log
echo "Running scalar bubble w/zero scalars"
$MICROHH_EXEC run scalar_bubble>>log
echo "Introducing the scalar bubble and new .ini"
cp scalar_bubble.ini scalar_bubble_restart.ini
cp scalar_bubble.prog scalar_bubble_restart.prof
$PYTHON_EXEC make_scalar_bubble.py
echo "Running from the scalar bubble onwards"
$MICROHH_EXEC run scalar_bubble_restart>>log
echo "Finding position of bubble at t=20"
$PYTHON_EXEC find_bubble.py



