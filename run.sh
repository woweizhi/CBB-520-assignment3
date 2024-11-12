#!/usr/bin/env python

echo "****************** Run the script of Group 6! ******************"
python src/protein1D_pattern_group6.py

echo "****************** Run the script of Group 7! ******************"
python src/protein1D_pattern_group7.py

echo "****************** Prepare the stride output for Group 8&9! ******************"
bash stride.sh

echo "****************** Run the script of Group 8! ******************"
python src/protein2D_pattern_group8.py

echo "****************** Run the script of Group 9! ******************"
python src/protein2D_pattern_group9.py

echo "****************** Run the script of Group 10! ******************"
python src/protein3D_pattern_group10.py

echo "****************** Finish running all the code! ******************"