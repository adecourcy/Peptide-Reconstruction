# Purpose

This is a program written to find process mass spectrometry data given MHC-I peptides. Essentially, it is a collection of heurstics to do a fast, efficient, thorough search of the mass spec data and figure out which peptide generated the spectrum. The basics of mass spectrometry are honestly difficult to explain, but I've attempted to give a brief introduction below.

This program is separated into 3 main parts: data pre-processing, data post-processing, and the search algorithm itself. 

The "main.py" file drives the program, taking the user input and executing each step as necessary.

The "backend" folder contains all pre-processing and post-processing functions, as well as functions for the data structures used by the program (each section separated into a separate subfolder).

The "scoring" folder contains the code for the main algorithm, written in C. It was translated from the original Python code in pretty short order, and it honestly has a lot of room for improvement.

Finally, the "analysis" folder contains a set of tools for analyzing the output of the program and evaluating its efficacy.


## Very Basics of Mass Spectrometry

The goal of this program is to discover new peptides. In this case, a peptide is a short protien, made up of a linear chain of amino acids. We are generally dealing with 20 amino acids. Each can be represented by a single letter and a mass:

[[https://github.com/adecourc/adecourcy/Peptide-Reconstruction/images/amino acids.png|alt=octocat]]
