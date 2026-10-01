Recurrent two-dimensional sequences over three and four letters
simulating von Neumann constructions of finite ordinals

Verification code and data accompanying the paper.

Contents
--------
code/  Recurrence definitions, verification programs, and figure generators.
data/  Finite tables, verification logs, and the origin-wall certificate.

Requirements
------------
Python 3 and GCC with OpenMP and undefined-behaviour sanitizer support.
The Python programs use the standard library. The C programs expect a
64-bit long, as on 64-bit Linux.

Verification checks
-------------------
Run these commands from this folder:

    python3 code/check_three_initial_local.py
    gcc -O2 -fsanitize=undefined -fno-sanitize-recover=undefined code/check_window_maps.c -o /tmp/check_window_maps
    /tmp/check_window_maps
    gcc -O2 -fopenmp -fsanitize=undefined -fno-sanitize-recover=undefined code/check_search_detector.c -o /tmp/check_search_detector
    /tmp/check_search_detector

These check the local catalogue for three initial symbols, the window maps,
and the search detector, including the four-letter positive control.
The recorded results are in:

    data/three_initial_local_verification.txt
    data/window_map_verification.txt
    data/search_detector_verification.txt

Exhaustive three-letter search
-----------------------------
    gcc -O3 -fopenmp code/search_three_letters.c -o /tmp/search_three_letters
    OMP_NUM_THREADS=8 /tmp/search_three_letters

This enumerates all 129,140,163 normalized symmetric rules. The recorded
completed run is in data/three_letter_search.txt.

Other verification programs
---------------------------
    python3 code/three_initial.py --size 2048 --rows 60000
    python3 code/check_three_initial_structure.py
    python3 code/check_counts_and_boundaries.py

For the two-initial-symbol sequence, verify.c and desc.h implement the grid
and local recurrence checks. coordcheck.py, tsquares.py, eclass.py,
envcheck.py, zones.py, and triples.py check geometry and argument triples.

The four-letter checks are in check_four_letters.c, check_virtual_boxes.c,
triple.c, desc4loc.c, vbox.c, whites.c, line1.c, and band.c. Build commands
and parameters are documented in their sources.

Verification logs and finite tables are in data/. The largest grid and
band checks can require substantial time and memory.
