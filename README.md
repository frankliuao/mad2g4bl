# mad2g4bl

**Security.** This tool has been unmaintained since 2015 and is Python 2 only.
This PR removes the remote `exec` in `getPDGdata`, `exec`/`eval` of MAD-X
sequence and TWISS input, and the shell `rm` of the output file. Other known
bugs (see review) remain; this is not a Python 3 port and it does not include
the other functional fixes.

Convert a MAD-X sequence and TWISS file to a G4beamline `.g4bl` file.
Copyright Ao Liu (frankliuao), 2015.
