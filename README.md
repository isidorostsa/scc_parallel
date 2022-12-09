# pdsProject

Project 1 for Parallel and Distributed Systems in A.U.Th.

The report is the "report.pdf" file contained here.

Build any of the files by going to the appropriate directory and running `make`.

To build the OpenCilk version edit the makefile to point to the location of clang++ in your system.

- A possible problem is clang++ being unable to locate the c++ standard library. This is solved by: [Clang++: fatal error: 'iostream' file not found - Stack Overflow](https://stackoverflow.com/questions/54521402/locating-iostream-in-clang-fatal-error-iostream-file-not-found)

stdc++20 is required, but this is trivial to change. If you need this to change please contact me at
isidtsao@ece.auth.gr

Run with the command `./colorSCC` to see the usage information.

Have fun testing!
