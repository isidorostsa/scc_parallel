# Strongly Connected Components using different parallelism backends

A Project for Parallel and Distributed Systems in A.U.Th.

There is a report of benchmarks and findings named "report.pdf" contained here.

Build any of the files by going to the appropriate directory and running `make`.

To build the OpenCilk version edit the makefile to point to the location of clang++ in your system.

- A possible problem is clang++ being unable to locate the c++ standard library. This is solved by: [Clang++: fatal error: 'iostream' file not found - Stack Overflow](https://stackoverflow.com/questions/54521402/locating-iostream-in-clang-fatal-error-iostream-file-not-found)

stdc++20 is required, but this is trivial to change. If you need this to change please contact me at
isidtsao@ece.auth.gr

Run with the command `./colorSCC` to see the usage information.

The report uploaded here has updated and more accurate measurements. Some benchmarking was done after the deadline due to no longer needing the computer to write the report.

Key diagrams can also be found here: https://imgur.com/a/KI79XMh

System Used: Lenovo Ideapad 5 pro, AMD 5600h 6 Core - 12 Thread CPU, 16GB 3200MHZ DDR4 RAM

Have fun testing!
