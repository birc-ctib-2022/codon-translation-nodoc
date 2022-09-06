"""
A command-line tool example.

This file shows how you can let the user of a tool provide
data for your program. The arguments that a user provides
on the command-line can always be found in the list sys.argv.
The first item in this list is the name of the program, and
if there are any additional arguments they come after that.
So, you can test how many arguments you got by taking
len(sys.argv) and subtracting one for the name of the program.

This program accepts zero, one, or two arguments, so we can check
which case we are in using

if len(sys.argv) == 1:
    # zero arguments
elif len(sys.argv) == 2:
    # one argument
elif len(sys.argv) == 3:
    # two arguments
else:
    # more than two arguments; that is an error

or we can use the new pattern matching syntax from Python 3.10
and do this, 

match len(sys.argv):
    case 1:
        # zero arguments
    case 2:
        # one argument
    case 3:
        # two arguments
    case _:
        # more than two arguments; that is an error

There is no difference between the two ways of doing it, so it
is entirely a question of taste.

The code below uses the second form to select the input and the
output file. The standard input file is found in sys.stdin and the
standard output file is found in sys.stdout. If we have a file name
we can open the file with open("filename") or open("filename", 'r')
if we want to read from it (for the input data) or with
open("filename", 'w') if we want to write to it (for the output data).

If we want to report erroneous situations, we should use sys.stderr. That way,
we do not mix error messages with the proper output that might go to sys.stdout.
If an error is so grievous that we need to inform the user about it more strongly,
we should report it as the exit status of the program. The exit status is something
the shell can check for, and the user can program around in shell-scripts. If a
program terminates normally, the exit status should be zero, and if something
went wrong, it should be non-zero. If you do nothing, Python will automatically
exit with a zero status, but you can use sys.exit(n) to exit with status n.

We are not checking if the files can be opened in this program, so errors there
are handled by Python (and it will set a proper exit status itself). We want
to report an error if we do not get the right number of arguments, though,
so there we write a message to stderr and exit with a non-zero status.

If we get past the argument parsing, we will iterate through each line in the
input, translate it, and write the result to the output. However, the translation
might fail, in which case the translation function you wrote should return None.
If that happens, we also write an error message and terminate the program with
a non-zero exit status.

"""

import sys

from codons import translate_dna

# You do not need to indent a command-line tool's code under
# a 'if __name__ == __main__:' construction, but it is common
# practise (if not exactly *best* practise). The idea behind
# it is this: when you import a Python file as a module, the
# __name__ variable is the name of the module, but if you run
# it as a program, __name__ will be '__main__'. So, this
# construction lets you specify code that should only be run
# when you use a file as a program, but not if you import
# the file for use as a module.
if __name__ == '__main__':
    infile, outfile = sys.stdin, sys.stdout
    match len(sys.argv):
        case 1:
            # zero arguments
            pass
        case 2:
            # one argument
            print("Feature not implemented yet.", file=sys.stderr)
            sys.exit(1)
        case 3:
            # two arguments
            print("Feature not implemented yet.", file=sys.stderr)
            sys.exit(1)
        case _:
            # more than two arguments; that is an error
            print("Too many arguments.", file=sys.stderr)
            sys.exit(1)

    # Process the input, one line at a time.
    for line in infile:
        # Remove the newline from `line`; we don't want to try to
        # translate that.
        line = line.strip()
        aa = translate_dna(line)
        if aa is None:
            # Something went wrong!
            print(f"Could not translate '{line}'.", file=sys.stderr)
            sys.exit(1)
        # If everthing went well, we write the result to the output
        print(aa, file=outfile)

    # It is polite to close files when we no longer need them.
    # It doesn't matter here because we are just about to terminate
    # the entire program, but it is best to get into the habit
    # so we don't forget later, where it might matter more.
    infile.close()
    outfile.close()
