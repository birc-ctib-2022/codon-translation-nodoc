# Translating protein-coding genes

*Protein-coding genes* are genes that cells first translate into RNA, then translate the RNA into amino acids, and finally use the amino acids as proteins that do whatever it is that proteins do.

We will not worry about the first step, from DNA to RNA. It is a straightforward process since each DNA nucleotide is translated to a corresponding RNA nucleotide, and we will skip it and instead translate directly from a DNA sequence to an amino acid sequence. If you can handle that, you can certainly also handle translating DNA into RNA. We will also add some simplifying assumptions to the project, and assume that a gene is a simple sequence of DNA, starting with the first code for an amino acid and ending at the last. This is close to how bacteria encode genes, but simpler than eukaryotic genes—we don’t want to go there this early in our programming career.

So, a gene, with these simplifications, is a DNA sequence where each nucleotide takes part in coding an amino acid. The way that DNA encodes amino acids is, luckily, simple. Each triplet of DNA encodes one amino acid, with the mapping from triplets from the table below:

```python
CODON_MAP = {'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
             'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S', 
             'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*', 
             'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W', 
             'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
             'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P', 
             'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q', 
             'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R', 
             'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M', 
             'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T', 
             'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K', 
             'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R', 
             'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V', 
             'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A', 
             'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E', 
             'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G'}
```

The table is close to universal for all life on Earth, with only tiny variation in some groups of species, so you rarely have to worry about other mappings.

So, if your DNA string is

```python
dna = "ATGACCGAACAGTAG"
```

consisting of triplets, 

```python
codons = ["ATG", "ACC", "GAA", "CAG", "TAG"]
```

you would expect to get the corresponding amino acids

```python
['M', 'T', 'E', 'Q', '*']
```

and you could build the following string from that:

```python
"MTEQ*"
```

## Setting up the template code

You sometimes need to install Python modules for use in your own programs. One way to do this is to use the `pip` tool distributed with Python. This tool can read the modules you need from a file and install them, including their dependencies. To do that, use the command

```sh
> python3 -m pip install -r requirements.txt
```

where `requirements.txt` contains the list of modules you want to install. There is an example `requirements.txt` in this repository that installs the `pytest` tool. Try installing it.

In the example `requirements.txt` we install a modules for testing code, `pytest`.

## Testing your code

We use two types of tests on projects in CTiB, `doctest` which is part of the standard Python library and `pytest` which we need to install. You can use them to test your code in your repository, but we also run tests every time you commit to GitHub, so we can keep track of your progress along the way.

The `doctest` module will validate if a function matches the behaviour we have described in its documentatiton string; if the documentation includes code that evaluates the function, `doctest` will test if it actually does what we claim.

You can run it on the included template code by running

```sh
> python3 -m doctest src/*.py
```

The `pytest` tool makes it easy for you to write so-called "unit tests"--small test functions that validates that different aspects of your code is working. It looks for Python files whose name starts with `test_`, and in them, it will locate all functions whose names start with `test_` and run them. There is already a test file in this repository, so you can try running `pytest` with the command

```sh
> python3 -m pytest src
```

Both tests will fail right now because there are tests to check if you have implemented the functions you need for this project, and for obvious reasons, you haven't done that yet. However, when you have implemented what you neeed to do, the tests will pass, and that will tell you that you have succeeded in your task.

Whenever you make changes to your code, you should run `doctest` and `pytest` to ensure that everything is still working. If you add new code, you can include documentation strings so `doctest` can keep you honest about the documentation. If you want more detailed testing, you can add another test, as a `test_*` function in an existing file or in a new `test_*.py` file.

When you push changes from your repository to GitHub, GitHub will also run tests on your code, and you can see the results in the `Feedback` pull request or in the `Actions` window on your repository.

## Template code

In the `src/` directory you can find two files relevant for the main task: `codons.py` and `test_codons.py`. The first is where you should write your code for this project, and the second is where you should write tests. There are already some tests, to get you started, but you may want to add more, if you wish to test your code in more detail. That is entirely up to you.

The `codons.py` file contains the table of codons to amino acid translations shown above and three functions. The first function, `split_codons()` splits a string of nucelotides into codons (triplets), the second, `translate_codons()` takes a list of these triplets and translate them into amino acids, and the last function, `translate_dna()` should translate all the way from a DNA sequence to an amino acid sequence. The last function can obviously use the first to achieve its goal, but you are free to implement it any way you want, as long as it can translate DNA strings to amino acids.

The task in this project is to implement the three functions in `codons.py`.

- [ ] Implement the `split_codons()` function according to the specification.
- [ ] Implement the `translate_codons()` function according to the specification.
- [ ] Implement the `translate_dna()` function according to the specification.

## Building a command-line tool

There is a third Python file in `src`, `src/main.py`, that will show you have to build a command-line program in Python. Generally, you can run any Python code by calling `python3` with the file that contains the code, e.g.

```sh
> python3 foo.py
```

will run any code in the file `foo.py`. However, there is a little more to writing command-line tools than this. For a tool to be useful, the user needs a mechanism to provide input to program, in some case also a way to provide optional flags, and the user needs to get the results of running the code back.

The way the user informs the program about which flags/options to use and where to find input and where to write output is through command-line arguments, and the actual data that goes into the program and comes out again goes through files (or "streams", which are essentially the same thing).

The file `src/main.py` shows you a very rudementary way of handling this in Python; in later projects we will see more advanced (and better) techniques.

You do not need to modify anything in this file, but I encourage you to read it, to get an idea about how you can turn your own code into something that works as a command-line tool. It won't be long before you will need to know how to do this.

For this project, though, once you have implemented the functions in `src/codons.py`, you can use `src/main.py` as a command-line tool. You can call it in three different ways--the source code explains how this works--depending on how you want to specify the input and where you want the output to go.

If you use the program without any arguments, it will read the input from `stdin` and write it to `stdout`, so you could use it in a pipeline as

```sh
> cat data/seqs.in | python3 src/main.py
```

or

```sh
> python3 src/main.py < data/seqs.in
```

to write the translated sequence from `data/seqs.in` to the terminal.

If you want the output to go to a file, you must redirect it

```sh
> python3 src/main.py < data/seqs.in > output.out
```

That is a very primitive interface, and a user will typically expect that a program that can read a file can also take a filename as an argument, and while not all programs accept a filename for the output, it is common enough that we should consider providing it.

You should be able to provide the input file directly to the program if you give it a single argument:

```sh
> python3 src/main.py data/seqs.in
```

and it will, again, write the translated sequences from `seqs.in` to the terminal.

If you provide two file-arguments to the tool, it should read the input data from the first file and write the output to the second file.

```sh
> python3 src/main.py data/seqs.in my-output.out
```

The arguments that a program gets are put in `sys.argv`, with the name of the program at `sys.argv[0]` and any following arguments after that. The number of real arguments this thus always `len(sys.argv) - 1`. You can test how many you have by checking `len(sys.argv)`:

```python
if len(sys.argv) == 1:
    # zero arguments -- use stdin and stdout
elif len(sys.argv) == 2:
    # one argument -- use sys.argv[1] instead of stdin
elif len(sys.argv) == 3:
    # two arguments -- use sys.argv[1] for input and sys.argv[2] for output
else:
    # more than two arguments; that is an error
```

To open an input file, a file you want to read from, use `open("filename")` or `open("filename", "r")`, and to open an output file, a file you want to write to, use `open("filename", "w")`.

 - [ ] Extend the program such that if `len(sys.argv) in [2, 3]`, the program should read from a file you open as `open(sys.argv[1])`.
 - [ ] Extend the program such that if `len(sys.argv) == 3`, the program should write to a file you open as `open(sys.argv[2], "w")`.

 If you have more than two arguments, terminate the program with `sys.exit(1)` to indicate an error.
 


We are not doing any sensible tests in the arguments in this project, so you will not get user-friendly error messages if you provide input files that do not exist, or try to write to a file you do not have permission to write to. That is something we will improve upon in later projects.

Writing `python3 program.py` to run the program `program.py` doesn't look like other programs on the command-line, where we would usually just write `program` to run `program`. There isn't anything wrong with that, it just says that we use Python to execute the program. There are, however, ways of making Python programs look more like other programs. Several, in fact. But, you guessed it, we leave that for future projects. We have already covered a lot in one project, and once you have everything up and running, you can pad yourself on the shoulders and take a short break before we continue.
