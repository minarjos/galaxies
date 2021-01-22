# Spiral Galaxies Solver

## Spiral Galaxies

Spiral galaxies is logical puzzle, you are supposed to partition given grid to _galaxies_. You are given the centers of galaxies, every galaxy has to contain exactly one center and has to be conected. Furthermore, the galaxies must be symetric by their center.

## Instalation

Requires ncurses, simply run `make galaxies`.

## Usage 

By default, solver reads problem from `stdin` and prints solution to `stdout`. The solver expects input in specific format - it is quite simple, can be easily deduced from several included examples. You can try `./galaxies < small`.

Option `-l <file>` loads puzzle from specified file (it is useful for interactive mode because you can't redirect input in that case).

There also is option `-i` which toggles interactive mode. This mode allows to edit the grid. There are several commands.

- _Arrow Keys_ - move cursor in given direction
- _Space_ - add or remove galaxy center
- _Enter_ - solve current board
- _Escape_ - exit the solver
- _Backspace_ - delete all centers
- _r/R_ - remove/add row
- _c/C_ - remove/add column
- _l_ - load board from file
- _s_ - save board to file

To start the interactive mode, run `./galaxies -i`.

## Performance 

Solver performes very well on standard puzzles (puzzles which are supposed to be solvable by humans and have one unique solution). There are examples of such puzzles included - `small`, `medium` and `big`, solver is able to solve all of them in a few miliseconds. It also performes well on quite difficult puzzles (at least difficult for humans), again, there is an example `hard`. Solution contains very big and conmplicated galaxies.

However, solver struggles on boards with few galaxies. Even though such boards have many solutions and humans can see them imedialtelly. Example of such board is `slow`, solver needs severeal seconds to solve it.
