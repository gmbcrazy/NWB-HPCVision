# NeurodataReHack Event Book

Creating a Jupyter book for projects from the NeurodataRehack event

## How to contribute


1. Fork+Git clone the repo.
2. Add a new folder inside docs
3. Convert your notebooks to myst+markdown using jupytext:
   `jupytext --to md:myst mynotebook.ipynb -o mynotebook.md`
4. Add the markdown file to the new folder
5. Add a pointer to the folder in the table of contents (toc.yml)
6. Send a PR
