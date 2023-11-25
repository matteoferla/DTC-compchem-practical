# Notebook basics
A notebook is composed of 'cells'.
This is a markdown cell, it is for words, the next is for code.
The code cells can be run by pressing the play button at the top or shift+Enter.
You can freely change things as you will have to in the next cell.

You have several varieties of notebooks:

* (original=vanilla) Jupyter notebook, runs on your machine
* Jupyter lab, newer version of the above
* Colab, run on Google cloud
* Jupyter hub, the multiuser version of Jupyter

Not sure how to use a Python object? Make a cell a type `help(🤖)` where 🤖 is the object.
(No, you cannot use emoji as variables in Python).
In this notebook, the imported module `dtc` will have `dtc.show_source(🤖)` will show the source code in a colourful way
(it calls `inspect.getsource` a handy function).

The output of the last command will be passed to `IPython.display.display`, which renders them based on the magic method
`_ipython_display_` or `_repr_html_` etc. and _in extremis_ `__repr__` then `__str__`. Hence why some cells end in `None` herein.


