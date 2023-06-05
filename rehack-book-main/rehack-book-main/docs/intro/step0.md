# Before we start: How to follow this tutorial

This tutorial contains a mix of lecture-based and interactive components.
The interactive components can be executed in your personal computer locally or using the [Binder](https://jupyter.org/binder) service.
You are welcome to follow along however you like, whether you just want to listen or code along with us.

```{attention}
Regardless of which setup method you choose, all of the Jupyter notebooks can be found in the `docs/notebook` folder.
```

## <i class="fa fa-rocket" aria-hidden="true"></i> Using Binder

Clicking the Binder button below will launch an interactive computing environment with all of the necessary software packages pre-installed.

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/nipreps/nipreps-book/main?urlpath=lab)

This is the easiest and quickest way to get started.

```{attention}
If using Binder, please be aware that the state of the computational environment it provides is not permanent.
If you are inactive for more than 10 minutes, the environment will timeout and all data will be lost.
If you would like to preserve any progress, please save the changed files to your computer.
```

## <i class="fas fa-hammer"></i> Local installation ("bare-metal")

If you would like to follow along using your own setup and you have a functional Python environment, you can run the following commands in your terminal:

```bash

# 1. clone this repository
git clone https://github.com/dandi/rehack-book

# 2. install the necessary python packages in your Python environment
cd rehack-book && pip install -r requirements.txt

# 3. launch a Jupyter lab instance
jupyter lab

```
