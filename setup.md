This work is being converted from an old state which used lots of Python 2 from the [mi-instrument package](https://github.com/oceanobservatories/mi-instrument) to a new state which uses Python 3. It's currently undergoing major changes at the moment, so please check back in Nov 2018 for a cleaner and more functional repo.

To create a conda environment for this project:

`conda env create -f environment.yml -n py3-sonar`

To make an ipython kernel hook, use:

`python -m ipykernel install --user --name py3-sonar`

To activate the environment run:

`source activate py3-sonar`

To run the notebooks run:
```shell
$ cd notebooks
$ jupyter notebook
```
and navigate to the notebooks of interest.

In the notebooks you will have to change the paths to where you have downloaded the dataset.

Once you are finished, to return to the base environment:

`source deactivate`.
