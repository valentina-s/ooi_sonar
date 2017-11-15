This work is Python 2 based. Currently we are using functionality from the [mi-instrument package](https://github.com/oceanobservatories/mi-instrument) which works only on Python 2. It is best to create a virtual environment for the project:

`conda create -n ooi_sonar2 python=2 anaconda`

To activate the environment run:

`source activate ooi_sonar2`

From within the environment install the dependencies:

`pip install -r requirements_clean.txt`

To run the notebooks run:

`cd notebooks`

`jupyter notebook`

and navigate to the notebooks of interest.

In the notebooks you will have to change the paths to where you have downloaded the dataset.

Once you are finished: 

`source deactivate`.



