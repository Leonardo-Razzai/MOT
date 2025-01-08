# AnMOT.py

This script is designed to handle and analyze data of photoemission of atoms loaded in a Magneto-Optical Trap (MOT). The main functionalities of the `MOTdata` class include:

- Importing data from CSV files.
- Selecting specific data ranges.
- Fitting Lorentzian functions to the data.
- Calculating the number of atoms.
- Plotting raw and processed data.
- Saving processed data and fit parameters to CSV files.

## Requirements

- Python 3.x
- numpy
- matplotlib
- pandas
- pathlib
- scipy

## Usage

### Importing the Script

To use the script, import it into your Python environment:

```python
from AnMOT import MOTdata
```

### Initializing the MOTdata Class

Create an instance of the `MOTdata` class by providing the date, magnetic field value, list of detuning values, and detuning limit:

```python
date = "2023-10-01"
B = 1.0 # magnetic field in A
BN_vals = [1000, 1020, 1040, 1060] # Cooler beat-notes for photo
D_L = 1030 # Cooler beat-note for mot loading

mot_data = MOTdata(date, B, Dvals, D_L)
```

### Importing Data

Import data from CSV files located in a directory named after the date:

```python
mot_data.ImportData()
```

### Selecting Data

Select data within a specified time range and optionally plot the data:

```python
mot_data.SelectData(tmin=3.8, tmax=3.9, plot=True)
```

### Fitting a Lorentzian Function

Fit a Lorentzian function to the selected data and optionally plot the fit:

```python
mot_data.FitLorentzian(plot=True)
```

### Calculating the Number of Atoms

Calculate the number of atoms based on the fitted Lorentzian parameters:

```python
mot_data.CalculateNumOfAtoms()
```

### Plotting Raw Data

Plot the raw data:

```python
mot_data.PlotRawData()
```

Plot the raw data within a specified time range:

```python
mot_data.PlotRawData_zoom(tmin=3.8, tmax=3.9)
```

### Saving Data

Save the photoemission data to a CSV file:

```python
mot_data.SaveData()
```

Save the Lorentzian fit parameters to a CSV file:

```python
mot_data.SaveParameters()
```

## Functions and Methods

### Lorentzian_in_MOT

```python
def Lorentzian_in_MOT(BN, V0, BN0, Gamma):
    """
    Lorentzian function for MOT.
    """
```

### Prob_photoemission

```python
def Prob_photoemission(Delta, s0, Gamma):
    """
    Photoemission probability for a single atom at detuning Delta.
    """
```

### NumOfAtoms

```python
def NumOfAtoms(V0, dV0, Delta, Gamma):
    """
    Calculate the number of atoms.
    """
```

### MOTdata Class

#### __init__

```python
def __init__(self, Date: str, B: float, BN_vals: list, D_L: float):
    """
    Initialize the MOTdata class.
    """
```

#### ImportData

```python
def ImportData(self):
    """
    Import data from CSV files.
    """
```

#### SelectData

```python
def SelectData(self, tmin=3.8, tmax=3.9, plot=False):
    """
    Select data within a specified time range.
    """
```

#### FitLorentzian

```python
def FitLorentzian(self, plot=False):
    """
    Fit a Lorentzian function to the data.
    """
```

#### CalculateNumOfAtoms

```python
def CalculateNumOfAtoms(self):
    """
    Calculate the number of atoms.
    """
```

#### PlotRawData

```python
def PlotRawData(self):
    """
    Plot the raw data.
    """
```

#### PlotRawData_zoom

```python
def PlotRawData_zoom(self, tmin=3.8, tmax=3.9):
    """
    Plot the raw data within a specified time range.
    """
```

#### SaveData

```python
def SaveData(self):
    """
    Save the photoemission data to a CSV file.
    """
```

#### SaveParameters

```python
def SaveParameters(self):
    """
    Save the Lorentzian fit parameters to a CSV file.
    """
```
