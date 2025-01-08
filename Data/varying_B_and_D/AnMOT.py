import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pathlib import Path
from scipy.optimize import curve_fit

pm = "\u00B1"
MEDIUM_SIZE = 11
BIGGER_SIZE = 13

plt.rc('font', size=MEDIUM_SIZE)          # controls default text sizes
plt.rc('axes', titlesize=BIGGER_SIZE)     # fontsize of the axes title
plt.rc('axes', labelsize=MEDIUM_SIZE)    # fontsize of the x and y labels
plt.rc('xtick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('ytick', labelsize=MEDIUM_SIZE)    # fontsize of the tick labels
plt.rc('legend', fontsize=MEDIUM_SIZE)    # legend fontsize

base_font = {'family': 'serif',
        'size': MEDIUM_SIZE,
        }

title_font = {
  'fontsize' : BIGGER_SIZE, 
  'font' : 'serif', 
  'weight' : 'bold'
}

def Lorentzian_in_MOT(BN, V0, BN0, Gamma):
    """
    Lorentzian function for MOT.

    Parameters
    ----------
    BN : float
        Cooler beat-note frequency in MHz.
    V0 : float
        Amplitude of the Lorentzian.
    BN0 : float
        Cooler beat-note frequency at resonance in MHz.
    Gamma : float
        FWHM of the Lorentzian in MHz.

    Returns
    -------
    float
        Value of the Lorentzian function.
    """
    s0 = 1.28
    return V0 * s0/2 / (1 + s0 + 4*(BN-BN0)**2 / Gamma**2)

def Prob_photoemission(Delta, s0, Gamma):
    """
    Photoemission probability for a single atom at detuning Delta.

    Parameters
    ----------
    Delta : float
        Detuning in MHz.
    s0 : float
        Saturation parameter.
    Gamma : float
        FWHM of the Lorentzian in MHz.

    Returns
    -------
    float
        Photoemission probability.
    """
    return s0/2 / (1 + s0 + 4 * Delta**2 / Gamma**2)                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                            

def NumOfAtoms(V0, dV0, Delta, Gamma):
    """
    Calculate the number of atoms.

    Parameters
    ----------
    V0 : float
        Signal voltage.
    dV0 : float
        Error in signal voltage.
    Delta : float
        Detuning in MHz.
    Gamma : float
        FWHM of the Lorentzian in MHz.

    Returns
    -------
    tuple
        Number of atoms and its error.
    """
    s0 = 1.28
    Decay_rate = 38.11 * 1e6 # Hz
    Ep = 1.589 * 1.602 * 1e-19 # J
    
    tau = 2.38 *1e6 # V/A 5% error
    eta = 0.5 # A/W 10 % error
    R_lens = 2 # cm
    dist_mot = 40 # cm
    sigma = np.pi * (R_lens/dist_mot)**2
    
    I = V0 / tau # A
    P = I / eta # W
    
    # relative errors
    dV0_rel = dV0/V0
    dtau_rel = 0.05
    deta_rel = 0.1
    ddist_mot_rel = 2 / 40
    dsigma_rel = 2 * ddist_mot_rel
    
    Rs = Prob_photoemission(Delta, s0, Gamma) * Decay_rate # Hz
    Na = P / (sigma * Rs * Ep)
    dNa = (dV0_rel + dtau_rel + deta_rel + dsigma_rel) * Na
    
    return (Na, dNa)

class MOTdata():
    """
    Class to handle MOT data.

    Attributes
    ----------
    Date : str
        Date of the data.
    B : float
        Magnetic field value.
    BN_vals : list
        List of detuning values.
    D_L : float
        Detuning limit.
    color_palette : list
        Color palette for plotting.
    DATA_RAW : list
        Raw data.
    ERRORS : list
        Errors in the data.
    PHOTOEMISSION_VALS : list
        Photoemission values.
    PHOTOEMISSION_TIMES : list
        Photoemission times.
    V0 : float
        Amplitude of the Lorentzian.
    BN0 : float
        Cooler beat-note frequency at resonance.
    Gamma : float
        FWHM of the Lorentzian.
    dV0 : float
        Error in amplitude of the Lorentzian.
    dBN0 : float
        Error in cooler beat-note frequency at resonance.
    dGamma : float
        Error in FWHM of the Lorentzian.
    """
    def __init__(self, Date: str, B: float, BN_vals: list, D_L: float):
        """
        Initialize the MOTdata class.

        Parameters
        ----------
        Date : str
            Date of the data.
        B : float
            Magnetic field value.
        BN_vals : list
            List of detuning values.
        D_L : float
            Detuning limit.
        """
        self.Date = Date
        self.B = B
        self.BN_vals = BN_vals
        self.D_L = D_L
        self.color_palette = color_palette = plt.cm.jet(np.linspace(0.6, 1, len(self.BN_vals)))
        
        self.DATA_RAW = []
        self.ERRORS = []
        self.PHOTOEMISSION_VALS = []
        self.PHOTOEMISSION_TIMES = []
        
        self.V0, self.BN0, self.Gamma = (0, 0, 0)
        self.dV0, self.dBN0, self.dGamma = (0, 0, 0)
        
    def ImportData(self):
        """
        Import data from CSV files.
        """
        directory = Path(f'./{self.Date}/')
        
        self.DATA_RAW = []
        self.ERRORS = []
        
        for i, D in enumerate(self.BN_vals):
            file_list = list(directory.glob(f"B={self.B}_D={D}.csv"))
            df = pd.read_csv(file_list[0], header=1, names=['Time [s]', 'Sig [V]'])
            
            x_data = df['Time [s]'].to_numpy()
            x_data = x_data - np.min(x_data)
            index_offset = (x_data > 0.2) * (x_data < 0.5)
            y_data = df['Sig [V]'].to_numpy()
            offset = np.mean(y_data[index_offset])
            
            y_data = y_data - offset
            dy_data = np.std(y_data[index_offset])
            
            self.DATA_RAW.append((x_data, y_data))
            self.ERRORS.append(dy_data)
    
    def SelectData(self, tmin = 3.8, tmax = 3.9, plot = False):
        """
        Select data within a specified time range.

        Parameters
        ----------
        tmin : float, optional
            Minimum time value, by default 3.8.
        tmax : float, optional
            Maximum time value, by default 3.9.
        plot : bool, optional
            Whether to plot the data, by default False.
        """
        self.PHOTOEMISSION_VALS = []
        self.PHOTOEMISSION_TIMES = []
        
        if plot:
            _, ax = plt.subplots(1, 2, figsize=(12, 5))

            ax[0].set_title(f'Derivative of MOT fluorescence', fontdict=title_font)
            ax[0].set_xlabel('Time [ms]', fontdict=base_font)
            ax[0].set_ylabel('Sig [V]', fontdict=base_font)
            ax[0].grid()
            
            ax[1].set_title(f'MOT fluorescence: zoom about peaks', fontdict=title_font)
            ax[1].set_xlabel('Time [ms]', fontdict=base_font)
            ax[1].set_ylabel('Sig [V]', fontdict=base_font)
            ax[1].grid()

        for i, D in enumerate(self.BN_vals):
            
            x_data = self.DATA_RAW[i][0]
            y_data = self.DATA_RAW[i][1]
            
            index_diff = (x_data > tmin) * (x_data < tmax)
            x_diff = x_data[index_diff][:-1]
            y_diff = np.diff(y_data[index_diff])
            y_data_diff = y_data[index_diff][:-1]
            
            if plot:
                ax[0].plot(x_diff * 1e3, y_diff, label=f'D = {D} MHz', color=self.color_palette[i])  
            
            delta_time = 0.00015
            
            if D < self.D_L:
                # select data where diff is max
                index_photo = np.argmax(y_diff)
            else:
                # select data where diff is min
                index_photo = np.argmin(y_diff)

            x_photo = x_diff[index_photo]
            y_photo = y_diff[index_photo]
            
            index_interval = (x_data > x_photo + delta_time) * (x_data < x_photo + 2*delta_time)
            x_interval = x_data[index_interval] - x_photo
            y_interval = y_data[index_interval]
            
            x_picked = np.mean(x_interval)
            y_picked = np.mean(y_interval)
            
            self.PHOTOEMISSION_TIMES.append(x_picked)
            self.PHOTOEMISSION_VALS.append(y_picked)
            self.ERRORS[i] = self.ERRORS[i] + np.std(y_interval)
            
            if plot:
                
                index_plot = (x_data > x_photo - 50*delta_time) * (x_data < x_photo + 50*delta_time)
                x_plot = x_data[index_plot] - x_photo
                y_plot = y_data[index_plot]
            
                ax[1].plot(x_plot * 1e3, y_plot, color=self.color_palette[i])
                ax[1].plot([x_picked * 1e3], [y_picked], '+', color='royalblue')

        if plot:
            ax[1].plot([x_picked * 1e3], [y_picked], '+', color='royalblue', label='Peaks')
            #plt.legend(prop=base_font)
            plt.show()
    
    def FitLorentzian(self, plot = False):
        """
        Fit a Lorentzian function to the data.

        Parameters
        ----------
        plot : bool, optional
            Whether to plot the fit, by default False.
        """
        x_data_fit = np.array(self.BN_vals)
        y_data_fit = np.array(self.PHOTOEMISSION_VALS)
        
        popt, pcov = curve_fit(Lorentzian_in_MOT, x_data_fit, y_data_fit, p0=[0.3, 1027, 1], sigma=np.array(self.ERRORS))
        
        self.V0, self.BN0, self.Gamma = popt
        self.dV0, self.dBN0, self.dGamma = np.diag(pcov)

        x_fit = np.linspace(990, 1060, 200)
        y_fit = Lorentzian_in_MOT(x_fit, *popt)
        
        if plot:
            _, ax = plt.subplots(1, 1, figsize=(10, 5))
            
            V0_val = r'$V_0$ = ' + f'({self.V0:.4f} {pm} {self.dV0:.4f}) V\n'
            BN0_val = r'$BN_0$ = ' + f'({self.BN0:.1f} {pm} {self.dBN0:.1f}) MHz\n'
            Gamma_val = r'$\Gamma$ = ' + f'({self.Gamma:.1f} {pm} {self.dGamma:.1f}) MHz\n'
            Label = 'Fit params:\n'+V0_val + BN0_val + Gamma_val
            
            ax.plot(x_fit, y_fit, '--', color='crimson', label=Label)
            ax.errorbar(x_data_fit, y_data_fit, yerr=self.ERRORS, fmt='o', ms=4, color='royalblue', label='\nData')
            
            ax.set_title(f'Fit of fluorescence signal vs Cooler Beat-Note ', fontdict=title_font)
            ax.set_xlabel(r'B-N Cooler [MHz]', fontdict=base_font)
            ax.set_ylabel('Sig [V]', fontdict=base_font)
            ax.grid()

            plt.legend()
            plt.show()

    def CalculateNumOfAtoms(self):
        """
        Calculate the number of atoms.
        
        For each detuning value, the number of atoms is calculated, then
        weighted average and error associated to these values are calculated.
        
        """
        BN = np.array(self.BN_vals)
        Delta = BN - self.BN0
        Na, dNa = NumOfAtoms(self.V0, self.dV0, Delta, self.Gamma)
        
        Na_mean, Na_std = np.average(Na, weights=1/dNa**2), np.sum(1/dNa**2)**-0.5
          
        print(f'Number of atoms = ({Na_mean/1e7:.2e} {pm} {Na_std/1e7:.2e}) * 10^7')
    
    def PlotRawData(self):
        """
        Plot the raw data.
        """
        _, ax = plt.subplots(1, 1, figsize=(12, 5))
        ax.set_title(f'MOT fluorescence signal B = {self.B} A', fontdict=title_font)

        for i, D in enumerate(self.BN_vals):
          
            x_data = self.DATA_RAW[i][0]
            y_data = self.DATA_RAW[i][1]
            
            ax.plot(x_data, y_data, color=self.color_palette[i])
            ax.set_xlabel('Time [s]', fontdict=base_font)
            ax.set_ylabel('Sig [V]', fontdict=base_font)

        plt.grid()
        plt.show()
    
    def PlotRawData_zoom(self, tmin = 3.8, tmax = 3.9):
        """
        Plot the raw data within a specified time range.

        Parameters
        ----------
        tmin : float, optional
            Minimum time value, by default 3.8.
        tmax : float, optional
            Maximum time value, by default 3.9.
        """
        _, ax = plt.subplots(1, 1, figsize=(12, 5))
        ax.set_title(f'MOT fluorescence at different detunings : zoom', fontdict=title_font)

        for i, D in enumerate(self.BN_vals):
            
            x_data = self.DATA_RAW[i][0]
            y_data = self.DATA_RAW[i][1]
            
            index_to_plot = (x_data > tmin) * (x_data < tmax)
            ax.plot(x_data[index_to_plot], y_data[index_to_plot], color=self.color_palette[i])
            ax.set_xlabel('Time [s]', fontdict=base_font)
            ax.set_ylabel('Sig [V]', fontdict=base_font)

        plt.grid()
        plt.show()
    
    def SaveData(self):
        """
        Save the photoemission data to a CSV file.
        """
        if len(self.PHOTOEMISSION_VALS) == 0:
            print('No data to save')
            return
        
        df = pd.DataFrame({
            'BN_cooler [MHz]': self.BN_vals,
            'PhotoEm [V]': self.PHOTOEMISSION_VALS,
            'Error [V]': self.ERRORS,
        })
        df.to_csv(f'Photoemission_Data_B={self.B}.csv')
    
    def SaveParameters(self):
        """
        Save the Lorentzian fit parameters to a CSV file.
        """
        if self.V0 == 0:
            print('No fit has been performed')
            return
        
        df = pd.DataFrame({
            'V0 [V]': [self.V0],
            'dV0 [V]': [self.dV0],
            'BN0 [MHz]': [self.BN0],
            'dBN0 [MHz]': [self.dBN0],
            'Gamma [MHz]': [self.Gamma],
            'dGamma [MHz]': [self.dGamma],
        })
        df.to_csv(f'Lorentzian_params_B={self.B}.csv')