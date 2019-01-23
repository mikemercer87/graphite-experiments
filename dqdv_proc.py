import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
from preformation_plot import Plotter
from scipy.interpolate import interp1d

'''
Script that handles postprocessing to get dQ/dV files.
'''

class DQDV(Plotter):
    def __init__(self):
        Plotter.__init__(self) # Inherit variables from Plotter.
        self.plotter_inst = Plotter() # Instantiate this class.
        self.dqdv_dis = pd.DataFrame()
        self.dqdv_ch = pd.DataFrame()

    def cubic_smoother(self, df, option = 'SOC', n_points = 5000):
        # Returns two smoothed and interpolated grids, x and the variable. Npoints is the number of points to smooth over.
        df_smoothed = pd.DataFrame() # Empty dataframe for the smoothed data.
        
        if option == 'SOC':
            df_sorted = df.sort_values(by = option)
            min_SOC = df_sorted['SOC'].min()
            max_SOC = df_sorted['SOC'].max()
            new_x = np.linspace(start = min_SOC, stop = max_SOC, num = n_points)
            soc_array = df_sorted['SOC'].values
            dqdv_array = df_sorted['dx/dE'].values
            smooth_fn = interp1d(x = soc_array.tolist(), y = dqdv_array.tolist(), kind='linear')
            smooth= smooth_fn(new_x)
            df_smoothed['SOC_%d' % n_points] = new_x
            df_smoothed['dx/dE_%d' % n_points] = smooth

        return(df_smoothed)

    def data_reduction(self, slice_window = 10000):
        # Cuts down the number of data points prior to smoothing.
        # Generates new dataframe, uses high resolution OCV channels.
        self.window = slice_window

        self.dqdv_dis['SOC'] = self.log_rolling_average(self.df2_dis['SOC'],slice_window)
        self.dqdv_dis['U[V]'] = self.log_rolling_average(self.df2_dis['OCV'],slice_window)
        self.dqdv_ch['SOC'] = self.log_rolling_average(self.df2_ch['SOC'],slice_window)
        self.dqdv_ch['U[V]'] = self.log_rolling_average(self.df2_ch['OCV'],slice_window)
        
    def log_rolling_average(self, df_col, slice_window = 10000):
        rolling_avg = df_col.rolling(int(slice_window),center=True).mean().dropna()
        return(rolling_avg)
            
    def differentiate(self, df, ir = False, state='dis'):
        # Performs differentiation of the dataset.
        dx = float(df['SOC'].iloc[1]) - float(df['SOC'].iloc[0])
        dV = float(df['U[V]'].iloc[1]) - float(df['U[V]'].iloc[0])
        if ir and state == 'dis':
            print(df['U[V]'].iloc[100], '\n\n*********discharge, before correction******\n\n')
            df['U[V]'] = df['U[V]'] + self.ir
            print(df['U[V]'].iloc[100], '\n\n*********discharge, after correction******\n\n')
        elif ir and state == 'ch':
            print(df['U[V]'].iloc[100], '\n\n*********charge, before correction******\n\n')            
            df['U[V]'] = df['U[V]'] - self.ir
            print(df['U[V]'].iloc[100], '\n\n*********charge, after correction******\n\n')            
        deriv = np.gradient(df['U[V]'],edge_order=2)
        deriv2 = np.gradient(df['SOC'],edge_order=2)             

        df['dE/dx'] = -deriv / dx
        df['dx/dE'] = 1 / (df['dE/dx'])
        return(df)

    def read_files(self, key='Galvan_C50', cells = 'all'):
        # Export processed preformation into csv. If cells is entered as a list, only those ones are processed.
        # Plots U as a function of time, for all cycles.
        self.preform_df_dict = {}

        if cells == 'all':
            self.all_cells = True
            self.cell_list = [i for i in range(256)]
        else:
            self.all_cells = False
            self.cell_list = cells
            
        for filename in self.file_list_dict[key]:
            self.filename = filename            
            print(filename)
            self.cell_ID = self.plotter_inst.extract_cell_ID(filename)
            self.cell_number = self.cell_ID[-1]
            print(self.cell_ID, type(self.cell_ID))
            if self.cell_number in self.cell_list or self.all_cells:
                self.plotter_inst.mpl_formatting(plot_type = 'dQdV')
                self.plotter_inst.extract_dataframe(filename, header = False) # Current working dataframe.
                self.preform_df = self.plotter_inst.df_cleanup(self.plotter_inst.preform_df)
                self.preform_df_dict[self.cell_ID] = self.preform_df # Put into the dictionary.
                self.plotter_inst.basytec_split(self.preform_df, show = False, csvs = False, plots = False)
                self.df_dict = self.plotter_inst.df_splits # Dictionary containing all the split data for the current cycle.
                self.df2_dis = self.df_dict['2D'] # Gets out split dataframes
                self.df2_ch = self.df_dict['2C']
                self.mpl_formatting(plot_type = 'dQdV')
                self.plotter_inst.get_ir_drop()
                self.ir = self.plotter_inst.ir
                
                for window in [50]:
                    kwindow = window /1000
                    self.data_reduction(slice_window = window)               
                
                    self.dqdv_dis = self.differentiate(self.dqdv_dis,ir=False,state='dis')
                    self.dqdv_ch = self.differentiate(self.dqdv_ch,ir=False,state='ch')
                    self.dqdv_dis_ir = self.differentiate(self.dqdv_dis,ir=True,state='dis')
                    self.dqdv_ch_ir = self.differentiate(self.dqdv_ch,ir=True,state='ch')


                    self.ax1.plot(self.dqdv_dis_ir['SOC'], np.log(self.dqdv_dis_ir['dx/dE']), label='iR corrected, dcharge')
                    self.ax1.plot(self.dqdv_ch_ir['SOC'], np.log(self.dqdv_ch_ir['dx/dE']), label='iR corrected, charge')
                    self.ax2.plot(self.dqdv_dis['SOC'], np.log(self.dqdv_dis['dx/dE']), label='uncorrected, dcharge')
                    self.ax2.plot(self.dqdv_ch['SOC'], np.log(self.dqdv_ch['dx/dE']), label='uncorrected, charge')
                    self.ax3.plot(self.dqdv_dis_ir['U[V]'], np.log(self.dqdv_dis_ir['dx/dE']), label='iR corrected, dcharge')
                    self.ax3.plot(self.dqdv_ch_ir['U[V]'], np.log(self.dqdv_ch_ir['dx/dE']), label='iR corrected, charge')
                    self.ax4.plot(self.dqdv_dis_ir['U[V]'], np.log(self.dqdv_dis_ir['dx/dE']), label='uncorrected, dcharge')
                    self.ax4.plot(self.dqdv_ch_ir['U[V]'], np.log(self.dqdv_ch_ir['dx/dE']), label='uncorrected, dcharge')                       

#                self.ax1.set_title('Discharge')
#                self.ax2.set_title('Charge')
#                self.ax3.set_title('Discharge')
#                self.ax4.set_title('Charge')
                
                for fig_label in self.plotter_inst.fig_labels:    
                    fig_label.legend(loc=0)
                    fig_label.tight_layout()
                 
                self.fig1.savefig('dQdV_%s_iR.png' % self.cell_number)
                self.fig2.savefig('dQdV_%s.png' % self.cell_number)
                self.fig3.savefig('dQdV_V_%s_iR.png' % self.cell_number)
                self.fig4.savefig('dQdV_V_%s.png' % self.cell_number)                 
                plt.show()
                plt.clf()     

                
if __name__ == '__main__':
    dqdv = DQDV()
    dqdv.read_files()
