import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
from preformation_plot import Plotter

'''
Script that handles postprocessing to get dQ/dV files.
'''

class DQDV(Plotter):
    def __init__(self):
        Plotter.__init__(self) # Inherit variables from Plotter.
        self.plotter_inst = Plotter() # Instantiate this class.
        self.dqdv_dis = pd.DataFrame()
        self.dqdv_ch = pd.DataFrame()

    def cubic_smoother(self, df_x, df_y, n_points=100):
        # Returns two smoothed and interpolated grids, x and the variable. Npoints is the number of points to smooth over.
        df_xred=df_x[df_x.notnull()].as_matrix()
        df_yred=df_y[df_y.notnull()].as_matrix()
        data_points=zip(df_xred,df_yred)
        data_points=sorted(data_points, key=lambda point: point[0])
        data_T = zip(*data_points)
        print('data_points', data_points[:][0])
#    print df_xred[-10:-1]
#    print df_yred[-10:-1]
        new_x=np.linspace(data_T[0][0], data_T[0][-1], n_points)
        print(new_x)
        smooth= interp1d(data_T[0], data_T[1], kind='cubic')(new_x)
        return(new_x, smooth)

    def data_reduction(self, slice_window = 1000):
        # Cuts down the number of data points prior to smoothing.
        self.window = slice_window
        dfs_in = [self.df2_dis, self.df2_ch]
        dfs_out = [self.dqdv_dis, self.dqdv_ch]
        
        for df_in, df_out in zip(dfs_in, dfs_out):
            df_out['SOC'] = df_in['SOC'].rolling(int(slice_window)).mean()
            df_out['U[V]'] = df_in['U[V]'].rolling(int(slice_window)).mean()
            
    def differentiate(self, df):
        # Performs differentiation of the dataset.
        # Improve the function! 
        
        dx = float(df['SOC'].iloc[1]) - float(df['SOC'].iloc[0])
        diff = np.diff(df['U[V]'])
        for i in range(len(df['SOC'])):
            if i == 0:
                df['diff'].iloc[i] = diff[i]
            else:
                df['diff'].iloc[i] = diff[i - 1]
        df['dx/dE'] = 1 / (df['diff'] / dx)
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
                self.plotter_inst.extract_dataframe(filename, header = False) # Current working dataframe.
                self.preform_df = self.plotter_inst.df_cleanup(self.plotter_inst.preform_df)
                self.preform_df_dict[self.cell_ID] = self.preform_df # Put into the dictionary.
                self.plotter_inst.basytec_split(self.preform_df, show = False, csvs = False, plots = False)
                self.df_dict = self.plotter_inst.df_splits
                self.df2_dis = self.df_dict['2D'] # Gets out split dataframes
                self.df2_ch = self.df_dict['2C']
                self.data_reduction()
                self.dqdv_dis = self.differentiate(self.dqdv_dis)
                self.dqdv_ch = self.differentiate(self.dqdv_ch)
                smooth1_dis = self.cubic_smoother(self.df_dis['SOC'],self.df_dis['dxdE'])
                smooth2_dis = self.cubic_smoother(self.df_dis['U[V]'],self.df_dis['dxdE'])
                smooth1_ch = self.cubic_smoother(self.df_ch['SOC'],self.df_ch['dxdE'])
                smooth2_ch = self.cubic_smoother(self.df_ch['SOC'],self.df_ch['dxdE'])
                plt.clf()
                plt.plot(smooth1_dis,smooth2_dis)
                plt.show()
                
if __name__ == '__main__':
    dqdv = DQDV()
    dqdv.read_files()
