import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
from preformation_plot import Plotter

'''
Script that handles postprocessing and plotting of all GITT files.
'''

class Gitt(Plotter):
    def __init__(self):
        Plotter.__init__(self) # Inherit variables from Plotter.
        self.plotter_inst = Plotter() # Instantiate this class.
        
    def gitt_splitter(self, filename, df, cells ='all', show = False, csvs = True, plots = True):
        # First stage split: determines charge and discharge in sequential profiles.
        # self.df_dis = df[df['I[A/kg]' < 0] # Excludes charging steps
        self.end_rows = df.index[df['State'] == 110].tolist()                 
        print('\n\n*********self.endrows=', self.end_rows,'**********\n\n')
        self.df_dis = df.iloc[:(self.end_rows[0] + 1)]
        self.df_ch = df.iloc[(self.end_rows[0] + 1):(self.end_rows[1] + 1)]
        self.mpl_formatting(plot_type = 'transient')        
        if plots:
            plt.plot(self.df_dis['Time[h]'],self.df_dis['U[V]'],label='discharge, cell '+ self.cell_ID)
            plt.plot(self.df_ch['Time[h]'],self.df_ch['U[V]'],label='charge, cell '+self.cell_ID)
            
            plt.legend()            
            plt.savefig('plot_output/' + filename.replace('.txt','') + '.png')
            plt.clf()
        if csvs:    
            self.df_dis.to_csv('csv_output/' + filename.replace('.txt','') + '_GITTdis.csv')
            self.df_ch.to_csv('csv_output/' + filename.replace('.txt','') + '_GITTch.csv')
        if show:
            plt.show()

    def gitt_process(self, filename, cells ='all', show = False, csvs = True, plots = True):
        # Performs splitting and analysis to get the GITT profiles (true OCV). 
        self.df_dis_ocv = self.df_dis[self.df_dis['I[A/kg]'] == 0] # First get OCV data.
        self.df_ch_ocv = self.df_ch[self.df_ch['I[A/kg]'] == 0]                                      
        self.gitt_dis_dict ={}
        self.gitt_ch_dict ={}
        self.gitt_dis_x = []
        self.gitt_ch_x = []
        self.gitt_dis_OCV = []
        self.gitt_ch_OCV = []
        self.gitt_dis_df = pd.DataFrame()
        self.gitt_ch_df = pd.DataFrame()
        self.ch_max_cap = self.df_ch['Ah-Ch[Ah/kg]'].max() # Maximum capacity obtained during charge.
                                    
        for n in range(1, self.df_dis_ocv['Count'].max() + 1):
            self.dis_current = self.df_dis_ocv[self.df_dis_ocv['Count'] == n]                        
            self.gitt_dis_dict[n] = self.dis_current
            self.gitt_dis_x.append(self.dis_current['Ah-Dis[Ah/kg]'].iloc[2]/self.max_cap) # Normalise to max cap (theoretical)
            self.mean_ocv = self.dis_current['U[V]'].iloc[-21:-1].mean() # Average last 20 points, excluding last.                        
            self.gitt_dis_OCV.append(self.mean_ocv)                                    
                                    
        for n in range(1, self.df_ch_ocv['Count'].max() + 1):
            self.ch_current = self.df_ch_ocv[self.df_ch_ocv['Count'] == n]                       
            self.gitt_ch_dict[n] = self.ch_current
            self.gitt_ch_x.append((self.ch_max_cap - self.ch_current['Ah-Ch[Ah/kg]'].iloc[2])/self.max_cap)
            self.mean_ocv = self.ch_current['U[V]'].iloc[-21:-1].mean() # Average last 20 points, excluding last.                        
            self.gitt_ch_OCV.append(self.mean_ocv)

        self.gitt_dis_df['SOC'] = self.gitt_dis_x
        self.gitt_dis_df['OCV'] = self.gitt_dis_OCV
        self.gitt_ch_df['SOC'] = self.gitt_ch_x
        self.gitt_ch_df['OCV'] = self.gitt_ch_OCV

        print('\n\n******Dataframe', self.gitt_dis_df, '********\n\n')
        self.mpl_formatting(plot_type='GITT')
                                    
        if plots:
            plt.plot(self.gitt_dis_df['SOC'],self.gitt_dis_df['OCV'],label='discharge GITT, cell '+ self.cell_ID)
            plt.plot(self.gitt_ch_df['SOC'],self.gitt_ch_df['OCV'],label='charge GITT, cell '+self.cell_ID)
            
            plt.legend()            
            plt.savefig('plot_output/' + filename.replace('.txt','') + '_GITT.png')
            plt.show()                        
            plt.clf()
                         
    def gitt(self, key, cells='all',show = False, csvs = True, plots = True, split = True, process = True, relaxation = False):
        # iterates through the experimental GITT files.
        self.show = show
        self.csvs = csvs
        self.plots = plots
        self.gitt_df_dict = {}
        plt.clf()
        
        if cells == 'all':
            self.all_cells = True
            self.cell_list = [i for i in range(256)]
        else:
            self.all_cells = False
            self.cell_list = cells

        for filename in self.file_list_dict[key]:
            self.filename = filename
            self.file_wo_dir = self.filename.split('/')[-1] # Omits directory structure.
            print(filename)
            self.cell_ID = self.extract_cell_ID(filename)
            self.cell_number = self.cell_ID[-1]
            print(self.cell_ID, type(self.cell_ID))
            if relaxation:
                self.fig1 = plt.figure(1)
                self.fig2 = plt.figure(2)
                self.ax1 = self.fig1.add_subplot(111)
                self.ax2 = self.fig2.add_subplot(111)                
                # Plot only results as a function of relaxation time/
                self.points = self.file_wo_dir.split('_')[-1].split('p')[0]
                print(self.points, ' = points')
                if self.points == '50': # Plot only the coarse grid results.
                    self.plotter_inst.extract_dataframe(filename = self.filename, header = False, data_type='GITT') # Current working dataframe.
                    self.gitt_df = self.df_cleanup(self.plotter_inst.gitt_df)
                    self.gitt_splitter(filename = filename, df = self.gitt_df, show = False, plots = False, csvs = False)
                    self.gitt_process(filename = filename, show = False, plots = False, csvs = False)                    
                    self.relaxation_time = self.file_wo_dir.split('_')[-1].replace('relax.txt','').split('p')[-1]
                    self.ax1.plot(self.gitt_dis_df['SOC'],self.gitt_dis_df['OCV'],label='discharge GITT, cell '+ self.cell_ID + ', %s min relax' % self.relaxation_time)
                    self.ax2.plot(self.gitt_ch_df['SOC'],self.gitt_ch_df['OCV'],label='charge GITT, cell '+ self.cell_ID+ ', %s min relax' % self.relaxation_time)

            else:         
                if self.cell_number in self.cell_list or self.all_cells:
                    self.extract_dataframe(filename, header = False, data_type='GITT') # Current working dataframe.
                    self.gitt_df = self.df_cleanup(self.gitt_df)
                    self.gitt_df_dict[self.cell_ID] = self.gitt_df # Put into the dictionary.
                    self.mpl_formatting(plot_type = 'transient')
                    if split:
                        self.gitt_splitter(filename = filename, df = self.gitt_df)
                    if process:    
                        self.gitt_process(filename = filename)
        if relaxation:
            self.fig1.legend()
            self.fig2.legend()
            self.plotter_inst.mpl_formatting(plot_type = 'GITT')            
            self.fig1.savefig('plot_output/dis_relax.png')
            self.plotter_inst.mpl_formatting(plot_type = 'GITT')
            self.fig2.savefig('plot_output/ch_relax.png')
            plt.show()
                
if __name__ == '__main__':
    gitt_obj = Gitt()
    gitt_obj.gitt(key = 'GITT', csvs = False, show = False, relaxation = True)
        
