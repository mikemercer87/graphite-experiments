import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os
from preformation_plot import Plotter

'''
Script that handles postprocessing and plotting of all entropy profile files.
'''

class Entropy(Plotter):
    def __init__(self):
        Plotter.__init__(self) # Inherit variables from Plotter.
        self.plotter_inst = Plotter() # Instantiate this class.
        self.methods = ['M1','M2','M3','M4']

    def entropy_cleanup(self):
        # Generates input dataframe generated from Matlab scripts and renames the columns to allow more meaningful plots.
        self.input_df['OCV'] = self.input_df['OCV [V]   ']
        self.input_df['deltaS'] = self.input_df['Bestfit Entropy [J mol-1 K-1]']
        self.input_df['errupper'] = self.input_df['Bestfit Entropy_Upper [J mol-1 K-1]']
        self.input_df['errlower'] = self.input_df['Bestfit Entropy_Lower [J mol-1 K-1]']
        for method in self.methods:
            self.input_df['deltaS_%s' % method] = self.input_df['%s Entropy [J mol-1 K-1]' % method]
            self.input_df['errupper_%s' % method] = self.input_df['%s Entropy_Upper [J mol-1 K-1]' % method]
            self.input_df['errlower_%s' % method] = self.input_df['%s Entropy_Lower [J mol-1 K-1]' % method]
              
        self.input_df['Capacity[mAh]'] = self.input_df['Charge/Discharge [mAh]'] # Note this is not normalised to active material mass!
        self.input_df = self.input_df[self.input_df['errupper'] < 4]
        self.input_df = self.input_df[self.input_df['errupper'] < 4]        
                    
    def entropy_plot(self, cells ='all', show = False, csvs = True, plots = True, methodplots = False, datatype = 'Dis_entropy'):
        # Datatype is 'discharge' by default. But can be charge!
        self.file_dict = {} # Will contain neatly organised file structure containing a zipped list of the relevant folders.        
        if datatype == 'Dis_entropy':
            self.directory = self.disentropy_dir
            self.name = 'discharge'
            self.input_name = 'disentropy'
        elif datatype == 'Ch_entropy':
            self.directory = self.chentropy_dir
            self.name = 'charge'
            self.input_name = 'chentropy'

        self.cell_folders = os.listdir(self.directory)
        for folder in self.cell_folders:
            self.cell_ID = folder.split('_')[-1]
            for filename in os.listdir(self.directory + folder):
                if filename.endswith('.txt'):
                    self.basyfile = filename
                elif filename.endswith('entropy.csv'):
                    self.matfile = filename
            self.file_dict[self.cell_ID] = (self.basyfile,self.matfile) # Put into dicionary
            
        for cell, (basyfile,matfile) in self.file_dict.items():
            self.cell_ID = cell
            self.basypath = self.directory + 'cell_' + self.cell_ID + '/' + basyfile
            self.matpath = self.directory + 'cell_' + self.cell_ID + '/' + matfile
                               
            self.input_df = pd.read_csv(self.matpath)
            self.entropy_cleanup()
            plt.clf()
            self.plotter_inst.mpl_formatting(plot_type = datatype)
            plt.errorbar(self.input_df['OCV'][::2],self.input_df['deltaS'][::2],yerr=[self.input_df['errlower'][::2],self.input_df['errupper'][::2]],label = 'cell %s, %s' % (self.cell_ID,self.name), linestyle = '',marker='o')
            plt.legend()
            plt.savefig('plot_output/' + self.directory + '_' + self.cell_ID + '_%s' % self.name)
            if show:
                plt.show()
            plt.clf()
            for method in self.methods:
                if methodplots:
                    self.plotter_inst.mpl_formatting(plot_type = datatype)                        
                    plt.errorbar(self.input_df['OCV'][::2],self.input_df['deltaS_%s'%method][::2],yerr=[self.input_df['errlower_%s'%method][::2],self.input_df['errupper_%s'%method][::2]],label = 'cell %s, %s, method %s' % (self.cell_ID,self.name,method),linestyle = '',marker='o')
                    plt.legend()
                    plt.savefig('plot_output/' + self.disentropy_dir + '_' + self.cell_ID + '_%s_%s.png' % (method,self.name))
                    if show:
                        plt.show()

            self.plotter_inst.extract_dataframe(self.basypath,data_type= datatype,skiprows=12) # Extract Basytec file.
            self.basy_df = self.plotter_inst.entropy_df
            self.active_mass = ((self.basy_df['Ah[Ah]'].iloc[-1])/(self.basy_df['Ah[Ah/kg]'].iloc[-1])) * 1000 # Convert kg to g.
            self.input_df['Capacity'] = np.abs(self.input_df['Capacity[mAh]']) / self.active_mass # Allow for charge and discharge.
            self.input_df['SOC'] = self.input_df['Capacity'] / self.max_cap
            self.plotter_inst.mpl_formatting(plot_type = datatype)
            plt.xlim([0.0,1.0])
            plt.xlabel('SOC')
            plt.errorbar(self.input_df['SOC'][::2],self.input_df['deltaS'][::2],yerr=[self.input_df['errlower'][::2],self.input_df['errupper'][::2]],label = 'cell %s, %s' % (self.cell_ID,self.name), linestyle = '',marker='o')
            plt.legend()
            plt.savefig('plot_output/' + self.directory + '_' + self.cell_ID + '_%s_vsSOC' % self.name)
            if show:
                plt.show()
                    
                            

if __name__ == '__main__':
    entropy_obj = Entropy()
    entropy_obj.entropy_plot(datatype = 'Dis_entropy', csvs = False,methodplots=False,show=True)
    entropy_obj.entropy_plot(datatype = 'Ch_entropy', csvs = False,methodplots=False,show=True)    
        
