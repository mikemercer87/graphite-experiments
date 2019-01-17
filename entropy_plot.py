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
        print(list(self.input_df),'**********keys***********')
        self.input_df['OCV'] = self.input_df['OCV [V]   ']
        self.input_df['deltaS'] = self.input_df['Bestfit Entropy [J mol-1 K-1]']
        self.input_df['errupper'] = self.input_df['Bestfit Entropy_Upper [J mol-1 K-1]']
        self.input_df['errlower'] = self.input_df['Bestfit Entropy_Lower [J mol-1 K-1]']
        for method in self.methods:
            self.input_df['deltaS_%s' % method] = self.input_df['%s Entropy [J mol-1 K-1]' % method]
            self.input_df['errupper_%s' % method] = self.input_df['%s Entropy_Upper [J mol-1 K-1]' % method]
            self.input_df['errlower_%s' % method] = self.input_df['%s Entropy_Lower [J mol-1 K-1]' % method]
              
        self.input_df['Capacity'] = self.input_df['Charge/Discharge [mAh]'] # Note this is not normalised to active material mass!
        self.input_df = self.input_df[self.input_df['errupper'] < 4]
        self.input_df = self.input_df[self.input_df['errupper'] < 4]        
        print(self.input_df,'************datafram************')
        self.cell_ID = self.plotter_inst.extract_cell_ID(self.full_path)
                    
    def entropy_plot(self, cells ='all', show = False, csvs = True, plots = True, methodplots = False, datatype = 'Dis_entropy'):
        # Datatype is 'discharge' by default. But can be charge!
        if datatype == 'Dis_entropy':
            for filename in os.listdir(self.disentropy_dir):
                self.full_path = self.disentropy_dir + filename
                if self.full_path.endswith('entropy.csv'):
                    self.input_df = pd.read_csv(self.full_path)
                    self.entropy_cleanup()
                    plt.clf()
                    self.plotter_inst.mpl_formatting(plot_type = 'Dis_entropy')
                    plt.errorbar(self.input_df['OCV'][::2],self.input_df['deltaS'][::2],yerr=[self.input_df['errlower'][::2],self.input_df['errupper'][::2]],label = 'cell %s, discharge' % self.cell_ID,linestyle = '',marker='o')
                    plt.legend()
                    plt.savefig('plot_output/' + self.disentropy_dir + '_' + self.cell_ID + '_discharge.png')
                    if show:
                        plt.show()
                    plt.clf()
                    for method in self.methods:
                        if methodplots:
                            self.plotter_inst.mpl_formatting(plot_type = 'Dis_entropy')                        
                            plt.errorbar(self.input_df['OCV'][::2],self.input_df['deltaS_%s'%method][::2],yerr=[self.input_df['errlower_%s'%method][::2],self.input_df['errupper_%s'%method][::2]],label = 'cell %s, discharge, method %s' % (self.cell_ID,method),linestyle = '',marker='o')
                            plt.legend()
                            plt.savefig('plot_output/' + self.disentropy_dir + '_' + self.cell_ID + '_discharge_%s.png' % method)
                        if show:
                            plt.show()
                    
        elif datatype == 'Ch_entropy':
            for filename in os.listdir(self.chentropy_dir):
                self.full_path = self.chentropy_dir + filename                
                if self.full_path.endswith('entropy.csv'):
                    self.input_df = pd.read_csv(self.full_path)
                    self.entropy_cleanup()
                    plt.clf()
                    if plots:
                        self.plotter_inst.mpl_formatting(plot_type = 'Ch_entropy')                    
                        plt.errorbar(self.input_df['OCV'][::2],self.input_df['deltaS'][::2],yerr=[self.input_df['errlower'][::2],self.input_df['errupper'][::2]],label = 'cell %s, charge' % self.cell_ID,linestyle='',marker='o')
                        plt.legend()
                        plt.savefig('plot_output/' + self.chentropy_dir + '_' + self.cell_ID + '_charge.png')
                    if show:
                         plt.show()
                    plt.clf()
                    for method in self.methods:
                        if methodplots:
                            self.plotter_inst.mpl_formatting(plot_type = 'Ch_entropy')                        
                            plt.errorbar(self.input_df['OCV'][::2],self.input_df['deltaS_%s'%method][::2],yerr=[self.input_df['errlower_%s'%method][::2],self.input_df['errupper_%s'%method][::2]],label = 'cell %s, charge, method %s' % (self.cell_ID,method),linestyle = '',marker='o')
                            plt.legend()
                            plt.savefig('plot_output/' + self.chentropy_dir + '_' + self.cell_ID + '_charge_%s.png' % method)
                        if show:
                            plt.show()

        

if __name__ == '__main__':
    entropy_obj = Entropy()
    entropy_obj.entropy_plot(datatype = 'Dis_entropy', csvs = False,methodplots=True)
    entropy_obj.entropy_plot(datatype = 'Ch_entropy', csvs = False,methodplots=True)    
        
