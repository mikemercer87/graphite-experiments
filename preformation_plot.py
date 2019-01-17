'''
Script for handling the importing of Basytec files from the graphite entropy project into an easily Python readable Pandas format and for generating Matplotlib plot files.

Assumes the directory structure has been preserved from the version in the Box folder.

Could serve as a good basis for batch processing of experimental data files.

Also available in my Github account.
'''

import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import os

class Plotter():
    # Handles plotting for galvan. Plus paths for all experimental data.
    def __init__(self):
        # Get out all files and put into a meaningful directory.'
        self.master_dir = '00_20micron/00_Electrochem_25degrees/'
        self.preform_dir = self.master_dir + '00_Preformation_C_10/'
        self.galvan_dir = self.master_dir + '01_Galvan_C_50/'
        self.gitt_dir = self.master_dir + '02_GITT/'
        self.disentropy_dir = self.master_dir + '03_Discharge_entropy/'
        self.chentropy_dir = self.master_dir + '04_Charge_entropy/'
        self.experimental_directories = [self.preform_dir,self.galvan_dir,self.gitt_dir,self.disentropy_dir,self.chentropy_dir]
        self.readable_keys = ['Preform','Galvan_C50','GITT','Dis_entropy','Ch_entropy']
        self.file_list_dict = {} # Creates empty dictionary. This will eventually be filled with keys with are the experimental data set and values are a list of the absolute path of each file.
        for key, directory in zip(self.readable_keys,self.experimental_directories):
            for filename in os.listdir(directory):
                self.file_list_dict[key] = [directory + filename for filename in os.listdir(directory)]
        self.mpl_formatting() # Default plot format.        
        self.df_splits = {} # Empty dictionary to put the split files into, from discharge. Keys 1D, 1C, 2D, 2C etc.
        self.df_gittrelaxation = {} # Put in gitt data from different relaxation times.
        self.split_keys = ['1D', '1C', '2D', '2C', '3D', '3C']
        for directory in self.experimental_directories:  # Check existence of output directories and update if necessary.  
            os.makedirs('csv_output/' + directory, exist_ok = True) 
            os.makedirs('plot_output/' + directory, exist_ok = True)
        self.max_cap = 372 # Sets theoretical capacity of graphite.
        
    def mpl_formatting(self,fontsize = 16, linewidth = 2, markersize = 6, plot_type = 'SOC'):
        # All the options for customising the appearance of the generated plots.
        plt.style.use('fivethirtyeight')
        font = {'weight' : 'normal',
                'size'   : fontsize}

        lines = {'linewidth': linewidth,
                 'markersize': markersize}
        
        mpl.rc('font', **font)
        mpl.rc('lines', **lines)
        if plot_type == 'transient':
            plt.xlabel('Time (h)')
            plt.ylabel('Voltage (V) vs. Li')
            plt.xlim()
            plt.ylim()
        elif plot_type == 'SOC':
            plt.xlabel('SOC')
            plt.ylabel('Voltage (V) vs. Li')
            plt.xlim([0.0,1.0])
            plt.ylim([0.0,0.6])
        elif plot_type == 'capacity':
            plt.xlabel('Capacity (mAh/g)')
            plt.ylabel('Voltage')
            plt.xlim()
            plt.ylim()
        elif plot_type == 'GITT':
            plt.xlabel('SOC')
            plt.ylabel('OCV (V) vs. Li')
            plt.xlim([0.0,1.0])
            plt.ylim([0.0,0.6])
        elif plot_type == 'Dis_entropy' or plot_type == 'Ch_entropy':
            plt.xlabel('OCV (V) vs. Li')
            plt.ylabel('DeltaS (J mol-1 K-1)')
            plt.xlim([0.0,0.6])
            plt.ylim([-15.0,20.0])

        plt.tight_layout()
        
    def cubic_smoother(self,df_x,df_y,n_points):
        # Returns two smoothed and interpolated grids, x and the variable. Npoints is the number of points to smooth over.
        df_xred=df_x[df_x.notnull()].as_matrix()
        df_yred=df_y[df_y.notnull()].as_matrix()
        data_points=zip(df_xred,df_yred)
        data_points=sorted(data_points, key=lambda point: point[0])
        data_T = zip(*data_points)
        print('data_points', data_points[:][0])
#    print df_xred[-10:-1]
#    print df_yred[-10:-1]
        new_x=np.linspace(data_T[0][0],data_T[0][-1],n_points)
        print(new_x)
        smooth= interp1d(data_T[0],data_T[1],kind='cubic')(new_x)
        return(new_x,smooth)

#    def dQdV(self):
        # Handles differentiation to get dQ/dV.
                
    def df_cleanup(self,df):
        # Generates extra column without the pesky comment character on the first line.
        for string in list(df):
            if string.endswith('Time[h]'):
                df['Time[h]'] = df[string]
            elif string.startswith('OCV'):
                df['OCV'] = df[string] / 1000 # Converts OCV from V to mV.
            elif string.startswith('MEM'):
                df['T' + string[4]] = df[string]
            elif string == 'Ah[Ah/kg]':
                df['Capacity'] = np.abs(df[string])
        return(df)
            
    def extract_cell_ID(self,filename,instance = 1):
        # Gets the cell ID through string splitting. Default instance extracts the second character.
        string = filename.split('/')[-1] # Separates the filename from the path.
        cell_ID = string.split('_')[instance] # Extract second element from the list.
        return(cell_ID)

    def extract_dataframe(self,filename,header=True,data_type='Preform'):
        if header == True:
            # File contains the Basytec header information.
            no_skipped = 13
        else:
            no_skipped = 0
            
        if data_type == 'Preform':    
            self.preform_df = pd.read_csv(filename,skiprows = no_skipped,encoding = "ISO-8859-1")
        elif data_type == 'GITT':
            self.gitt_df = pd.read_csv(filename,skiprows = no_skipped,encoding = "ISO-8859-1")
            
    def basytec_split(self, df, no_cycles = 2, show = False, csvs = True, plots = True):
        # Handles all the cycle splits. Charge and discharge combined.
        self.show = show
        self.df = df
        self.no_cycles = no_cycles
        self.csvs = csvs
        self.plots = plots
        self.df_dict = {}
        self.used_keys = self.split_keys[0 : (2 * no_cycles + 1)]
        self.n = 1
        
        while self.n <= self.no_cycles:
            self.df_split = self.df[self.df['Count'] == self.n]
            self.df_split['Capacity'] = self.df_split['Capacity'] - self.df_split['Capacity'].iloc[0] # Correction to only compute capacity on that cycle.
#            df_split['SOC'] = df_split['Capacity']/df_split['Max_cap'] # Normalise to theoretical capacity.
            if plots:
                plt.plot(self.df_split['Time[h]'],self.df_split['U[V]'],label='cycle number = ' + str(self.n))
                self.mpl_formatting(plot_type = 'transient')
                plt.legend()
                plt.savefig('plot_output/' + self.filename.replace('.txt','') + 'Vvst_%d.png' % self.n)
                plt.clf()
    
                plt.plot(self.df_split['Capacity'],self.df_split['U[V]'],label='cycle number = ' + str(self.n))
                self.mpl_formatting(plot_type = 'capacity')
                plt.legend()
                plt.savefig('plot_output/' + self.filename.replace('.txt','') + 'Vvscap_%d.png' % self.n)
                plt.clf()
            if csvs:    
                self.df_split.to_csv('csv_output/' + self.filename.replace('.txt','') + '_%d.csv' % self.n)
                
            if show:    
                plt.show()
                
            self.basytec_second_split(df = self.df_split, n = self.n, show = self.show, csvs = self.csvs, plots = self.plots)
            self.n += 1

    def basytec_second_split(self, df, n, show = False, csvs = True, plots = True):
        # Uses a split dataframe (1st cycle, 2nd cycle) as input, splits again into separate charge and discharge cycles.
        self.df_dis = df[df['I[A/kg]'] < 0]
        self.df_ch = df[df['I[A/kg]'] > 0]

        cap_min_dis = self.df_dis.min()['Capacity'] # Defines SOC = 0
        cap_min_ch = self.df_ch.min()['Capacity'] # Same for charge.

        self.df_dis['SOC'] = (self.df_dis['Capacity'] - cap_min_dis) / self.max_cap
        self.df_ch['SOC'] = (self.df_ch['Capacity'] - cap_min_ch) / self.max_cap

        if plots:
            plt.plot(self.df_dis['SOC'],self.df_dis['U[V]'],label='Discharge no %d' % n)
            self.mpl_formatting(plot_type = 'SOC')
            plt.legend()
            plt.savefig('plot_output/' + self.filename.replace('.txt','') + 'Vvssocdis_%d.png' % n)
            plt.clf()
            self.mpl_formatting(plot_type = 'SOC')            
            plt.plot(self.df_ch['SOC'],self.df_ch['U[V]'],label='Charge no %d' % n)
            plt.legend()
            plt.savefig('plot_output/' + self.filename.replace('.txt','') + 'Vvssocch_%d.png' % n)
            plt.clf()
            self.mpl_formatting(plot_type = 'SOC')
            plt.plot(self.df_dis['SOC'],self.df_dis['U[V]'],label='Discharge no %d' % n)            
            plt.plot(self.df_ch['SOC'],self.df_ch['U[V]'],label='Charge no %d' % n)
            plt.legend()
            plt.savefig('plot_output/' + self.filename.replace('.txt','') + 'Ch+Dis_%d.png' % n)
            plt.clf()            
        if show:
            plt.show()
        if csvs:
            self.df_dis.to_csv('csv_output/' + filename.replace('.txt','') + '_D%d.csv' % n)            
            self.df_ch.to_csv('csv_output/' + filename.replace('.txt','') + '_C%d.csv' % n)            

    def cycling(self, key, cells = 'all', show = False, csvs = True, plots = True):
        # Export processed preformation into csv. If cells is entered as a list, only those ones are processed.
        # Plots U as a function of time, for all cycles.
        self.show = show
        self.csvs = csvs
        self.plots = plots
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
            self.cell_ID = self.extract_cell_ID(filename)
            self.cell_number = self.cell_ID[-1]
            print(self.cell_ID, type(self.cell_ID))
            if self.cell_number in self.cell_list or self.all_cells:
                self.extract_dataframe(filename, header = False) # Current working dataframe.
                self.preform_df = self.df_cleanup(self.preform_df)
                self.preform_df_dict[self.cell_ID] = self.preform_df # Put into the dictionary.
                plt.plot(self.preform_df['Time[h]'],self.preform_df['U[V]'],label='cell '+ self.cell_ID)
                self.mpl_formatting(plot_type = 'transient')
                plt.legend()
                if plots:
                    plt.savefig('plot_output/' + filename.replace('.txt','') + '.png')
                
                if show:
                    plt.show()

                if csvs:
                    self.preform_df.to_csv('csv_output/' + filename.replace('.txt','') + '.csv')
                plt.clf()
                self.basytec_split(self.preform_df, show = self.show, csvs = self.csvs, plots = self.plots)

#    def galvan_C50(self, show = False, csvs = True, plots = True):
          
        # Extract all information from GITT files.
        # Get C rate.
#    def entropy_discharge(self, csvs=True, plots=True):                
#    def create_all_files(self,csvs=True,plots=True):
        # Catch all function that will plot absolutely everything!

if __name__ == '__main__':
    plots = Plotter()
#    plots.cycling(key = 'Preform', csvs=False)
#    plots.cycling(key = 'Galvan_C50', csvs=False)

