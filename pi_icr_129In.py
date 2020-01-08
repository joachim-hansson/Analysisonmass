#A script reads number from a .csv file. Zhihao Gao
#from PyQt5 import QtWidgets
from pyqtgraph.Qt import QtGui,QtCore
import pyqtgraph as pg
import pyqtgraph.exporters
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from scipy.optimize import curve_fit 
from scipy import interpolate
import sys
import re
import helpers
import ast
import pyqtgraph.opengl as gl
import gaussfitter2
import math
from sklearn import mixture
from sklearn.utils import resample

x_velocity = 2*0.00066
y_velocity = 2*0.00066
ps2ns = 1000
ps2us = 1000000
h_x_cent = -1.364
h_y_cent = 0.909
x_center = -0.513
y_center = -0.134
#QtGui.QApplication.setGraphicsSystem('raster')
app = QtGui.QApplication([])
win =pg.GraphicsWindow(title="Basic plotting")
win.resize(1600,800)
win.setWindowTitle('TOF histograms')
win_1 =pg.GraphicsWindow(title="Basic plotting")
win_1.resize(1000,800)
win_1.setWindowTitle('Position plots')
pg.setConfigOptions(antialias=True)

bunch_data_raw = None
ions_data_raw = None

def load_csv():
    #filepath = "../Measurement_IGISOL_2019/28/20_14_i134_t_acc_226ms_1/"# data of 134I
    filepath = "../Measurement_IGISOL_2019/30/01_25_in129_T2_157ms_1/"# data of 129In
    bunch_filename = filepath + "bunches.csv"
    ions_filename = filepath + "ions.csv"
    bunch_data_raw = pd.read_csv(bunch_filename)
    if 'Unnamed: 0' in bunch_data_raw.columns:
        bunch_data_raw= bunch_data_raw.drop(['Unnamed: 0'],axis=1)
    bunch_empty = bunch_data_raw[bunch_data_raw['ions per bunch']==0]
    print("Total number of bunch =",len(bunch_data_raw))

    ions_data_raw = pd.read_csv(ions_filename)
    if 'Unnamed: 0' in ions_data_raw.columns:
        ions_data_raw = ions_data_raw.drop(['Unnamed: 0'],axis=1)
    
    merged_raw = pd.merge(bunch_data_raw,ions_data_raw,on="bunch number")
 #   filepath_1 = "../Measurement_IGISOL_2019/28/21_01_i134_t_acc_226ms_1/"# Path of first extra data
#    merged_raw = merged_raw.append(read_csv_ex(filepath_1))  # Mergering extra data
    print("Total number of ions = ", len(merged_raw))

    max_ions_per_bunch = merged_raw['ions per bunch'].max()
#    print(merged_raw[merged_raw['ions per bunch']==max_ions_per_bunch])
    print("Maximum of ions per bunch",max_ions_per_bunch)
    #Total_ions = len(merged_raw)
    merged_raw['delay_x1'] = merged_raw['tof_2']-merged_raw['tof_1']
    merged_raw['delay_x2'] = merged_raw['tof_3']-merged_raw['tof_1']
    merged_raw['delay_y1'] = merged_raw['tof_4']-merged_raw['tof_1']
    merged_raw['delay_y2'] = merged_raw['tof_5']-merged_raw['tof_1']
    merged_raw['delay_tof_sum_x'] = merged_raw['delay_x1']+merged_raw['delay_x2']
    merged_raw['delay_tof_sum_y'] = merged_raw['delay_y1']+merged_raw['delay_y2']
    merged_raw['tof_1'] = merged_raw['tof_1']/ps2us
        #print(merged.at[1,'delay_x2'])
#    print("Shape of the analysing data",merged_raw.shape[0])
        #print(merged.head(merged.shape[0]))
        
    plt_tof = win.addPlot(title="TOF spectrum (us)")
    plt_tof_x = win.addPlot(title="TOF of x coordinate (ns)")
    plt_tof_x1 = win.addPlot(title="TOF of x1 coordinate (ns)")
    plt_tof_x2 = win.addPlot(title="TOF of x2 coordinate (ns)")
    plt_bunch = win.addPlot(title="Ions per bunch")
    win.nextRow()
    plt_tof_cut = win.addPlot(title="TOF spectrum with cut (us)")
    plt_tof_y = win.addPlot(title="TOF of y coordinate (ns)")
    plt_tof_y1 = win.addPlot(title="TOF of y1 coordinate (ns)")
    plt_tof_y2 = win.addPlot(title="TOF of y2 coordinate (ns)")

    tof_y,tof_x = np.histogram(merged_raw['tof_1'],bins=np.linspace(30,55,500))
    plt_tof.plot(tof_x,tof_y,stepMode=True,fillLevel=0,brush=(0,0,255,150))
    lr = pg.LinearRegionItem([44,48])
    lr.setZValue(-10)
    plt_tof.addItem(lr)
    plt_tof_cut.plot(tof_x,tof_y,stepMode=True,fillLevel=0,brush=(0,0,255,150))
    def updatePlot():
        plt_tof_cut.setXRange(*lr.getRegion(),padding=0)
    def updateRegion():
        lr.setRegion(plt_tof_cut.getViewBox().viewRange()[0])
    lr.sigRegionChanged.connect(updatePlot)
    plt_tof_cut.sigXRangeChanged.connect(updateRegion)
    updatePlot()

    tof_sumx_y,tof_sumx_x = np.histogram(merged_raw['delay_tof_sum_x']/ps2ns,bins=np.linspace(40,48,100))
    plt_tof_x.plot(tof_sumx_x,tof_sumx_y,stepMode=True,fillLevel=0,brush=(0,0,255,150))
    tof_sumx_y1,tof_sumx_x1 = np.histogram(merged_raw['delay_x1']/ps2ns,bins=np.linspace(0,50,200))
    plt_tof_x1.plot(tof_sumx_x1,tof_sumx_y1,stepMode=True,fillLevel=0,brush=(0,0,255,150))
    tof_sumx_y2,tof_sumx_x2 = np.histogram(merged_raw['delay_x2']/ps2ns,bins=np.linspace(0,50,200))
    plt_tof_x2.plot(tof_sumx_x2,tof_sumx_y2,stepMode=True,fillLevel=0,brush=(0,0,255,150))
    tof_sumy_y,tof_sumy_x = np.histogram(merged_raw['delay_tof_sum_y']/ps2ns,bins=np.linspace(40,48,100))
    plt_tof_y.plot(tof_sumy_x,tof_sumy_y,stepMode=True,fillLevel=0,brush=(0,0,255,150))
    tof_sumy_y1,tof_sumy_x1 = np.histogram(merged_raw['delay_y1']/ps2ns,bins=np.linspace(0,50,200))
    plt_tof_y1.plot(tof_sumy_x1,tof_sumy_y1,stepMode=True,fillLevel=0,brush=(0,0,255,150))
    tof_sumy_y2,tof_sumy_x2 = np.histogram(merged_raw['delay_y2']/ps2ns,bins=np.linspace(0,50,200))
    plt_tof_y2.plot(tof_sumy_x2,tof_sumy_y2,stepMode=True,fillLevel=0,brush=(0,0,255,150))
    print("Total counts before tof cut",len(merged_raw))
    #TOF cut
    merged_raw = merged_raw[merged_raw['tof_1']>43]
    merged_raw = merged_raw[merged_raw['tof_1']<47]
    merged_raw = merged_raw.reset_index(drop=True)
    Total_ions = len(merged_raw)
    bunch_nums =[]
    bunch_i =[]
    data_name = []
    for i in range(max_ions_per_bunch+1):
        data_name.append(merged_raw[merged_raw['ions per bunch']==i])
        bunch_nums.append(len(data_name[i]))
        bunch_i.append(i)
    #    print("Number of ions bunch_\n",i,bunch_nums[i])
    bunch_i.append((max_ions_per_bunch+1))
    plt_bunch.plot(bunch_i,bunch_nums,stepMode=True,fillLevel=0,brush=(0,0,255,150))
    bunch_nums.append(Total_ions)
#    bunch_i.append((max_ions_per_bunch+1))
    print("Total number of ions after a primary cut,tof,scanning pars",bunch_nums[max_ions_per_bunch+1])
#    plt_bunch.plot(bunch_i,bunch_nums,stepMode=True,fillLevel=0,brush=(0,0,255,150))
    
    #for j in range(1,max_ions_per_bunch+1):
#        print("*****Starting analysis of bunch",j)
#        ana_process(data_name[j])
#        print("*****Analysis of bunch",j, "is over*****")
    ana_process(merged_raw)

def read_csv_ex(filepath_extra):
#    filepath_extra = "../Measurement_IGISOL_2019/28/21_01_i134_t_acc_226ms_1/"
    bunch_filename_ex = filepath_extra + "bunches.csv"
    ions_filename_ex = filepath_extra + "ions.csv"
    bunch_data_raw_ex = pd.read_csv(bunch_filename_ex)
    if 'Unnamed: 0' in bunch_data_raw_ex.columns:
        bunch_data_raw_ex= bunch_data_raw_ex.drop(['Unnamed: 0'],axis=1)
    ions_data_raw_ex = pd.read_csv(ions_filename_ex)
    if 'Unnamed: 0' in ions_data_raw_ex.columns:
        ions_data_raw_ex = ions_data_raw_ex.drop(['Unnamed: 0'],axis=1)

    raw_data_ex = pd.merge(bunch_data_raw_ex,ions_data_raw_ex,on="bunch number")
    return raw_data_ex
#    print(merged_raw)

def load_homo_data(filename):

    homo_in_data = pd.read_csv(filename)
    if 'Unnamed: 0' in homo_in_data.columns:
        homo_in_data = homo_in_data.drop(['Unnamed: 0'],axis=1)

    return homo_in_data

def ana_process(merged):

   # merged = pd.merge(bunch_data_raw,ions_data_raw,on="bunch number",how="left")
    empty_data = pd.DataFrame({'delay_x_sum':[np.nan],'delay_y_sum':[np.nan],'x_pos_d':[np.nan],'y_pos_d':[np.nan],'x_pos_s':[np.nan],'y_pos_s':[np.nan],'x_pos':[np.nan],'y_pos':[np.nan]})
    merged = merged.join(empty_data,on="bunch number")
#    merged.reset_index(drop=True)
        
    sum_x = []
    sum_y = []
    pos_x = []
    pos_y = []
    pos_x_d = []
    pos_y_d = []
    pos_x_s = []
    pos_y_s = []
    tof_region = []
    tof_c=40000
    x_sum_center = 43745/2#For 129In
    y_sum_center = 43511/2#For 129In
#    x_sum_center = 43863/2#For Mirror 136I
 #   y_sum_center =435582/2#For Mirror 136I
#    print(merged)
    n_double_s =0
    n_err=0
#    print("Index = ",merged.index) 
    for k in merged.index:
#        print("k = ",k)
#        if merged.at[k,'tof_1']>43 and merged.at[k,'tof_1']<48:
        if merged.at[k,'delay_x1'] in range(10000,tof_c) and merged.at[k,'delay_x2'] in range(10000,tof_c) and merged.at[k,'delay_y1'] in range(10000,tof_c) and merged.at[k,'delay_y2'] in range(10000,tof_c):
            merged.at[k,'x_pos_d'] = (merged.at[k,'delay_x1']-merged.at[k,'delay_tof_sum_x']/2)*x_velocity
            merged.at[k,'y_pos_d'] = (merged.at[k,'delay_y1']-merged.at[k,'delay_tof_sum_y']/2)*y_velocity
            merged.at[k,'x_pos'] = (merged.at[k,'delay_x1']-merged.at[k,'delay_tof_sum_x']/2)*x_velocity
            merged.at[k,'y_pos'] = (merged.at[k,'delay_y1']-merged.at[k,'delay_tof_sum_y']/2)*y_velocity
            merged.at[k,'delay_x_sum'] = merged.at[k,'delay_tof_sum_x']
            merged.at[k,'delay_y_sum'] = merged.at[k,'delay_tof_sum_y']
            n_double_s +=1
        elif (merged.at[k,'delay_x1'] in range(10000,tof_c) or merged.at[k,'delay_x2'] in range(10000,tof_c)) and (merged.at[k,'delay_y1'] in range(10000,tof_c) or merged.at[k,'delay_y2'] in range(10000,tof_c)):
            n_err +=1
            if merged.at[k,'delay_x1'] in range(10000,tof_c) and merged.at[k,'delay_x2'] in range(10000,tof_c):
                merged.at[k,'delay_x_sum'] = merged.at[k,'delay_tof_sum_x']
                merged.at[k,'x_pos_s'] = (merged.at[k,'delay_x1']-merged.at[k,'delay_tof_sum_x']/2)*x_velocity
                merged.at[k,'x_pos'] = (merged.at[k,'delay_x1']-merged.at[k,'delay_tof_sum_x']/2)*x_velocity
                if merged.at[k,'delay_y1'] not in range(10000,tof_c) and merged.at[k,'delay_y2'] in range(10000,tof_c):
                    merged.at[k,'y_pos'] = (y_sum_center-merged.at[k,'delay_y2'])*y_velocity
                    merged.at[k,'y_pos_s'] = (y_sum_center-merged.at[k,'delay_y2'])*y_velocity
                elif merged.at[k,'delay_y1'] in range(10000,tof_c) and merged.at[k,'delay_y2'] not in range(10000,tof_c):
                    merged.at[k,'y_pos'] = (merged.at[k,'delay_y1']-y_sum_center)*y_velocity
                    merged.at[k,'y_pos_s'] = (merged.at[k,'delay_y1']-y_sum_center)*y_velocity
            if merged.at[k,'delay_x1'] in range(10000,tof_c) and merged.at[k,'delay_x2'] not in range(10000,tof_c):
                merged.at[k,'x_pos'] = (merged.at[k,'delay_x1']-x_sum_center)*x_velocity
                merged.at[k,'x_pos_s'] = (merged.at[k,'delay_x1']-x_sum_center)*x_velocity
                if merged.at[k,'delay_y1'] in range(10000,tof_c) and merged.at[k,'delay_y2'] in range(10000,tof_c):
                    merged.at[k,'y_pos_s'] = (merged.at[k,'delay_y1']-merged.at[k,'delay_tof_sum_y']/2)*y_velocity
                    merged.at[k,'y_pos'] = (merged.at[k,'delay_y1']-merged.at[k,'delay_tof_sum_y']/2)*y_velocity
                    merged.at[k,'delay_y_sum'] = merged.at[k,'delay_tof_sum_y']
                elif merged.at[k,'delay_y1'] not in range(10000,tof_c) and merged.at[k,'delay_y2'] in range(10000,tof_c):
                    merged.at[k,'y_pos'] = (y_sum_center-merged.at[k,'delay_y2'])*y_velocity
                    merged.at[k,'y_pos_s'] = (y_sum_center-merged.at[k,'delay_y2'])*y_velocity
                elif merged.at[k,'delay_y1'] in range(10000,tof_c) and merged.at[k,'delay_y2'] not in range(10000,tof_c):
                    merged.at[k,'y_pos'] = (merged.at[k,'delay_y1']-y_sum_center)*y_velocity
                    merged.at[k,'y_pos_s'] = (merged.at[k,'delay_y1']-y_sum_center)*y_velocity
            if merged.at[k,'delay_x2'] in range(10000,tof_c) and merged.at[k,'delay_x1'] not in range(10000,tof_c):
                merged.at[k,'x_pos'] = (x_sum_center-merged.at[k,'delay_x2'])*x_velocity
                merged.at[k,'x_pos_s'] = (x_sum_center-merged.at[k,'delay_x2'])*x_velocity
                if merged.at[k,'delay_y1'] in range(10000,tof_c) and merged.at[k,'delay_y2'] in range(10000,tof_c):
                    merged.at[k,'y_pos_s'] = (merged.at[k,'delay_y1']-merged.at[k,'delay_tof_sum_y']/2)*y_velocity
                    merged.at[k,'y_pos'] = (merged.at[k,'delay_y1']-merged.at[k,'delay_tof_sum_y']/2)*y_velocity
                    merged.at[k,'delay_y_sum'] = merged.at[k,'delay_tof_sum_y']
                elif merged.at[k,'delay_y1'] not in range(10000,tof_c) and merged.at[k,'delay_y2'] in range(10000,tof_c):
                    merged.at[k,'y_pos'] = (y_sum_center-merged.at[k,'delay_y2'])*y_velocity
                    merged.at[k,'y_pos_s'] = (y_sum_center-merged.at[k,'delay_y2'])*y_velocity
                elif merged.at[k,'delay_y1'] in range(10000,tof_c) and merged.at[k,'delay_y2'] not in range(10000,tof_c):
                    merged.at[k,'y_pos'] = (merged.at[k,'delay_y1']-y_sum_center)*y_velocity
                    merged.at[k,'y_pos_s'] = (merged.at[k,'delay_y1']-y_sum_center)*y_velocity

    tof_x_sum = merged['delay_x_sum']#.dropna().tolist()
    fitting_hist(tof_x_sum,40,50)
    tof_y_sum = merged['delay_y_sum']#.dropna().tolist()
    fitting_hist(tof_y_sum,40,50)
    print("Number of all ions",len(merged))
    print("Double signals for both of X and Y",n_double_s)
    print("Losing one signal for X or Y or Both",n_err)
    
#    if len(merged)>0:
#        for j in merged.index:
#            if merged.at[j,'y_pos_d']>=0:
#                merged.at[j,'angle'] = np.arctan2(merged.at[j,'y_pos_d']-y_center,merged.at[j,'x_pos_d']-x_center)*180/np.pi
#            else:
#                merged.at[j,'angle'] = np.arctan2(merged.at[j,'y_pos_d']-y_center,merged.at[j,'x_pos_d']-x_center)*180/np.pi+360
#            merged.at[j,'radius'] = np.sqrt( np.square(x_center-merged.at[j,'x_pos_d'])+np.square(y_center-merged.at[j,'y_pos_d']))
#            if merged.at[j,'angle']<60:
#                merged.at[j,'angle'] = merged.at[j,'angle']-60+360
#            else:
#                merged.at[j,'anlge'] = merged.at[j,'angle']-60
#    else:
#        merged.at[j,'angle'] = np.nan
#        merged.at[j,'radius'] = np.nan
    if len(merged)>0:
        merged['angle'] = merged.apply(
                lambda row: np.arctan2(row.y_pos - y_center,row.x_pos-x_center)*180/np.pi if (row.y_pos - y_center >=0) else np.arctan2(row.y_pos-y_center, row.x_pos-x_center)*180/np.pi + 360,axis =1)
        merged['radius'] = np.sqrt(np.square(merged['x_pos']-x_center)+np.square(y_center-merged['y_pos']))
        merged['radius_cor'] = np.sqrt(np.square(merged['x_pos']-h_x_cent)+np.square(h_y_cent-merged['y_pos']))
        merged['angle_cor'] = merged.apply(
                lambda row: np.arctan2(row.y_pos -h_y_cent,row.x_pos-h_x_cent)*180/np.pi if (row.y_pos - h_y_cent >=0) else np.arctan2( row.y_pos-h_y_cent, row.x_pos-h_x_cent)*180/np.pi + 360,axis =1)
    else:
        merged['angle'] = np.nan
        merged['radius'] = np.nan
        merged['angle_cor'] = np.nan
        merged['radius_cor'] = np.nan
        
#    print(merged)
    plt_pos_d = win_1.addPlot(title="Position on MCP using double signals")
    plt_pos_d.setRange(xRange=(-25,25),yRange=(-25,25))
    plt_pos_d.setLabel('left',text="Y position",units='mm')
    plt_pos_d.setLabel('bottom',text="X position",units='mm')
    pos_d = pg.ScatterPlotItem(size=2,pen=pg.mkPen(None),brush=pg.mkBrush(0,0,255,255))
    pos_d.addPoints(merged['x_pos_d'],merged['y_pos_d'])
    plt_pos_d.addItem(pos_d)
    
    plt_pos = win_1.addPlot(title="Position on MCP using all data")
    plt_pos.setRange(xRange=(-25,25),yRange=(-25,25))
    plt_pos.setLabel('left',text="Y position",units='mm')
    plt_pos.setLabel('bottom',text="X position",units='mm')
    pos = pg.ScatterPlotItem(size=2,pen=pg.mkPen(None),brush=pg.mkBrush(0,0,255,255))
    pos.addPoints(merged['x_pos'],merged['y_pos'])
    plt_pos.addItem(pos)
    
    plt_pos_s = win_1.addPlot(title="Position on MCP using single signals")
    plt_pos_s.setRange(xRange=(-25,25),yRange=(-25,25))
    plt_pos_s.setLabel('left',text="Y position",units='mm')
    plt_pos_s.setLabel('bottom',text="X position",units='mm')
    pos_s = pg.ScatterPlotItem(size=2,pen=pg.mkPen(None),brush=pg.mkBrush(0,0,255,255))
    pos_s.addPoints(merged['x_pos_s'],merged['y_pos_s'])
    plt_pos_s.addItem(pos_s)


    win_1.nextRow()
    plt_x = win_1.addPlot(title="X coordinate on MCP")
    y,x = np.histogram(merged['x_pos'],bins=np.linspace(-20,20,120))
    plt_x.plot(x,y,stepMode=True,fillLevel=0,brush=(0,0,255,150))
    plt_y = win_1.addPlot(title="Y corrdinate on MCP")
    y,x = np.histogram(merged['y_pos'],bins=np.linspace(-20,20,120))
    plt_y.plot(x,y,stepMode=True,fillLevel=0,brush=(0,0,255,150))

    win_1.nextRow()
    plt_r = win_1.addPlot(title="R of ions on MCP")
    y,x = np.histogram(merged['radius'],bins=np.linspace(0,25,100))
    plt_r.plot(x,y,stepMode=True,fillLevel=0,brush=(0,0,255,150))
    plt_theta = win_1.addPlot(title="Theta of ions on MCP")
    y,x = np.histogram(merged['angle'],bins=np.linspace(0,360,180))
    plt_theta.plot(x,y,stepMode=True,fillLevel=0,brush=(0,0,255,150))
    lr = pg.LinearRegionItem([140,240],brush=(0,0,255,100))
    lr.setZValue(-10)
    plt_theta.addItem(lr)
    lr1 = pg.LinearRegionItem([50,140],brush=(0,255,0,100))
    lr1.setZValue(-10)
    plt_theta.addItem(lr1)
   
    spot_1_fit = pos_cut(merged,4,20,80,160)
    spot_2_fit = pos_cut(merged,4,20,220,340)
    print("Suggesting angle deviation after separated fitting=",spot_1_fit[2]-spot_2_fit[2])
    print(spot_1_fit)
    print(spot_2_fit)
    means_prior = []
    prior_weight =1
    if len(spot_1_fit)>0:
        if spot_1_fit[0]>=spot_2_fit[0]:
            prior_weight = spot_1_fit[0]/(spot_1_fit[0]+spot_2_fit[0])
            means_prior.append([spot_1_fit[2],spot_1_fit[5]])
            means_prior.append([spot_1_fit[5]*np.cos(spot_1_fit[2]*np.pi/180),spot_1_fit[5]*np.sin(spot_1_fit[2]*np.pi/180)])
        else:
            prior_weight = spot_2_fit[0]/(spot_1_fit[0]+spot_2_fit[0])
            means_prior.append([spot_2_fit[2],spot_2_fit[5]])
            means_prior.append([spot_2_fit[5]*np.cos(spot_2_fit[2]*np.pi/180),spot_2_fit[5]*np.sin(spot_2_fit[2]*np.pi/180)])
    print(prior_weight)
    print(means_prior)

    #pos_mix = merged[['x_pos_d','y_pos_d','angle_cor','radius_cor']]
    pos_mix = merged[['x_pos_d','y_pos_d','angle_cor','radius_cor']] #All data
    pos_mix = pos_mix.dropna()
    peak_x = pos_mix['x_pos_d'].values
    peak_y = pos_mix['y_pos_d'].values
    peak_data_angle = pos_mix['angle_cor'].values
    peak_data_r = pos_mix['radius_cor'].values
    pos_data_xy = np.stack((peak_x,peak_y),axis=-1)
    pos_data_polar = np.stack((peak_data_angle,peak_data_r),axis=-1)
    #print(pos_data_polar)
    pos_data_raw = np.hstack((pos_data_xy,pos_data_polar))
    n_ions = len(pos_data_raw)
   # print(pos_data_raw)
#    print("Total ions to be used in Gaussian mixture = ",n_ions))
    
    peak_counts,data_x, data_y = np.histogram2d(peak_x,peak_y,bins=np.arange(-20,20,1))
    y,angle_edges = np.histogram(peak_data_angle,bins = np.arange(0,360,2))
    y,radius_edges = np.histogram(peak_data_r,bins = np.arange(0,30,0.25))
    peak_counts, data_angle, data_r = np.histogram2d(peak_data_angle,peak_data_r,bins=(angle_edges,radius_edges))
    #print(data_angle)
    data_x = np.delete(data_x,-1)
    data_y = np.delete(data_y,-1)
    data_angle = np.delete(data_angle,-1)
    data_r = np.delete(data_r,-1)
    n_class = 3
    n_fittings = 10
    
    x_fit = np.linspace(-20,20,40)
    y_fit = np.linspace(-20,20,40)
    X_fit, Y_fit = np.meshgrid(x_fit,y_fit)
    XX = np.array([X_fit.ravel(),Y_fit.ravel()]).T

    score_r =-1e5
    params_xy_r = []
    Proba_raw = np.zeros(n_class)
    for i in range(n_fittings):
        clf = mixture.BayesianGaussianMixture(n_components=n_class,covariance_type='full',tol=1e-12,max_iter=1000,weight_concentration_prior=prior_weight,mean_prior= means_prior[1])
        Z_pre = clf.fit_predict(pos_data_xy)
        score_tem = clf.score(pos_data_xy)
        Z_proba_tem = clf.predict_proba(pos_data_xy)
        if clf.converged_:
            if score_tem >score_r:
                score_r = score_tem
                weights_r = clf.weights_
                means_r = clf.means_
                covariances_r = clf.covariances_
                Proba_r = Z_proba_tem
                score_sam = -clf.score_samples(XX)
        else:
            print("Original: One fitting is not converged")
    for j in range(n_class):
        if weights_r[j] == np.amax(weights_r):
            w_max = weights_r[j]
            p1 = j
        elif weights_r[j] == np.amin(weights_r):
            w_min = weights_r[j]
            p3 = j
        else:
            w_mid = weights_r[j]
            p2 = j
        for N in range(len(Proba_r)):
            Proba_raw[j] = Proba_raw[j] + Proba_r[N][j]
        print("Probability belong to spot%s"%j,Proba_raw[j])
        print("Weight belong to spot%s"%j,weights_r[j])
        print("Mean position belong to spot%s"%j,means_r[j])
        print("Covariance matrix belong to spot%s"%j,covariances_r[j])
    print(Proba_r)
    Proba_s = np.zeros(n_class)
    Proba_s_max = []
    Proba_var = np.zeros(n_class)
    score_b = -1e5
    re_samples_n = 10
    weight_var = np.zeros(n_class)
    mean_var = np.zeros((n_class,2))
    cov_var = np.zeros((n_class,2,2))
    print(len(pos_data_xy))
    for i in range(re_samples_n):
        re_sample = resample(pos_data_xy,n_samples=len(pos_data_xy))
    #    print(len(re_sample))
        for j in range(n_fittings):
            clf_s = mixture.BayesianGaussianMixture(n_components=n_class,covariance_type='full',tol=1e-12,max_iter=1000,weight_concentration_prior=prior_weight,mean_prior= means_prior[1])
            Z_sample = clf_s.fit_predict(re_sample)
            score_tem = clf_s.score(re_sample)
            Proba_tem = clf_s.predict_proba(re_sample)
            if clf_s.converged_:
                if score_tem >score_b:
                    score_b = score_tem
                    weights_b = clf_s.weights_
                    means_b = clf_s.means_
                    covariances_b = clf_s.covariances_
                    Proba_b = Proba_tem
            else:
                weights_b = weights_r
                means_b = means_r
                covariances_b = covariances_r
                Proba_b = Proba_r
                print("Sample: one fitting is not converged")
        for l in range(n_class):
            Proba_s[l] = 0
            for N in range(len(Proba_b)):
                Proba_s[l] = Proba_s[l] + Proba_b[N][l]
            if weights_b[l] == np.amax(weights_b):
                weight_var[0] = weight_var[0]+(weights_b[l]-w_max)**2
                Proba_var[0] = Proba_var[0] + (Proba_s[l]-Proba_raw[p1])**2
                Proba_s_max.append(Proba_s[l])
                for k in range(2):
                    mean_var[0][k] = mean_var[0][k] + (means_b[l][k]-means_r[p1][k])**2
                    for m in range(2):
                        cov_var[0][k][m] = cov_var[0][k][m] + (covariances_b[l][k][m]-covariances_r[p1][k][m])**2
            elif weights_b[l] == np.amin(weights_b):
                weight_var[2] = weight_var[2]+(weights_b[l]-w_min)**2
                Proba_var[2] = Proba_var[2] + (Proba_s[l]-Proba_raw[p3])**2
                for k in range(2):
                    mean_var[2][k] = mean_var[2][k] + (means_b[l][k]-means_r[p3][k])**2
                    for m in range(2):
                        cov_var[2][k][m] = cov_var[2][k][m] + (covariances_b[l][k][m]-covariances_r[p3][k][m])**2
            else:
                weight_var[1] = weight_var[1]+(weights_b[l]-w_mid)**2
                Proba_var[1] = Proba_var[1] + (Proba_s[l]-Proba_raw[p2])**2
                for k in range(2):
                    mean_var[1][k] = mean_var[1][k] + (means_b[l][k]-means_r[p2][k])**2
                    for m in range(2):
                        cov_var[1][k][m] = cov_var[1][k][m] + (covariances_b[l][k][m]-covariances_r[p2][k][m])**2
#            print("Proba  of each sample = ",Proba_s[l])
    for l in range(n_class):
        print("Variance for spot%s"%l,Proba_var[l]/re_samples_n)
        print("Sqrt Variance for spot%s"%l,np.sqrt(Proba_var[l]/re_samples_n))
        print("Variance for weight of spot%s"%l, weight_var[l]/re_samples_n)
        print("Sqrt Variance for weight of spot%s"%l, np.sqrt(weight_var[l]/re_samples_n))
    print("Mean of number of samplings",np.mean(Proba_s_max))
    print(Proba_b)
    Z_fit_0 = score_sam.reshape(X_fit.shape)
    Cs = plt.contour(X_fit,Y_fit,Z_fit_0,norm = colors.LogNorm(vmin=1.0,vmax= 100.0),levels=np.logspace(0,3,40))
    CB = plt.colorbar(Cs,shrink=0.8,extend='both')

    proba = np.zeros(n_class)
    counts_cor = np.zeros(n_class)
    uncer = np.zeros(n_class)
    theta_mixture = []
    theta_un =[]
    r_mixture = []
    relative_uncer = np.zeros(n_class)
    homo_counts =[]
    for j in range(n_class):
        #sigma_x_fit_error = np.sqrt(params_xy[2][j][0][0]/2/n_ions)
        #sigma_y_fit_error = np.sqrt(params_xy[2][j][1][1]/2/n_ions)
        x_fit_error = 2.355*np.sqrt(covariances_r[j][0][0])
        y_fit_error = 2.355*np.sqrt(covariances_r[j][1][1])
        theta_error = 180*(np.sqrt((y_fit_error*means_r[j][0]/(means_r[j][0]**2+means_r[j][1]**2))**2+(x_fit_error*means_r[j][1]/(means_r[j][0]**2+means_r[j][1]**2))**2))/np.pi
        theta_un.append(theta_error)
        theta_mixture.append(np.arctan2(means_r[j][1]-y_center,means_r[j][0]-x_center)*180/np.pi if (means_r[j][1] - y_center >=0) else np.arctan2(means_r[j][1]-y_center,means_r[j][0]-x_center)*180/np.pi + 360)
        r_mixture.append(np.sqrt( np.square(x_center-means_r[j][0])+np.square(y_center-means_r[j][1])))   
        for i in range(len(Proba_r)):
            #counts_cor[j] = counts_cor[j] + Proba[i][j]/homo_cut(pos_data_raw[i][2],pos_data_raw[i][3])
            #homo_counts.append(homo_cut(pos_data_raw[i][2],pos_data_raw[i][3]))
            proba[j] = proba[j]+Proba_r[i][j]
        uncer[j] = math.sqrt(proba[j])
        relative_uncer[j] = 1/math.sqrt(proba[j])
        #print("Standard deviation",np.std(homo_counts))
        print("XY: Probability belong to Spot %s ="%j, proba[j])
        #print("XY: After correction, probability belong to Spot %s ="%j, counts_cor[j])
        print("XY: Radius %s"%j,r_mixture[j])
        print("XY: Theta %s"%j,theta_mixture[j])
        print("XY: Uncertainty of Angle%s"%j,theta_un[j])
        print("Uncertainty of spot %s ="%j, uncer[j])
        print("Relative uncertainty of spot %s ="%j, relative_uncer[j])

    spot_color = Proba_r#np.insert(Proba_r,2,0,axis=1)
   # spot_color = np.delete(Z_proba,3,axis=1)
    plt.scatter(pos_data_xy[:,0],pos_data_xy[:,1],.8,c=spot_color)
    plt.xlabel("X position (mm)")
    plt.ylabel("Y position (mm)")
    #plt.savefig("129In.pdf")
 #   plt.subplot(222)
  #  spot_color = Z_proba_polar#np.insert(Z_proba_polar,2,0,axis=1)
    #spot_color = np.delete(Z_proba_polar,3,axis=1)
   # plt.scatter(pos_data_xy[:,0],pos_data_xy[:,1],.8,c=spot_color)
    #plt.xlabel("X position (mm)")

#plt.show()
#    theta_data = merged['angle'].dropna().tolist()
#    if len(theta_data)>0:
#        y_va,x_borders = np.histogram(theta_data,bins=np.arange(0,360,2))
#        x_borders = np.delete(x_borders,-1)
#        x_cen = [i + 1 for i in x_borders]        
#        y_errors = np.sqrt(y_va)
#        y_errors[y_errors ==0] = 100000
#        try:
#            popt,pcov = curve_fit(doublegauss,x_cen,y_va,p0=[spot_1_fit[1],spot_1_fit[2],spot_1_fit[3],spot_2_fit[1],spot_2_fit[2],spot_2_fit[3]])
#        except:
#            popt =[]
#            popt.append(0)
#    print("Parameters after double gauss fitting",popt)
#    print("Suggesting angle deviation after an uniform fitting=",popt[4]-popt[1])
#    theta_point = np.arange(1,360,2)
#    plt_theta.plot(theta_point,doublegauss(theta_point,popt[0],popt[1],popt[2],popt[3],popt[4],popt[5]),stepMode=False,fillLevel=0,pen=pg.mkPen(color=(255,0,0,150),width=2))
    #plt_3d_fit = gl.GLSurfacePlotItem(x=X_fit,y=Y_fit,z=Z_fit,shader='shaded',color=(.6,0,0,.8))
    plt.show()

def decay_correction(path_s,Y_is,Y_gs):
    half_live = np.array([0.57,1.23]) #s g.s. i.s.
    setting_file = path_s + "settings.ini"
    f_setting = open(setting_file)
    settings = f_settings.readlines()
    filling_f = re.findall(r'\d+',settings[537][1])
    filling_t = re.findall(r'\d+',settings[561][1])

def homo_cut(mid_theta,mid_radius):
    homo_data = load_homo_data("homo_sel.csv")
    bining = 5
    radius_min = mid_radius-0.5
    radius_max = mid_radius+0.5
    homo_data = homo_data[homo_data['radius']>= radius_min]
    homo_data = homo_data[homo_data['radius']< radius_max]
    if mid_theta + bining/2 > 360:
        n_eff = len(homo_data[homo_data['angle']>= mid_theta - bining/2])
        n_eff = len(homo_data[homo_data['angle']< mid_theta + bining/2-360])+n_eff
    elif mid_theta-bining/2 < 0:
        n_eff = len(homo_data[homo_data['angle']< mid_theta + bining/2])
        n_eff = len(homo_data[homo_data['angle']>= mid_theta - bining/2+360])+n_eff
    else:
        homo_data = homo_data[homo_data['angle']>=mid_theta-bining/2]
        homo_data = homo_data[homo_data['angle']<mid_theta+bining/2]
        n_eff = len(homo_data)
#        print("Okay?")

    return n_eff

def pos_cut(pos_data,r_min,r_max,the_min,the_max):
    if pos_data is None:
        return
    spot = pos_data[pos_data['radius']>r_min]
    spot = spot[spot['radius']<r_max]
    spot = spot[spot['angle']>the_min]
    spot = spot[spot['angle']<the_max]
    print("Number of ions in selected spot = ",len(spot))

    plt_pos_selected = win_1.addPlot(row=1,col=2,title="Selected ions")
    plt_pos_selected.setRange(xRange=(-20,20),yRange=(-20,20))
    plt_pos_selected.setLabel('left',text="Y position",units='mm')
    plt_pos_selected.setLabel('bottom',text="X position",units='mm')
    pos_selected = pg.ScatterPlotItem(size=2,pen=pg.mkPen(None),brush=pg.mkBrush(0,0,255,255))
    pos_selected.addPoints(spot['x_pos_d'],spot['y_pos_d'])
    plt_pos_selected.addItem(pos_selected)
    plt_the_selected = win_1.addPlot(row=2,col=2,title="Angles of selected ions on MCP")
    y,x = np.histogram(spot['angle'],bins=np.linspace(0,360,180))
    plt_the_selected.plot(x,y,stepMode=True,fillLevel=0,brush=(0,0,255,150))
    plt_r_selected = win_1.addPlot(row=0,col=3,title="Radius of selected ions on MCP")
    y,x = np.histogram(spot['radius'],bins=np.linspace(0,30,120))
    plt_r_selected.plot(x,y,stepMode=True,fillLevel=0,brush=(0,0,255,150))

    the_mean = spot['angle'].mean(skipna=True)
    the_sigma = spot['angle'].std(skipna=True)
    n_ions = spot['angle'].count()
    print("Number of fitting ions=",n_ions)
#    print("mean = ",the_mean)
#    print("sigma = ",the_sigma)
    std_of_mean = the_sigma/(np.sqrt(n_ions))
    mu_r= spot['radius'].mean(skipna=True)
    sigma_r = spot['radius'].std(skipna=True)
    
    the_position = spot['angle'].dropna().tolist()
#    print(x_position)
    if len(the_position)> 0:
        y,x_bin_borders = np.histogram(the_position,bins=np.arange(the_min,the_max,2))
        x_bin_borders = np.delete(x_bin_borders,-1)
        x_cent = [i + 0.5 for i in x_bin_borders]        
        y_errors = np.sqrt(y)
        y_errors[y_errors ==0] = 100000
        height = np.amax(y)
        x_inival = x_cent[np.argmax(y)]
        try:
            popt, pcov = curve_fit(helpers.gauss,x_cent,y,p0=[height,x_inival,the_sigma],sigma = y_errors)
        except:
            popt =[]
            popt.append(0)
            popt.append(0)
            popt.append(1e6)
        the_fit_plot_point = np.arange(the_min+1,the_max,2)

        plt_the_selected.plot(the_fit_plot_point,helpers.gauss(the_fit_plot_point,popt[0],popt[1],popt[2]),stepMode = False,fillLevel = 0,pen=pg.mkPen(color=(255,0,0,255),width =2))

#        print("Fitting position",popt[1])
#        print("Fitting width",popt[2])i
    r_position = spot['radius'].dropna().tolist()
    if len(r_position)> 0:
        r_y,r_x_bin_borders = np.histogram(r_position,bins=np.arange(r_min,r_max,0.25))
        r_x_bin_borders = np.delete(r_x_bin_borders,-1)
        r_x_cent = [i + 0.125 for i in r_x_bin_borders]        
        r_y_errors = np.sqrt(r_y)
        r_y_errors[r_y_errors ==0] = 100000
        r_height = np.amax(r_y)
        r_x_inival = r_x_cent[np.argmax(r_y)]
        try:
            r_popt, r_pcov = curve_fit(helpers.gauss,r_x_cent,r_y,p0=[r_height,r_x_inival,sigma_r],sigma = r_y_errors)
        except:
            r_popt =[]
            r_popt.append(0)
            r_popt.append(0)
            r_popt.append(1e6)
        r_fit_plot_point = np.arange(r_min+0.125,r_max,0.25)

        plt_r_selected.plot(r_fit_plot_point,helpers.gauss(r_fit_plot_point,r_popt[0],r_popt[1],r_popt[2]),stepMode = False,fillLevel = 0,pen=pg.mkPen(color=(255,0,0,255),width =2))
    #    print (len(spot),r_popt[1],r_popt[2])
        return [len(spot),popt[0],popt[1],popt[2],r_popt[0],r_popt[1],r_popt[2]]
    else:
        return [len(spot),0,0,0,0,0,0]

def xy_2_polar(data_xy):
    if len(data_xy) > 0:
        data_xy['angle'] = data_xy.apply(
                lambda row: np.arctan2(row.y_pos_d - y_center,row.x_pos_d-x_center)*180/np.pi if (row.y_pos_d - y_center >=0) else np.arctan2( row.y_pos_d-y_center, row.x_pos_d-x_center)*180/np.pi + 360,axis =1)
        data_xy['radius'] = np.sqrt( np.square(x_center-data_xy['x_pos_d'])+np.square(y_center-data_xy['y_pos_d']))
    return data_xy

def polar_2_xy(data_polar):
    if len(data_polar) > 0:
        data_polar['y_pos_d'] = data_polar.apply(
                lambda row: row.radius*np.sin(row.angle*np.pi/180))
        data_polar['x_pos_d'] = data_polar.apply(
                lambda row: row.radius*np.cos(row.angle*np.pi/180))
    return data_polar

def doublegauss(axis_data,amp_1,cent_1,width_1,amp_2,cent_2,width_2):
    y = []
    for x in axis_data:
        gauss_1 = amp_1*math.exp(-((x-cent_1)/width_1/2)**2)
        gauss_2 = amp_2*math.exp(-((x-cent_2)/width_2/2)**2)
        y.append(gauss_1+gauss_2)
    return y


def gauss2d(params_data):
    mu_x = params_data[0]
    mu_y = params_data[2]
    sigma_x = params_data[1]
    sigma_y = params_data[3]
    rho = params_data[4]
    amp = 20   # Should come from fitting
    offset = 0 # background
    x = np.arange(-20,20,0.5)
    y = np.arange(-20,20,0.5)
    cross_x =[]
    cross_y =[]
    g = np.zeros([len(x),len(x)])
    for i in range(len(x)):
        for j in range(len(y)):
            a = ((x[i]-mu_x)/sigma_x)**2
            b = ((y[j]-mu_y)/sigma_x)**2
            c = 2*rho*(x[i]-mu_x)*(y[j]-mu_y)/sigma_x/sigma_y
            g[i][j] = offset + amp*np.exp((a+b-c)/(rho**2-1)/2)
            if g[i][j] > amp/4:#-5 and g[i][j]<amp/2+5:
                cross_x.append(x[i])
                cross_y.append(y[j])
    print(len(cross_x))
    plt_gauss_3d = gl.GLSurfacePlotItem(x=x,y=y,z=g,shader='shaded',color=(0.5,0,0,0.5))
#    plt_3d = gl.GLSurfacePlotItem(x=data_x,y=data_y,z=peak_counts,shader='normalColor')
    plt_gauss_3d.translate(-10,10,0)
    win_2.addItem(plt_gauss_3d)
#    return g.ravel()


def fitting_hist(fit_data,fit_min,fit_max):
    
    ini_mean = fit_data.mean(skipna=True)/ps2ns
    ini_sigma = fit_data.std(skipna=True)/ps2ns
    print("Mean of center tof = ",ini_mean)
    print("Sigma of center tof = ",ini_sigma)
    ini_sigma = 1
    if len(fit_data) >0:
        y,x_bin_borders = np.histogram(fit_data/ps2ns,bins=np.linspace(fit_min,fit_max,fit_max-fit_min))
        x_bin_borders = np.delete(x_bin_borders,-1)
        x_cent = [i + 0.5 for i in x_bin_borders]
        y_errors = np.sqrt(y)
        y_errors[y_errors==0] =100000
        height = np.amax(y)
        x_inival = x_cent[np.argmax(y)]
        try:
            popt,pcov = curve_fit(helpers.gauss,x_cent,y,p0=[height,x_inival,ini_sigma],sigma=y_errors)
        except:
            popt =[]
            popt.append(0)
            popt.append(0)
            popt.append(1e6)
        x_fit_point = np.arange(fit_min,fit_max,1)
        #plt_fitting.plot(x_fit_ponit,helpers.gauss(x_fit_point,popt[0],popt[1],popt[2]),stepMode = False,fillLevel=0,pen=pg.mkPen(255,0,0,255),width=2)
        print("Fitting Center",popt[1])
        print("Fitting width",popt[2])
        center_pos=popt[1]

#    return center_pos

load_csv()
#ana_process(merged_raw)

if __name__ == '__main__':
    import sys
    if (sys.flags.interactive != 1) or not hasattr(QtCore, 'PYQT_VERSION'):
        QtGui.QApplication.instance().exec_()
