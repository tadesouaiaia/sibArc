#!/usr/bin/python3
import sys
from collections import defaultdict as dd
import numpy as np
from scipy import stats
import random
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib
from mpl_toolkits.axes_grid1.inset_locator import inset_axes      
matplotlib.rcParams['xtick.labelsize'] = 24                                                                                                                                                                                 
matplotlib.rcParams['ytick.labelsize'] = 24                                                                                                                                                                                 
#matplotlib.rcParams['xtick.major.size'] = 18                                                                                                                                                                               
#matplotlib.rcParams['ytick.major.size'] = 18                                                                                                                                                                               
matplotlib.rcParams['axes.linewidth'] = 1.85  


plt.rc('text', usetex=True)
mpl.use("Cairo")
import warnings
warnings.filterwarnings("ignore")


# printi

def rage_color_bar(hmap, ax, minv=0,maxv=99,LOW=False,MID=False):                                                                                                                                                                               
    CBAR = plt.cm.ScalarMappable(cmap=hmap, norm=plt.Normalize(vmin = minv, vmax=maxv))                                                                                                                                     
    CBAR._A = []                                                                                                                                                                                                            
    if LOW: axins = inset_axes(ax, width = "200%", height="35%", loc = 'upper left',bbox_to_anchor=(0.66,-1.0,0.4,0.3),bbox_transform=ax.transAxes,borderpad=0,)                                                                    
    elif MID:   axins = inset_axes(ax, width = "200%", height="35%", loc = 'upper left',bbox_to_anchor=(0.09,-0.80,0.4,0.3),bbox_transform=ax.transAxes,borderpad=0,)                                                                    
    else:   axins = inset_axes(ax, width = "200%", height="35%", loc = 'upper left',bbox_to_anchor=(0.10,-0.6,0.4,0.3),bbox_transform=ax.transAxes,borderpad=0,)                                                                    
    plt.colorbar(CBAR, cax = axins, orientation='horizontal', ticks=[])                                                                                                                                                     
    yp = -22                                                                     
    if LOW: 
        ax.text(85,yp-42,'Polygenic Heritability ($h_p^2$)',clip_on=False, fontsize=25)                                                                                                                                          
        ax.text(65,yp-42,'0\%',clip_on=False, fontsize=22)                                                                                                                                                                         
        ax.text(140,yp-41,'100\%',clip_on=False, fontsize=22)     
    elif MID: 
        ax.text(35,yp-23,'Polygenic Heritability ($h_p^2$)',clip_on=False, fontsize=21)                                                                                                                                          
        ax.text(10,yp-26,'0\%',clip_on=False, fontsize=20)                                                                                                                                                                         
        ax.text(85,yp-26,'100\%',clip_on=False, fontsize=20)     
    else: 
        ax.text(28,yp-2,'Polygenic Heritability ($h_p^2$)',clip_on=False, fontsize=25)                                                                                                                                          
        ax.text(10,yp,'0\%',clip_on=False, fontsize=22)                                                                                                                                                                         
        ax.text(83,yp,'100\%',clip_on=False, fontsize=22)     



class SibTable:
    def __init__(self, args, ax,  clr='white', PETT=False): 
        self.args, self.ax = args, ax 
        self.rows, self.cols = [], [] 

    def get_loc(self,X,Y):
        x1,x2 = X[0]/100.0 , X[1]/100.0
        y1,y2 = Y[0]/100.0 , Y[1]/100.0
        return (x1,y1,x2-x1,y2-y1)

    def add_row(self,row_data,X=None,Y=None,COLORS=[],WIDTHS=[],FS=13,ALPHA=0,TITLE=False, CL = 'center',CLEAR=False): 
        if X == None: X = self.X_SPAN
        if Y == None: Y = self.Y_SPAN
        cL,rL,rD,xL = CL,None,[row_data],len(row_data)  
        bl = self.get_loc(X,Y) 
        while len(WIDTHS) < len(row_data): WIDTHS.append(10) 
        while len(COLORS) < len(row_data): COLORS.append('white') 
        if CL != 'center': row = self.ax.table(cellText=rD,rowLabels=rL,cellColours = [COLORS[0:xL]],colWidths=WIDTHS[0:xL],bbox=bl,loc = cL, cellLoc=cL, alpha=ALPHA, clip_on=False) 
        else:              row = self.ax.table(cellText=rD,rowLabels=rL,cellColours = [COLORS[0:xL]],colWidths=WIDTHS[0:xL],bbox=bl,cellLoc=cL, alpha=ALPHA, clip_on=False) 
        row.auto_set_font_size(False)
        row.set_fontsize(FS)
        table_props = row.properties()
        if ALPHA > 0: 
            for cell in row._cells: row._cells[cell].set_alpha(ALPHA)
        self.rows.append(row) 


    def make_res(self, tt, mode, flen): 
        
        if   mode == 1:   fs0,fs1,fs2,fs3 = 16,18,16,14 
        elif mode == 2:   fs0,fs1,fs2,fs3 = 17,14,13,11 
        elif mode == 3:   fs0,fs1,fs2,fs3 = 15,13,11,9 
        else:             fs0,fs1,fs2,fs3 = 16,14,12,10 
        
        c1,c2 = 'lightgrey','whitesmoke'
        SW = [18,14,14]
        th = 9 
        self.add_row(['Lower 1\%','Tail Result','Upper 1\%'],COLORS=[c1,c1,c1],X=(0,100),Y=(88,100), FS = fs0, WIDTHS=[44,12,44], TITLE=True) 
        
        MM = [['Summary'], ['De Novo'], ['Mendelian'], ['Optimal\nRates']]
        ST = [['Samples','Index Mean','Dist-P-val'],['Expected Mean($s_2$)','Observed','P-val'],['Exp Count($s_2 \in$ 1\%)','Observed','P-val'],['De Novo','Mendelian','Polygenic']]
        yL,yS = 88, 22
        for M,S in zip(MM,ST): 
            self.add_row(M,X=(44,56),Y=(yL-yS,yL), WIDTHS=[12], COLORS=[c2],FS=fs1,TITLE=True) 
            self.add_row(S,X=(0,44),Y=(yL-th,yL), WIDTHS=SW, COLORS=[c2,c2,c2],FS=fs2-1,TITLE=True) 
            self.add_row(S,X=(56,100),Y=(yL-th,yL), WIDTHS=SW, COLORS=[c2,c2,c2],FS=fs2-1,TITLE=True) 
            yL -= yS
        rd = [tt['0'],tt['99']] 
        D1 = [[r.size,round(r.index_avg,2), r.d_str] for r in rd] 
        D2 = [[round(r.n_exp,3), round(r.n_obs,3), r.n_str] for r in rd] 
        D3 = [[round(r.m_exp,2), int(r.m_obs), r.m_str] for r in rd] 
        D4 = [[r.rates['novo'],r.rates['mend'],r.rates['poly']] for r in rd] 
        yL,yS = 88, 22
        for D in [D1,D2,D3,D4]: 
            self.add_row(D[0],X=(0,44),Y=(yL-yS,yL-th), WIDTHS=SW, FS=fs2+1,TITLE=False) 
            self.add_row(D[1],X=(56,100),Y=(yL-yS,yL-th), WIDTHS=SW, FS=fs2+1,TITLE=False) 
            yL-=yS 
        ek = 'darkslategray'
        self.ax.plot([0,100],[101,101],color=ek,linewidth=3,clip_on=False)
        self.ax.plot([0,100],[0,0],color=ek,linewidth=3,clip_on=False)
        self.ax.plot([0,0],[0,101],color=ek,linewidth=3,clip_on=False)
        self.ax.plot([100,100],[0,101],color=ek,linewidth=3,clip_on=False)
        self.ax.set_ylim(0,100) 
        self.ax.set_xlim(0,100) 
        self.ax.axis('off') 
        return self 
    


    def make_h2(self, h2, mode, flen):


        if   mode == 1:   fs0,fs1,fs2,fs3 = 16,15,13.5,12 
        elif mode == 2:   fs0,fs1,fs2,fs3 = 17,14,13,11 
        elif mode == 3:   fs0,fs1,fs2,fs3 = 14,13,11,9 
        else:             fs0,fs1,fs2,fs3 = 16.5,13,11.5,10 
        c1,c2 = 'lightgrey','whitesmoke'
        self.add_row(['Conditional\nHeritability'],COLORS=[c1,c2,c1],X=(0,100),Y=(90,100), FS = fs0, WIDTHS=[100], TITLE=True) 
        self.add_row(['Range\nCovered\n[Samples]','Iterative\nEstimate\n'+'(95\% CI)'],COLORS=[c2,c2],X=(0,100),Y=(77,90), FS = fs1, WIDTHS=[50,50], TITLE=True)  
        h_names, h_locs =    ['All','Bod','Mid','LoH','HiH','LoT','HiT'], ['0-100','5-95','35-65','5-40','60-95','0-4','96-100']
        h_smart = ['Full','Dist Body','Middle','Lower','Upper','Low Tail','Top Tail'] 
        yL,yS = 77, 11
        for i,(k,l,s) in enumerate(zip(h_names, h_locs, h_smart)):
            if i % 2 != 0: myclr = c2 
            else:          myclr = 'white' 
            hd = h2.CI[k] 
            hz = str(h2.PW[k][0])
            h_str = str(hd[0])+'\n('+str(hd[1])+','+str(hd[2])+')'
            self.add_row([s+':\n'+l+'\n['+hz+']',h_str],COLORS=[myclr,myclr],X=(0,100),Y=(yL-yS,yL), FS = fs2, WIDTHS=[50,50], TITLE=True) 
            yL -= yS 
        self.ax.axis('off') 
        return




#############################################################################################
##################################  PROGRESS STREAM  ########################################
#############################################################################################


class SibProgress:
    def __init__(self, args, command_line): 
        self.args = args
        if args.silent: self.ACTIVE = False 
        else:           self.ACTIVE = True 
        #self.out1  = open(self.args.out+'.progress.log','w') 
        self.out2  = sys.stderr 
        self.space = '' 
        self.show('\nSibArc Begins:  '+command_line+'\n')
        self.show('   Input Files: '+",".join([sf.name.split('/')[-1] for sf in args.sibfiles])+'\n') 
        self.show('Output Prefix: '+self.args.out+'\n\n')
        self.loc = None


    def initialize(self, f_name, t): 
        self.show('Beginning: '+f_name+'\n') 




    def show(self, msg, space=''):
        
        if space == 'NA': myspace = '' 
        else:             myspace = self.space

        #self.out1.write(myspace+msg) 
        if self.ACTIVE: 
            self.out2.write(myspace+msg)
            self.out2.flush() 

    def start(self,msg):
        if self.loc is not None: self.show('...Finished\n','NA') 
        self.space = '' 
        self.show('\n'+msg+':\n') 
        self.loc = None 
        self.space = '  '

    def end(self, msg): 
        if self.loc is not None: self.show('...Finished ('+msg+')\n',space='NA') 
        self.loc = None 
        return

    def finish(self): 
        if self.loc is not None: self.show('...Finished\n','NA') 
        self.show('\nSibArc Completed\n','NA') 

    def update(self, msg): 
        if self.loc is not None: self.show('...Finished\n','NA') 
        self.loc = msg 
        self.show(msg+'...') 
     
























#############################################################################################
##################################  TRAIT PLOT     ##########################################
#############################################################################################



class SibPlot:
    def __init__(self, args, progress):
        self.args, self.progress = args, progress 
        self.fig_prefix = self.args.out+'.fig' 
        
        if self.args.savePlotdata:
            self.out = open(self.args.out+'.plotData','w') 
            self.out.write('%-30s %40s %15s %s\n' % ('---', 'name', 'dataType', 'values')) 
        self.names = args.names 
        self.color_rainbow = plt.get_cmap('coolwarm')
        self.herits = [round(i * 0.05, 2) for i in range(21)]
        self.h_colors = [self.color_rainbow(1.*i/float(len(self.herits))) for i in range(len(self.herits))]
        self.color_key = {h: c for h, c in zip(self.herits, self.h_colors)}
        self.rounds = 1 
        self.flen = min(6, len(self.names)) 
        self.setup() 

    def setup(self):
        self.fig = mpl.pyplot.gcf()
        self.axes, self.ax_index = [], 0  
        xs1,xs2,xs3 = 16, 4,7
        ys1,ys2,ys3 = 16, 4,4
        if   self.flen == 1:   self.mode, self.WD, self.HT, self.rows, self.cols = 1, 12, 11, 25, 22 
        elif self.flen == 2: self.mode, self.WD, self.HT, self.rows, self.cols = 2, 19, 10, 25, 44
        elif self.flen <= 4: self.mode, self.WD, self.HT, self.rows, self.cols = 3, 19, 15, 53, 44 
        else:           self.mode, self.WD, self.HT, self.rows, self.cols = 4, 26.5, 18, 53, 68
        self.fig.set_size_inches(self.WD, self.HT)
        if self.mode == 1:   self.fs0,self.fs1, self.fs2, self.fs3 = 16,18, 16, 14 
        elif self.mode == 2: self.fs0,self.fs1, self.fs2, self.fs3 = 17,14, 13, 11 
        elif self.mode == 3: self.fs0,self.fs1, self.fs2, self.fs3 = 15,14, 11.5, 9 
        else:                self.fs0,self.fs1, self.fs2, self.fs3 = 15,14, 12, 10 
        

        cnt, row, col = 0, 0, 1 
        while cnt < self.flen: 
            self.axes.append(plt.subplot2grid((self.rows, self.cols), (row,   col),     rowspan=xs1, colspan=ys1))
            self.axes.append(plt.subplot2grid((self.rows, self.cols), (row+xs1+2, col-1), rowspan=xs3, colspan=ys1+ys2+1))
            self.axes.append(plt.subplot2grid((self.rows, self.cols), (row, col+ys1),   rowspan=xs1, colspan=ys2))
            if self.mode == 1: return 
            elif self.mode > 1: 
                self.axes.append(plt.subplot2grid((self.rows, self.cols), (row,   col+(ys1+ys2+3)), rowspan=xs1, colspan=ys1))
                self.axes.append(plt.subplot2grid((self.rows, self.cols), (row+xs1+2, col+(ys1+ys2+3)-1), rowspan=xs3, colspan=ys1+ys2+1))
                self.axes.append(plt.subplot2grid((self.rows, self.cols), (row, col+(ys1+ys1+ys2+2)+1), rowspan=xs1, colspan=ys2))
                
                if self.mode < 4: cnt+= 2 
                else: 
                    self.axes.append(plt.subplot2grid((self.rows, self.cols), (row,   col+(ys1+ys1+ys2+ys2+6)), rowspan=xs1, colspan=ys1))
                    self.axes.append(plt.subplot2grid((self.rows, self.cols), (row+xs1+2, col+ys1+ys1+ys2+ys2+6-1), rowspan=xs3, colspan=ys1+ys2+1))
                    self.axes.append(plt.subplot2grid((self.rows, self.cols), (row, col+(ys1+ys1+ys1+ys2+ys2+ys2+2)), rowspan=xs1, colspan=ys2))
                    cnt += 3 
                row += (xs1+xs2+3+6) 

        


    def draw_summary(self, f_name, name, pairs, h2, tt): 
        self.pairs, self.h2, self.tt = pairs, h2, tt 
        
        self.yMax, self.X = 0.5, sorted(self.pairs.keys())  
        self.m1 = [np.mean([p[0] for p in self.pairs[x]]) for x in self.X]
        self.m2 = [np.mean([p[1] for p in self.pairs[x]]) for x in self.X]
        self.s2 = [stats.sem([p[1] for p in self.pairs[x]]) for x in self.X]
        
        self.ax, ax2, ax3 = self.axes[self.ax_index], self.axes[self.ax_index + 1], self.axes[self.ax_index+2] 
        
        self.draw_herit_lines() 
        self.draw_curves(name, self.h2.body) 
        self.ax.text(-1, self.yMax*0.98, " ".join(name.split('_')), va='top', ha='left', fontsize=self.fs0+12) 
        if self.args.savePlotdata: 
            self.out.write('%-30s %40s %15s %s\n' % (f_name, name, 'X', ','.join([str(x) for x in self.X]))) 
            self.out.write('%-30s %40s %15s %s\n' % (f_name, name, 'm1', ','.join([str(round(x,4)) for x in self.m1]))) 
            self.out.write('%-30s %40s %15s %s\n' % (f_name, name, 'm2', ','.join([str(round(x,4)) for x in self.m2]))) 
            self.out.write('%-30s %40s %15s %s\n' % (f_name, name, 'sem', ','.join([str(round(x,4)) for x in self.s2]))) 
        
        t1 = SibTable(self.args, self.axes[self.ax_index+1]).make_res(tt,  self.mode, self.flen)  
        t2 = SibTable(self.args, self.axes[self.ax_index+2]).make_h2(h2,   self.mode, self.flen)  
        self.ax_index += 3


    def draw_herit_lines(self): 
        for i, h in enumerate(self.herits):
            sh = [(s * h)/2.0 for s in self.m1]
            if i == 0: self.ax.plot(self.X, sh, color=self.color_key[h], linewidth=3, alpha=0.7)
            else:
                z1 = [(a+b+b)/3.0 for a, b in zip(sh, lp)]
                z2 = [(a+b)/2.0 for a, b in zip(sh, lp)]
                self.ax.plot(self.X, z1, color=self.color_key[h], linewidth=4, alpha=0.3)
                self.ax.plot(self.X, z2, color=self.color_key[h], linewidth=4, alpha=0.3)
                self.ax.plot(self.X, sh, color=self.color_key[h], linewidth=4, alpha=0.3)
            lp = [z for z in sh]
            if max(lp) > self.yMax: self.yMax = max(lp) 



    def draw_curves(self, name, bodyH):
        sE = [(s*bodyH)/2.0 for s in self.m1]
        self.ax.plot(self.X, sE, color='k', linewidth=2, alpha=1.0, zorder=2)
        for i, (x, y) in enumerate(zip(self.X[1:-1], self.m2[1:-1])): self.ax.scatter(x, y, s=25, alpha=0.75,marker='o', color='grey', edgecolor='k', zorder=50)
        

        self.yMax = max([self.yMax,abs(self.tt['0'].n_obs), abs(self.tt['99'].n_obs)]) + 0.1 
        


        for loc,mark in zip(['0','99'],['v','^']):            
            tt = self.tt[loc] 
            self.ax.scatter(int(loc), tt.n_exp, edgecolor='k', marker=mark, zorder=90,facecolor='whitesmoke', linewidth=1.5, s=200, clip_on=False)
            self.ax.scatter(int(loc), tt.n_obs, edgecolor='k', marker='h', zorder=100,facecolor=tt.clr, linewidth=1.2, s=200, clip_on=False)
        self.ax.set_yticks([-2, -1, 0, 1, 2])
        self.ax.set_xticks([0,20,80,100])
        self.ax.set_xlim(-1.5, 100.5)
        self.ax.set_ylim(-self.yMax, self.yMax)
        self.ax.set_xlabel('Index Sibling Rank \% ($s_{(1)}$)', fontsize=self.fs1, labelpad=-10)
        self.ax.set_ylabel('Conditional Sibling Z-Value ($s_2$)', fontsize=self.fs1, labelpad=-2)
        return


    def complete_plot(self): 
        fs1, fs2 = 20, 18 
        labs = ['Expected Mean\n(Lower Tail)','Expected Mean\nUpper Tail','Observed\nMean'] 
        marks = ['v','^','h'] 
        if self.mode < 3: 
            ax, ax2 = self.axes[1], self.axes[1] 
            yMin,yMax = ax.get_ylim() 
            ys = (yMax - yMin) / 10.0
            if self.mode == 1:   
                plt.subplots_adjust(left=0.02, bottom=0.15, right=1.036,top=0.985, wspace=0.15, hspace=0.02) 
                fs1, fs2 = 16, 13 
                xp,yp = 0.5, yMin - 2*ys 
                xJ = 19                                                                                                                                                                
                for lab,mark,step in zip(labs, marks, [2.5,2.5,2.5]): 
                    ax.scatter(xp,yp,marker=mark,s=233,clip_on=False,edgecolor='k', zorder=100, facecolor='whitesmoke',linewidth=2.5) #,clip_on=False)                                                                                   
                    ax.text(xp+step,yp,lab,fontsize=fs1,va='center')                                                                                                                                            
                    xp+= xJ 
                xp = 55                                                                                                                                                                                                  
                for c,n in zip(['gray','lime','red'],['Polygenic\nArchitecture','De Novo\nArchitecture','Mendelian\nArchitecture']):                                                                                                
                    if c != 'black':                                                                                                                                                                                                
                        ax.scatter(xp,yp,marker='o',s=333, color=c,clip_on=False)                                                                                                                                                   
                    ax.text(xp+2,yp,n,ha='left',va='center',fontsize=fs1,fontweight='bold')                                                                                                                                          
                    xp += xJ*0.90                                                                                                                                                                                             
                rage_color_bar('coolwarm', ax2, minv=0,maxv=99,MID=True) 
                                                                                                                      
            elif self.mode == 2: 
                plt.subplots_adjust(left=0.02, bottom=0.2, right=0.97,top=0.95, wspace=0.15, hspace=0.02) 
                xp,yp = 10, yMin - 2*ys 
                xJ = 30                                                                                                                                                                                                      
                for lab,mark,step in zip(labs, marks, [4.5,5.5,6]): 
                    ax.scatter(xp,yp,marker=mark,s=500,clip_on=False,edgecolor='k', zorder=100, facecolor='whitesmoke',linewidth=2.5) #,clip_on=False)                                                                                   
                    ax.text(xp+step,yp,lab,fontsize=fs1,va='center')                                                                                                                                            
                    xp+= xJ 
                xp = 125                                                                                                                                                                                                  
                for c,n in zip(['gray','lime','red'],['Polygenic\nArchitecture','De Novo\nArchitecture','Mendelian\nArchitecture']):                                                                                                
                    if c != 'black':                                                                                                                                                                                                
                        ax.scatter(xp,yp,marker='o',s=700, color=c,clip_on=False)                                                                                                                                                   
                    ax.text(xp+5,yp,n,ha='left',va='center',fontsize=24,fontweight='bold')                                                                                                                                          
                    xp += xJ*0.90                                                                                                                                                                                             
                rage_color_bar('coolwarm', ax2, minv=0,maxv=99, LOW=True)                                                                                                                                                            
        elif self.mode == 3: 
            plt.subplots_adjust(left=0.02, bottom=0.01, right=0.99,top=0.98, wspace=0.15, hspace=0.02) 
            ax, ax2 = self.axes[4], self.axes[4] 
            yMin,yMax = ax.get_ylim() 
            ys = (yMax - yMin) / 10.0 
            xp,yp = 10, yMin - (12*ys)  
            xJ = 30 
            for lab,mark,step in zip(labs, marks, [4.5,5.5,6]): 
                ax.scatter(xp,yp,marker=mark,s=500,clip_on=False,edgecolor='k', zorder=100, facecolor='whitesmoke',linewidth=2.5) #,clip_on=False)                                                                                   
                ax.text(xp+step,yp,lab,fontsize=fs1,va='center')                                                                                                                                            
                xp+= xJ 
            xp = 10   
            yp -= ys*10
            for c,n in zip(['gray','lime','red'],['Polygenic\nArchitecture','De Novo\nArchitecture','Mendelian\nArchitecture']):                                                                                                
                if c != 'black':                                                                                                                                                                                                
                    ax.scatter(xp,yp,marker='o',s=700, color=c,clip_on=False)                                                                                                                                                   
                ax.text(xp+5,yp,n,ha='left',va='center',fontsize=24,fontweight='bold')                                                                                                                                          
                xp += xJ*0.90                                                                                                                                                                                             
            rage_color_bar('coolwarm', ax2, minv=0,maxv=99)                                                                                                                                                            
        else:                
            plt.subplots_adjust(left=0.015, bottom=0.01, right=0.999,top=0.98, wspace=0.12, hspace=-0.2) 
            ax, ax2 = self.axes[9], self.axes[4] 
            yMin,yMax = ax.get_ylim() 
            ys = (yMax - yMin) / 10.0 
            xp,yp = 5, yMax + ys 
            xJ = 46                                                                                                                                                                                                        
            for lab,mark,step in zip(labs, marks, [4.5,5.5,6]): 
                ax.scatter(xp,yp,marker=mark,s=500,clip_on=False,edgecolor='k', zorder=100, facecolor='whitesmoke',linewidth=2.5) #,clip_on=False)                                                                                   
                ax.text(xp+step,yp,lab,fontsize=fs1,va='center')                                                                                                                                            
                xp+= xJ 
            xp = 300                                                                                                                                                                                                     
            for c,n in zip(['gray','lime','red'],['Polygenic\nArchitecture','De Novo\nArchitecture','Mendelian\nArchitecture']):                                                                                                
                if c != 'black':                                                                                                                                                                                                
                    ax.scatter(xp,yp,marker='o',s=700, color=c,clip_on=False)                                                                                                                                                   
                ax.text(xp+5,yp,n,ha='left',va='center',fontsize=24,fontweight='bold')                                                                                                                                          
                xp += xJ*0.90                                                                                                                                                                                             
            rage_color_bar('coolwarm', ax2, minv=0,maxv=99)                                                                                                                                                            
                                                                                                                  



        

    def reset(self,flen): 
        self.complete_plot()  
        plt.savefig(self.fig_prefix+str(self.rounds)+'.png', dpi=300) 
        plt.savefig(self.fig_prefix+str(self.rounds)+'.pdf', dpi=300) 
        plt.clf() 
        self.rounds+=1 
        self.flen = flen 
        self.setup() 
        


    def finish(self): 
        if self.flen in [3,5]: 
            for j in [-1,-2,-3]: self.axes[j].axis('off') 
        self.complete_plot()  
        if self.rounds > 1: fig_name = self.fig_prefix+str(self.rounds) 
        else:               fig_name = self.fig_prefix 
        plt.savefig(fig_name+'.png', dpi=300) 
        plt.savefig(fig_name+'.pdf', dpi=300) 

    


#############################################################################################
##################################  FAMILY DATA  ############################################
#############################################################################################


class SibData:
    def __init__(self, file_handle, args):
        self.args = args 
        self.fams, self.members = [], []
        self.read_file(file_handle)
        self.collate()
        #self.summarize() 

    def read_file(self, file_handle):
        first_line = file_handle.readline()
        if len(first_line.split()) == 2:      
            SPLIT='WHITESPACE' 
            line = first_line.split() 
        elif len(first_line.split(',')) == 2: 
            SPLIT='COMMA'
            line = first_line.split(',') 
        else:
            sys.stderr.write('Incorrect File Format\nUnrecognized Line: '+first_line) 
            sys.exit() 
        try:     
            k1, k2 = Kid(float(line[0]), str(float(i))), Kid(float(line[1]), str(i+0.5))
            self.fams.append([k1, k2])
            self.members.extend([k1, k2])
        except: pass 

        for i,line in enumerate(file_handle): 
            if     SPLIT == 'COMMA': line = line.split(',') 
            else:                    line = line.split() 
            try:               k1, k2 = Kid(float(line[0]), str(float(i))), Kid(float(line[1]), str(i+0.5))
            except ValueError: raise Exception('Incorrect File Format')
            self.fams.append([k1, k2])
            self.members.extend([k1, k2])
        self.total = len(self.members)
        self.size = len(self.fams)
        return

    def collate(self):
        self.members.sort(key=lambda X: X.val)
        for i, k in enumerate(self.members):
            k.z = stats.norm.ppf((i+0.5)/self.total)
            k.rank = int(100*(i+0.5) / self.total)

        self.pairs = dd(lambda: dd(list))
        for k1, k2 in self.fams:
            kids = [k1, k2]
            if self.args.randomize: random.shuffle(kids)
            rank = kids[0].rank
            self.pairs['vals'][kids[0].rank].append([kids[0].val, kids[1].val])
            self.pairs['norm'][kids[0].rank].append([kids[0].z, kids[1].z])
        return 

    def summarize(self): 
        self.data = dd(lambda: dd(list)) 
        for k, p in self.pairs.items():
            S1,S2 = [], []
            for pt in p.keys():
                if len(p[pt]) > 1: 
                    S1 = [DZ[0] for DZ in p[pt]]
                    S2 = [DZ[1] for DZ in p[pt]]
                    m1, m2, sm  = np.mean(S1), np.mean(S2), stats.sem(S2) 
                else:
                    m1,m2,sm = stats.norm.ppf((pt+0.05)/100.0), stats.norm.ppf((pt+0.05)/100.0),0 
                self.data[k][pt] = [m1, m2, sm] 
        return



class Kid:
    def __init__(self, val, idx):
        self.val, self.idx = val, idx
        self.z, self.rank = 'NA', 'NA'

    def __str__(self):
        return 'Individual-'+self.idx

    def __repr__(self):
        return 'Individual-'+self.idx


#############################################################################################
################################## TAIL TESTING   ###########################################
#############################################################################################


class TailTests:
    def __init__(self, args, pairs, h2):
        self.args, self.pairs, self.h2 = args, pairs, h2

    def get_data(self, pt, locs=[], ALL=False):
        my_data = []
        for loc in locs:
            my_data.extend(self.pairs[loc])
        if len(my_data) > 1:
            return my_data
        p = 1
        while len(my_data) < 2:
            if pt == 0:
                my_data.extend(self.pairs[locs[-1]+p])
            else:
                my_data.extend(self.pairs[locs[0]-p])
            p += 1
        return my_data

    def summarize(self):        
        self.raw_result = self.calculate() 
        return self




    def calculate(self):
        lower_tail = [[0],   [0, 1], [0, 1, 2],  [0, 1, 3],  [0, 1, 2, 3, 4], [0, 1, 2, 3, 4, 5]]
        upper_tail = [[99], [98, 99], [97, 98, 99], [96, 97, 98, 99], [95, 96, 97, 98, 99], [94, 95, 96, 97, 98, 99]]
        tail_locs = [0 for lt in lower_tail] + [1 for ut in upper_tail]
        my_var = 1 - (self.h2*self.h2)/4.0
        my_std = my_var ** 0.5
        
        self.results = {} 

        ###
        #print()
        #print() 
        my_data = sorted(self.get_data(0,[0]))
        index_sib = sorted([md[0] for md in my_data])[-1]
        #print(index_sib) 
        #mendTest  = TailTest(0,index_sib).run_mend(my_data, self.h2, my_var, my_std)
        ### 

        #sys.exit()  

        for vals, loc in zip(lower_tail+upper_tail, tail_locs):
            vs = str(vals[0])
            if vals[-1] != vals[0]:
                vs += '-'+str(vals[-1])
            my_data = sorted(self.get_data(loc, vals)) 
            index_sibs = sorted([md[0] for md in my_data]) 
            if loc == 0: index_sib = index_sibs[-1] 
            else:        index_sib = index_sibs[0] 
             

            novoTest  = TailTest(loc,index_sib).run_novo(my_data, self.h2, my_var, my_std)
            mendTest  = TailTest(loc,index_sib).run_mend(my_data, self.h2, my_var, my_std)
            distTest  = TailTest(loc,index_sib).run_dist(my_data, self.h2, my_var, my_std)

            self.results[vs] = TailResult(args, vs, novoTest, mendTest, distTest)  

        return self





class TailResult:
    def __init__(self, args, loc, nv, md, dt): 
        self.args, self.loc, self.size, self.index_avg = args, loc, nv.size, nv.index_sib
        self.rates = {'novo': dt.n_rate, 'mend': dt.m_rate, 'poly': dt.p_rate} 
        self.n_pv, self.m_pv, self.d_pv = nv.pv, md.pv, dt.pv 
        self.n_exp, self.n_obs = nv.exp, nv.obs
        self.m_exp, self.m_obs = md.exp, md.obs
        if self.d_pv < 0.0001: self.d_str = '%9.2e' % self.d_pv 
        else:                  self.d_str = '%9.4f' % self.d_pv 
        if self.n_pv < 0.0001: self.n_str = '%9.2e' % self.n_pv 
        else:                  self.n_str = '%9.4f' % self.n_pv 
        if self.m_pv < 0.0001: self.m_str = '%9.2e' % self.m_pv 
        else:                  self.m_str = '%9.4f' % self.m_pv 

        pvs = sorted([[nv.pv, 'novo'],[md.pv, 'mend'],[dt.pv, 'dist']])  

        if pvs[0][0]   > self.args.alpha:                            self.arch = 'polygenic'
        elif pvs[0][1] == 'dist': 
            if pvs[1][0] > self.args.alpha:    self.arch = 'undetermined' 
            elif pvs[1][1] == 'novo': self.arch = 'undetermined,denovo'
            elif pvs[1][1] == 'mend': self.arch = 'undetermined,mendelian'
        elif pvs[0][1] == 'novo':     self.arch = 'denovo'
        elif pvs[0][1] == 'mend':     self.arch = 'mendelian'
        else:                         self.arch = 'polygenic' 
        if self.arch == 'polygenic': self.clr = 'grey' 
        elif self.arch == 'denovo': self.clr = 'lime' 
        elif self.arch == 'mendelian': self.clr = 'red' 
        elif self.arch.split(',')[-1] == 'denovo': self.clr = 'green'
        elif self.arch.split(',')[-1] == 'mendelian': self.clr = 'maroon'
        else: self.clr = 'purple' 
        return









class TailTest:
    def __init__(self, loc, index_sib):
        if loc == 0: self.side = 'lower'
        else:        self.side = 'upper'
        self.index_sib = round(index_sib,5) 
        self.n_rate, self.m_rate, self.p_rate = 'NA', 'NA', 'NA' 



    def run_dist(self, tail_data, h2, h_var, h_std): 
        self.size = str(len(tail_data))
        S1, S2 = [s[0] for s in tail_data], [s[1] for s in tail_data] 
        expected = [(s1*0.5*h2) for s1,s2 in tail_data] 
        null_standardized = [(s2 - (s1*0.5*h2))/h_std for s1, s2 in tail_data] 
        self.Z, self.pv = stats.kstest(null_standardized, "norm") 
        self.exp = round(np.mean(expected),4) 
        self.obs = np.mean(S2) 
        tM = TailMax(self.side, tail_data, h2, h_var, h_std) 
        self.rates = tM.brute_force() 
        self.n_rate, self.m_rate, self.p_rate = round(self.rates[0],2), round(self.rates[1],2), round(self.rates[2],2) 
        return self 

    def run_novo(self, tail_data, h2, h_var, h_std):
        self.size = str(len(tail_data))
        SS = sum([s2 - (s1 * 0.5 * h2) for s1, s2 in tail_data])
        self.exp = round(np.mean([(s1*0.5*h2) for s1, s2 in tail_data]),4) 
        self.obs = np.mean([s2 for s1, s2 in tail_data])
        self.U = SS / h_var
        self.I = (-1*len(tail_data)) / h_var
        self.Z = SS / (len(tail_data)*h_var)**0.5
        if self.side != 'lower': 
            self.pv = stats.norm.cdf(self.Z)
        else:                    
            self.pv = 1 - stats.norm.cdf(self.Z)
            if self.pv == 0.0: self.pv = stats.norm.cdf(-1*self.Z) 
        return self

    


    def get_pi(self, concord_probs): 
        pi, pt, n = np.mean(concord_probs), sum(concord_probs), len(concord_probs)

    def run_mend(self, tail_data, h2, h_var, h_std):
        self.size = str(len(tail_data))
        n = float(len(tail_data))
        idxMin, idxMax = tail_data[0][0], tail_data[-1][0]
        if self.side == 'lower':
            concord_probs = [stats.norm.cdf(idxMax, (s1*0.5*h2), h_std) for s1, s2 in tail_data]
            r = sum([s2 <= idxMax for s1, s2 in tail_data])
        else:
            concord_probs = [1-stats.norm.cdf(idxMin, (s1*0.5*h2), h_std) for s1, s2 in tail_data]
            r = sum([s2 >= idxMin for s1, s2 in tail_data])
                     
        pi = np.mean(concord_probs)
        self.U = (r - n*pi) / (pi*(1-pi))
        self.I = -n / (pi*(1-pi))
        self.Z = (r - n*pi) / (n*pi*(1-pi))**0.5
        self.exp = round(n*pi,4) 
        self.exp = sum(concord_probs) 
        self.obs = r
        
        self.pv = stats.norm.cdf(-1*self.Z)
        if self.pv < 0.01 and self.obs - self.exp < 2: self.pv = 0.1 + random.random()/2.0
        return self



class TailMax: 
    def __init__(self, side, tail_data, h2, h_var, h_std, mend_std = 0.5): 
        self.side = side 
        self.tail_data, self.h2, self.h_var, self.h_std = tail_data, h2, h_var, h_std 
        self.cnts, self.probs = [0.0,0.0,0.0], [] 
        for s1,s2 in tail_data: 
            probs = [stats.norm.pdf(s2, 0, 1), stats.norm.pdf(s2, s1, mend_std), stats.norm.pdf(s2, (s1*0.5*h2), h_std)]
            self.probs.append(probs) 
            if probs[0] > probs[1] and probs[0] > probs[2]: self.cnts[0] += 1.0 
            elif probs[1] > probs[0]:                       self.cnts[1] += 1.0 
            else:                                           self.cnts[2] += 1.0 
        self.rates = [round(c/sum(self.cnts),2) for c in self.cnts] 
        
    def brute_force(self):  
        AD, P = [], [i/100.0 for i in range(100)] 
        for i,nv in enumerate(P): 
            ld = [] 
            for j,mv in enumerate([p for p in P if p < 1-nv]): 
                rates = [nv, mv, round(1 - (nv+mv),2)]  
                like = self.get_log_like(rates) 
                ld.append([like, rates]) 
            AD.append(sorted(ld)[0]) 
        scr, rates = sorted(AD)[0] 
        return rates  
            

    def get_log_like(self, rates): 
        log_like = 0 
        for p1,p2,p3 in self.probs: 
            log_like +=  -math.log((p1*rates[0]) + (p2*rates[1]) + (p3*rates[2]))
        return log_like  





#############################################################################################
##################################  H2 ESIMATION   ##########################################
#############################################################################################


class ConditionalHeritability:
    def __init__(self, pairs):
        self.pairs = pairs
        self.keys = sorted(self.pairs.keys())
        self.h_range = [round(0.0 + (i*0.01), 2) for i in range(101)]

    def estimate(self):
        self.PW = dd(lambda: dd(list))
        self.CI = dd(lambda: dd(list))
        self.set_ptwise()
        self.run_ptwise('All')
        self.run_ptwise('Bod', A=4, B=95)
        self.run_ptwise('Mid', A=35, B=65)
        self.run_ptwise('LoH', A=5, B=40)
        self.run_ptwise('HiH', A=60, B=95)
        self.run_ptwise('LoT', A=-1, B=3)
        self.run_ptwise('HiT', A=96, B=101)
        self.body = self.PW['Bod'][2] 
        
        for k in self.PW.keys(): 
            x,j = self.PW[k][2], self.PW[k][-1] 
            self.CI[k] = [round(x,2), round(max(x-j,0),2), round(min(x+j,1),2)]
        return self 

        


    def set_ptwise(self):
        my_keys = sorted(self.pairs.keys())
        self.h_key = {h: self.set_log_like(h) for h in self.h_range}
        self.estimate1 = sorted([[sum([sum(self.h_key[h][k]) for k in my_keys]), h] for h in self.h_key])[0][1]
        for hn in [round(self.estimate1 + (i+1.0) / 200, 3) for i in range(-60, 60, 2)]:
            if hn > 0.01 and hn < 0.99:  self.h_key[hn] = self.set_log_like(hn)
        return

    def set_log_like(self, h):
        h_likes = dd(list)
        h_var = (1 - (h*h) / 4.0)
        h_frac = 1 / ((h_var ** 0.5) * ((2*math.pi)**0.5))
        for k in self.pairs.keys():
            for b, a in self.pairs[k]:
                Ne = ((b - ((h*a)/2.0)) ** 2) / (2*h_var)
                LP = math.exp(-Ne) * h_frac
                h_likes[k].append(-math.log(LP))
        return h_likes

    def run_ptwise(self, name='All', A=-1, B=101, iterations=50):
        my_keys = sorted([k for k in self.pairs.keys() if k >= A and k <= B])
        my_size = str(int(sum([len(self.pairs[k]) for k in my_keys])))
        if my_size == '0': 
            self.PW[name] = [my_size, 0.0, 0.0, 0.5]
            return  
        
        
        my_estimate = sorted([[sum([sum(self.h_key[h][k]) for k in my_keys]), h] for h in self.h_key])[0][1]
        my_obs, my_lists, my_scores = [], dd(list), dd(list)
        for h in self.h_key.keys():
            for k in my_keys:
                my_scores[h].extend(self.h_key[h][k])
                my_lists[h].append(self.h_key[h][k])
        my_lens = [len(sL) for sL in my_lists[h]]

        if min(my_lens) > 40:
            for itr in range(iterations):
                my_idxs = [sorted(random.sample(range(len(sL)), int(len(sL)/3.0) + 1)) for sL in my_lists[h]]
                my_obs.append(sorted([[sum([sum([L[j] for j in idxs]) for L, idxs in zip(SL, my_idxs)]), k] for k, SL in my_lists.items()])[0][1])
            my_mean, my_std, my_var = round(np.mean(my_obs), 3),  np.std(my_obs), np.var(my_obs)
            self.PW[name] = [my_size, my_estimate, my_mean, round(my_std, 3)]

        elif sum(my_lens) > 1000:
            idxs = range(len(my_scores[h]))
            self.subset_size = int(len(idxs)/3.0)
            for itr in range(iterations):
                random_indexes = sorted(random.sample(idxs, self.subset_size))
                my_obs.append(sorted([[sum([L[ri] for ri in random_indexes]), k] for k, L in my_scores.items()])[0][1])
            my_mean, my_std, my_var = round(np.mean(my_obs), 3),  np.std(my_obs), np.var(my_obs)
            self.PW[name] = [my_size, my_estimate, my_mean, round(my_std, 3)]
        elif my_estimate < 0.01:
            self.PW[name] = [my_size, 0.01, 0.01, 0.01]
        else:
            self.PW[name] = [my_size, my_estimate, my_estimate, 0.45]
        return


#############################################################################################
##################################       MAIN      ##########################################
#############################################################################################




        
class SibAnalysis:
    def __init__(self, args, progress): 
        self.args, self.progress = args, progress 
        
        if self.args.normalize: self.w = open(args.out+'.normalized.result.out','w') 
        else:                    self.w = open(args.out+'.result.out','w') 


    def go(self,f_name, name, plot, sd, key = 'vals'): 
        
        if key == 'vals': self.progress.start('Analyzing Input Values') 
        else:             self.progress.start('Analyzing Normalized Values') 
        #self.name, self.pairs, self.pts = name, sd.pairs[key], sd.data[key]  
        self.name, self.pairs = name, sd.pairs[key]
        self.h2 = self.estimate_heritability() 
        self.tt = self.infer_arch(self.h2.body) 
        self.save_output(name) 
        if self.args.skipPlots: return 
        self.progress.update('Drawing Plot') 
        plot.draw_summary(f_name, name, self.pairs, self.h2, self.tt.results) 
        
        #X = sorted(self.pts.keys())
        #print(X) 
        #print([self.pts[x][1] for x in X]) 
        
        #if self.args.savePlotdata:  self.save_plot_data() 

    def estimate_heritability(self): 
        self.progress.update('Estimating heritability...') 
        h2 = ConditionalHeritability(self.pairs).estimate()
        self.progress.end('h2_full,h2_body = '+str(h2.CI['All'][0])+','+str(h2.CI['Bod'][0]))         
        return h2 

    def infer_arch(self,bodyH2): 
        self.progress.update('Inferring Tail Architecture...') 
        tt = TailTests(self.args, self.pairs, bodyH2).calculate()
        return tt 


    def save_output(self,name):
        self.progress.update('Saving Output') 
        h_names, h_locs =    ['All','Bod','Mid','LoH','HiH','LoT','HiT'], ['0-99','5-95','35-65','5-40','60-95','0-4','96-100']
        self.w.write('--------------------------------------- Conditional Heritability Estimates: '+name+'  ----------------------------------------------------------------------------------------\n') 
        self.w.write('%-25s %12s %7s %9s %10s %10s %9s %9s\n' % ('traitName','sibH2','range', 'size', 'h2Init', 'h2Iter', 'h2CiLo','h2CiHi')) 
        for k,r in zip(h_names, h_locs): 
            size, h2_init, h2_iter, h2_err = self.h2.PW[k] 
            self.w.write('%-25s %12s %7s %9s %10s %10s %9s %9s\n' % (name,'sibH2',r,size,h2_init,h2_iter,self.h2.CI[k][1],self.h2.CI[k][2])) 
        self.w.write('--------------------------------------- Inferred Tail Architecture: '+name+' ------------------------------------------------------------------------------------------------\n') 
        self.w.write('%-25s %12s %7s %9s %10s %10s ' % ('traitName','tailTests','tail', 'size', 'idxAvg', 'distPv')) 
        self.w.write('%9s %9s %9s %9s %9s %9s ' % ('novoPv','novoObs','novoExp','mendPv','mendObs','mendExp')) 
        self.w.write('%9s %9s %9s\n' % ('novoRate','mendRate','polyRate'))         
        tails = ['0', '0-1', '0-2', '0-3', '0-4','95-99','96-99','97-99','98-99','99'] 
        for t in tails: 
            r = self.tt.results[t] 
            self.w.write('%-25s %12s %7s %9s %10.3f %10s ' % (name,'tailTests',t, r.size, r.index_avg, r.d_str)) 
            self.w.write('%9s %9.3f %9.3f %9s %9d %9.3f ' % (r.n_str,r.n_obs, r.n_exp, r.m_str,r.m_obs, r.m_exp))
            self.w.write('%9.3f %9.3f %9.3f\n' % (r.rates['novo'], r.rates['mend'], r.rates['poly'])) 
        self.w.write('\n') 
        return self 










def run_script(file_handles, args, command_line):

    #NAMES = ['48','30030','30070','78','30720','20015','30610'] 
    #NK = {'30030': 'Reticulocyte_Count','30610': 'Alkaline_Phosphatase', '78': 'Heel_Bone_Mineral_Density_(BMD)', '30720': 'Cystatin_C','48': 'Waist_Circumference', '20015': 'Sitting_Height', '30070': 'RBC_Distribution_Width'}
    progress = SibProgress(args,command_line) 

    if len(args.names) < len(args.sibfiles): 
        for i,sf in enumerate(args.sibfiles): 
            if i < len(args.names): continue 
            #elif sf.name.split('/')[-1].split('-')[0] in NAMES: 
            #    args.names.append(NK[sf.name.split('/')[-1].split('-')[0]]) 
            #elif 5 < 8: args.names.append(NK[NAMES[i]]) 
            elif len(sf.name.split('/')[-1].split('.')) > 1: args.names.append(".".join(sf.name.split('/')[-1].split('.')[0:-1]).split('-')[0]) 
            else:                                            args.names.append(sf.name.split('/')[-1].split('-')[0]) 
    

    sibAnalysis = SibAnalysis(args, progress)  
    sibPlot     = SibPlot(args, progress) 
    
    for i,(f_handle,f_name) in enumerate(zip(args.sibfiles, args.names)): 
        progress.update('Reading Sibling Data for Trait '+f_name.lower().capitalize()) 
        fd = SibData(f_handle, args)
        progress.end(str(len(fd.fams))+' siblings pairs found')     
        if not args.normalize: vA = sibAnalysis.go(f_handle.name.split('/')[-1], f_name, sibPlot, fd, key='vals') 
        else:                  nA = sibAnalysis.go(f_handle.name.split('/')[-1], f_name, sibPlot, fd, key = 'norm') 
        
        if not args.skipPlots and i > 0 and i % 5 == 0: sibPlot.reset(len(args.names) - (i+1)) 
    if not args.skipPlots and (i%5 != 0 or i ==0): sibPlot.finish() 
    progress.finish() 
    



if __name__ == '__main__':
    import argparse, sys
    usage = "usage: ./%prog [options] data_file"
    parser = argparse.ArgumentParser()
    parser.allow_abbrev=True
    parser.add_argument('sibfiles',nargs='+',type=argparse.FileType('r')) 
    parser.add_argument("--names", nargs='+',default=[],type=str,help="Trait Name(s)")
    parser.add_argument("--out", type=str,default='tailArc',help="Output Prefix")
    parser.add_argument("--nn", type=str,help="Trait Nick Name",metavar='')
    parser.add_argument("--alpha", type=float,default=0.01,help="Cutoff to Identify Complex Architecture",metavar='') 
    parser.add_argument("--skipPlots", action='store_true', default=False,help="Skip Plotting")
    parser.add_argument("--savePlotdata", action='store_true', default=False,help="Save Plotting Data") 
    parser.add_argument("--normalize", action='store_true', default=False,help="Rank Inverse Normal Transform Data (For NonNormal Input)") 
    parser.add_argument("--silent", action='store_true', default=False,help="Suppress Output Stream") 
    parser.add_argument("--randomize", action='store_true', default=False,help="Randomize Siblings") 
    args = parser.parse_args() 
    run_script(args.sibfiles, args, ' '.join(sys.argv)) 






