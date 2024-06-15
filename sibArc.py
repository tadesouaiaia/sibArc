#!/usr/bin/python3
import sys
from collections import defaultdict as dd
import numpy as np
from scipy import stats
import random
import math
import matplotlib as mpl
import matplotlib.pyplot as plt
import pandas as pd
plt.rc('text', usetex=True)
mpl.use("Cairo")

#############################################################################################
##################################  PROGRESS STREAM  ########################################
#############################################################################################


class SibProgress:
    def __init__(self, args, command_line, ACTIVE=True):
        self.args = args
        self.ACTIVE = ACTIVE 
        self.initialize(command_line) 

    def initialize(self, command_line):  
        self.out1  = open(self.args.out+'.progress.log','w') 
        self.out2  = sys.stderr 
        self.space = '' 
        self.show('\nSibArc Begins: '+command_line+'\n')
        self.space = ''
        self.show('   Input File: '+args.sibfile.name+'\n') 
        self.show('Output Prefix: '+args.out+'\n\n')
        self.loc = None

    def show(self, msg, space=''):
        
        if space == 'NA': myspace = '' 
        else:             myspace = self.space

        self.out1.write(myspace+msg) 
        if self.ACTIVE: 
            self.out2.write(myspace+msg)
            self.out2.flush() 

    def start(self,msg):
        if self.loc is not None: self.show('...Finished\n','NA') 
        self.space = '' 
        self.show('\n'+msg+':\n') 
        self.loc = None 
        self.space = '  '

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
    def __init__(self, args, out_prefix):

        self.args, self.figname, self.WD, self.HT = args, out_prefix+'.fig.png', 12, 8
        self.fig = mpl.pyplot.gcf()
        self.fig.set_size_inches(self.WD, self.HT)
        self.setup()
        self.ax = plt.subplot2grid((1, 1), (0, 0), rowspan=1, colspan=1)

    def setup(self):
        self.color_rainbow = plt.get_cmap('coolwarm')
        self.herits = [round(i * 0.05, 2) for i in range(21)]
        self.h_colors = [self.color_rainbow(
            1.*i/float(len(self.herits))) for i in range(len(self.herits))]
        self.color_key = {h: c for h, c in zip(self.herits, self.h_colors)}

    def draw_summary_table(self, pairs, bH, tests):
        self.bH, s1, s2, s_len = bH, [], [], [] 
        X = sorted(pairs.keys())
        for x in X:
            s1.append(np.mean([p[0] for p in pairs[x]]))
            s2.append(np.mean([p[1] for p in pairs[x]]))
            s_len.append(len(pairs[x]))
        for i, h in enumerate(self.herits):
            sh = [(s * h)/2.0 for s in s1]
            if i == 0:
                self.ax.plot(
                    X, sh, color=self.color_key[h], linewidth=3, alpha=0.7)
            else:
                z1 = [(a+b+b)/3.0 for a, b in zip(sh, lp)]
                z2 = [(a+b)/2.0 for a, b in zip(sh, lp)]
                self.ax.plot(X, z1, color=self.color_key[h], linewidth=4, alpha=0.3)
                self.ax.plot(X, z2, color=self.color_key[h], linewidth=4, alpha=0.3)
                self.ax.plot(X, sh, color=self.color_key[h], linewidth=4, alpha=0.3)
            lp = [z for z in sh]

        if round(max(abs(lp[0]), lp[-1])+0.05, 1) < 1.5:   yMax = 1.5
        elif round(max(abs(lp[0]), lp[-1])+0.05, 1) < 1.8: yMax = 1.8
        else:                                              yMax = 2
        if self.args.name is not None: self.ax.text(0, yMax*0.96, self.name, ha='left',va='top', fontsize=30, fontweight='bold')

        sE = [(s*self.bH)/2.0 for s in s1]
        self.ax.plot(X, sE, color='k', linewidth=2, alpha=1.0, zorder=2)
        
        for i, (x, y) in enumerate(zip(X, s2)):
            if x > 1 and x < 99:
                self.ax.scatter(x, y, s=25, alpha=0.75,marker='o', color='grey', edgecolor='k', zorder=50)
        self.ax.scatter(0, tests['0']['NOVO'].exp, edgecolor='k', marker='v', zorder=100,facecolor='whitesmoke', linewidth=2.5, s=200, clip_on=False)
        self.ax.scatter(99, tests['99']['NOVO'].exp, edgecolor='k', marker='v', zorder=100,facecolor='whitesmoke', linewidth=2.5, s=200, clip_on=False)
        
        for loc in ['0','99']: 
            nT,mT,dT = tests[loc]['NOVO'], tests[loc]['MEND'], tests[loc]['DIST'] 
            pvs = sorted([[nT.pv,nT,'NOVO'],[mT.pv,mT,'MEND'],[dT.pv,dT,'DIST']])
            if pvs[0][0] > 0.05:  self.ax.scatter(int(loc), nT.obs, edgecolor='k', marker='h', zorder=100,facecolor='grey', linewidth=2.5, s=200, clip_on=False)
            else: 
                if pvs[0][-1] == 'NOVO': 
                    if pvs[1][0] > 0.05: myClr = 'lime' 
                    else:                myClr = 'green' 
                elif pvs[1][0] > 0.05 and pvs[0][-1] == 'MEND': myClr = 'red' 
                else:                                           myClr = 'purple' 
                self.ax.scatter(int(loc), nT.obs, edgecolor='k', marker='h', zorder=100,facecolor=myClr, linewidth=2.5, s=200, clip_on=False)
        self.reset_fig(yMax) 
        return



    

    def reset_fig(self, yMax):
        self.ax.set_yticks([-2, -1, 0, 1, 2])
        self.ax.set_xlim(-1.5, 100.5)
        self.ax.set_ylim(-yMax, yMax)
        self.ax.set_xlabel('Index Sibling Rank \% ($s_{(1)}$)', fontsize=22)
        self.ax.set_ylabel('Conditional Sibling Z-Value ($s_2$)', fontsize=22)
        plt.subplots_adjust(left=0.08, bottom=0.09, right=0.98,
                            top=0.99, wspace=0.14, hspace=0.46)
        plt.savefig(self.figname, dpi=300)
        plt.clf()
        self.fig = mpl.pyplot.gcf()
        self.fig.set_size_inches(self.WD, self.HT)
        self.ax = plt.subplot2grid((1, 1), (0, 0), rowspan=1, colspan=1)


#############################################################################################
##################################  FAMILY DATA  ############################################
#############################################################################################


class SibData:
    def __init__(self, file_handle):
        self.fams, self.members = [], []
        self.read_file(file_handle)
        self.collate()
        self.summarize() 


    def read_file(self, file_handle):
        for i, line in enumerate(file_handle):
            line = line.split(",")
            try:
                k1, k2 = Kid(float(line[0]), str(float(i))), Kid(float(line[1]), str(i+0.5))
            except ValueError:
                if i == 0:
                    continue
                else:
                    raise Exception('Incorrect File Format')

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
            random.shuffle(kids)
            rank = kids[0].rank
            self.pairs['vals'][kids[0].rank].append([kids[0].val, kids[1].val])
            self.pairs['norm'][kids[0].rank].append([kids[0].z, kids[1].z])
        return 

    def summarize(self): 
        self.data = dd(lambda: dd(list)) 
        for k, p in self.pairs.items():
            S1,S2 = [], []
            for pt in p.keys(): 
                S1 = [DZ[0] for DZ in p[pt]]
                S2 = [DZ[1] for DZ in p[pt]]
                m1, m2, sm  = np.mean(S1), np.mean(S2), stats.sem(S2) 
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


#class TailTests:
#    def __init__(self, name, pairs, h2):
#        self.name, self.pairs, self.h2 = name, pairs, h2

class TailTests:
    def __init__(self, pairs, h2):
        self.pairs, self.h2 = pairs, h2

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

    def calculate(self):
        lower_tail = [[0],   [0, 1], [0, 1, 2],  [
            0, 1, 3],  [0, 1, 2, 3, 4], [0, 1, 2, 3, 4, 5]]
        upper_tail = [[99], [98, 99], [97, 98, 99], [96, 97, 98, 99], [
            95, 96, 97, 98, 99], [94, 95, 96, 97, 98, 99]]
        tail_locs = [0 for lt in lower_tail] + [1 for ut in upper_tail]

        #bodyH = self.h2['Bod'][1]
        my_var = 1 - (self.h2*self.h2)/4.0
        my_std = my_var ** 0.5
        RES = dd(lambda: {})

        for vals, loc in zip(lower_tail+upper_tail, tail_locs):

            vs = str(vals[0])
            if vals[-1] != vals[0]:
                vs += '-'+str(vals[-1])
            my_data = self.get_data(loc, vals)
            
            index_sibs = sorted([md[0] for md in my_data]) 
            
            if loc == 0: index_sib = index_sibs[-1] 
            else:        index_sib = index_sibs[0] 
            
            
            RES[vs]['NOVO'] = TailTest(loc,index_sib).run_novo(my_data, self.h2, my_var, my_std)
            RES[vs]['MEND'] = TailTest(loc,index_sib).run_mend(my_data, self.h2, my_var, my_std)
            RES[vs]['DIST'] = TailTest(loc,index_sib).run_dist(my_data, self.h2, my_var, my_std)

        return RES





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
        return self

    def run_mend(self, tail_data, h2, h_var, h_std):
        self.size = str(len(tail_data))
        n = float(len(tail_data))
        idxMin, idxMax = tail_data[0][0], tail_data[-1][0]
        if self.side == 'lower':
            concord_probs = [stats.norm.cdf(
                idxMax, (s1*0.5*h2), h_std) for s1, s2 in tail_data]
            r = sum([s2 <= idxMax for s1, s2 in tail_data])
        else:
            concord_probs = [
                1-stats.norm.cdf(idxMin, (s1*0.5*h2), h_std) for s1, s2 in tail_data]
            r = sum([s2 >= idxMin for s1, s2 in tail_data])
        pi = np.mean(concord_probs)
        self.U = (r - n*pi) / (pi*(1-pi))
        self.I = -n / (pi*(1-pi))
        self.Z = (r - n*pi) / (n*pi*(1-pi))**0.5
        self.exp = round(n*pi,4) 
        self.obs = r
        self.pv = 1 - stats.norm.cdf(self.Z)
        if self.pv < 0.01 and self.obs - self.exp < 2:
            self.pv = 0.1 + random.random()/2.0
        return self



class TailMax: 
    def __init__(self, side, tail_data, h2, h_var, h_std, mend_std = 0.1): 
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
        self.set_ptwise()
        self.run_ptwise('All')
        self.run_ptwise('Bod', A=4, B=95)
        self.run_ptwise('Mid', A=35, B=65)
        self.run_ptwise('LoH', A=5, B=40)
        self.run_ptwise('HiH', A=60, B=95)
        self.run_ptwise('LoT', A=-1, B=3)
        self.run_ptwise('HiT', A=96, B=101)
        self.body = self.PW['Bod'][1] 
        return self 

        

        #self.bH, s1, s2, s_len = round(h2['Bod'][1], 2), [], [], []
        #   return self.PW

    def set_ptwise(self):
        my_keys = sorted(self.pairs.keys())
        self.h_key = {h: self.set_log_like(h) for h in self.h_range}
        self.estimate1 = sorted(
            [[sum([sum(self.h_key[h][k]) for k in my_keys]), h] for h in self.h_key])[0][1]
        for hn in [round(self.estimate1 + (i+1.0) / 200, 3) for i in range(-60, 60, 2)]:
            if hn > 0.01 and hn < 0.99:
                self.h_key[hn] = self.set_log_like(hn)
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
        
        
        my_estimate = sorted(
            [[sum([sum(self.h_key[h][k]) for k in my_keys]), h] for h in self.h_key])[0][1]
        my_obs, my_lists, my_scores = [], dd(list), dd(list)

        for h in self.h_key.keys():
            for k in my_keys:
                my_scores[h].extend(self.h_key[h][k])
                my_lists[h].append(self.h_key[h][k])
        my_lens = [len(sL) for sL in my_lists[h]]
        

        if min(my_lens) > 40:
            for itr in range(iterations):
                my_idxs = [sorted(random.sample(range(len(sL)), int(
                    len(sL)/3.0) + 1)) for sL in my_lists[h]]
                my_obs.append(sorted([[sum([sum([L[j] for j in idxs]) for L, idxs in zip(
                    SL, my_idxs)]), k] for k, SL in my_lists.items()])[0][1])
            my_mean, my_std, my_var = round(
                np.mean(my_obs), 3),  np.std(my_obs), np.var(my_obs)
            self.PW[name] = [my_size, my_estimate, my_mean, round(my_std, 3)]

        elif sum(my_lens) > 1000:
            idxs = range(len(my_scores[h]))
            self.subset_size = int(len(idxs)/3.0)
            for itr in range(iterations):
                random_indexes = sorted(random.sample(idxs, self.subset_size))
                my_obs.append(sorted(
                    [[sum([L[ri] for ri in random_indexes]), k] for k, L in my_scores.items()])[0][1])
            my_mean, my_std, my_var = round(
                np.mean(my_obs), 3),  np.std(my_obs), np.var(my_obs)
            self.PW[name] = [my_size, my_estimate, my_mean, round(my_std, 3)]

        elif my_estimate < 0.01:
            self.PW[name] = [my_size, 0.01, 0.01, 0.01]
        else:
            self.PW[name] = [my_size, my_estimate, my_estimate, 0.5]
        return


#############################################################################################
##################################       MAIN      ##########################################
#############################################################################################




class SibAnalysis:
    def __init__(self, name, args, progress, pairs, pts): 
        self.out_prefix, self.args, self.progress, self.pairs, self.pts = args.out+'.'+name, args, progress, pairs, pts
        
         
         

    def go(self): 
        self.progress.update('Estimating heritability...') 
        self.h2 = ConditionalHeritability(self.pairs).estimate()
        self.progress.update('Inferring Tail Architecture...') 
        self.tt = TailTests(self.pairs, self.h2.body).calculate()
        self.progress.update('Saving Output') 
        self.save_output() 
        if self.args.savePlot: 
            self.progress.update('Creating Plot') 
            tP = SibPlot(self.args, self.out_prefix).draw_summary_table(self.pairs, self.h2.body, self.tt) 
        if self.args.savePlotData: 
            self.progress.update('Saving Plot Data') 
            self.save_plot_data() 


    def save_plot_data(self): 
        X = [x for x in range(100)] 
        w = open(self.out_prefix+'.plotData.out','w')
        w.write('%-7s %24s %24s %24s %24s\n' % ('---', 'sib1', 'sib2_exp', 'sib2_obs', 'sib2_err')) 
        for x in X:
            try: 
                s1, s2, s_err = self.pts[x] 
                s_exp = (s1*self.h2.body)/2.0            
            except: 
                s1, s_exp, s2, s_err = 'NA', 'NA', 'NA', 'NA'  
            w.write('%-7s %24s %24s %24s %24s\n' % (x, s1, s_exp, s2, s_err)) 
        w.close() 
        return 



    def save_output(self):
        w = open(self.out_prefix+'.h2.out','w')
        h2_locs =    ['All','Bod','Mid','LoH','HiH','LoT','HiT'] 
        h2_ranges =  ['0-99','5-95','35-65','5-40','60-95','0-3','97-100']
        w.write('%-7s %24s %24s %24s %24s\n' % ('range', 'size', 'h2_est', 'h2_iter', 'h2_err')) 
        for k,r in zip(h2_locs, h2_ranges):
            size, h2_init, h2_iter, h2_err = self.h2.PW[k] 
            w.write('%-7s %24s %24s %24s %24s\n' % (r,size,h2_init,h2_iter,h2_err)) 
        w.close() 
        w = open(self.out_prefix+'.result.out','w')
        tails = ['0', '0-1', '0-2', '0-3', '96-99','97-99','98-99','99'] 
        w.write('%-7s %9s %10s %10s ' % ('tail', 'size', 'index_sib', 'ks-pv')) 
        w.write('%14s %12s %12s ' % ('novo_pv','novo_obs','novo_exp')) 
        w.write('%14s %12s %12s ' % ('mend_pv','mend_obs','mend_exp')) 
        w.write('%12s %12s %12s ' % ('novo_rate','mend_rate','poly_rate')) 
        w.write('\n') 
        for t in tails: 
            nv, md, dt = self.tt[t]['NOVO'], self.tt[t]['MEND'], self.tt[t]['DIST'] 
            tail_size, index_sib   = nv.size, nv.index_sib
            n_rate, m_rate, p_rate = dt.n_rate, dt.m_rate, dt.p_rate 
            n_pv, m_pv, d_pv = nv.pv, md.pv, dt.pv 
            n_exp, n_obs = nv.exp, nv.obs
            m_exp, m_obs = md.exp, md.obs
            if dt.pv < 0.5:     w.write('%-7s %9s %10.5f %10.2e ' % (t, tail_size, index_sib, dt.pv)) 
            else:               w.write('%-7s %9s %10.5f %10.5f ' % (t, tail_size, index_sib, dt.pv))
            if nv.pv < 0.05:    w.write('%14.4e %12.5f %12.5f ' % (nv.pv, nv.obs, nv.exp)) 
            else:               w.write('%14.5f %12.5f %12.5f ' % (nv.pv, nv.obs, nv.exp)) 
            if md.pv < 0.05:    w.write('%14.2e %12d %12.5f ' % (md.pv, md.obs, md.exp)) 
            else:               w.write('%14.5f %12d %12.5f ' % (md.pv, md.obs, md.exp)) 
            w.write('%12.3f %12.3f %12.3f' % (dt.n_rate, dt.m_rate, dt.p_rate)) 
            w.write('\n') 
        w.close() 
        return self 


def run_script(file_handle, args, command_line):
    if args.out is None: args.out = 'tailArc.'+file_handle.name.split('/')[-1].split('.')[0].split('-')[0]
    progress = SibProgress(args,command_line) 
    progress.update('Reading Sibling Data') 
    fd = SibData(file_handle)
    if not args.normalize: 
        progress.start('Analyzing Input Values') 
        vA = SibAnalysis('input',args, progress, fd.pairs['vals'], fd.data['vals']).go() 
    else: 
        progress.start('Analyzing Normalized Values') 
        nA = SibAnalysis('normalized',args, progress, fd.pairs['norm'], fd.data['norm']).go() 
    progress.finish() 
    



if __name__ == '__main__':
    import argparse, sys
    usage = "usage: ./%prog [options] data_file"
    parser = argparse.ArgumentParser() 
    parser.add_argument('sibfile',type=argparse.FileType('r')) 
    parser.add_argument("--out", type=str,help="Output file prefix")
    parser.add_argument("--name", type=str,help="Output file prefix")
    parser.add_argument("--savePlot", action='store_true', default=False,help="Whether we should generate the plots")
    parser.add_argument("--savePlotData", action='store_true', default=False,help="Save Plotting Data") 
    parser.add_argument("--normalize", action='store_true', default=False,help="Renormalize Data") 
    args = parser.parse_args() 
    run_script(args.sibfile, args, ' '.join(sys.argv)) 






