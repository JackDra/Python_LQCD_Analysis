#!/usr/bin/env python

import time,datetime
import numpy as np
import socket

str_buffer = ' '*70
str_buffer_small = ' '*20
myeps = np.finfo(0.0).eps

class Timer(object):

    """

    there are probably default classes that can do this, but meh

    """



    ## dont pass in linklist if you just want total time measured
    def __init__(self,linklist=[],name='',rewrite=True,prod_or_sum='prod'):
        socket_name = socket.gethostname()
        if 'JackLappy' in socket_name:
            this_machine = 'jack_lappy'
        elif 'juqueen' in socket_name:
            this_machine = 'juqueen'
        elif 'dev' in socket_name or 'gateway' in socket_name or 'lac-' in socket_name:
            this_machine = 'laconia'
        else:
            this_machine = 'jack_lappy'

        if isinstance(linklist,(int,float)):
            self.llsize = float(int(linklist))
        else:
            if len(linklist) > 0 and isinstance(linklist[0],(list,tuple,np.ndarray)):
                if 'sum' in prod_or_sum:
                    self.llsize = float(sum([len(ival) for ival in linklist]))
                elif 'prod' in prod_or_sum:
                    self.llsize = float(np.prod([len(ival) for ival in linklist]))
            else:
                self.llsize = float(len(linklist))
        self.rewrite = rewrite and this_machine != 'laconia'
        self.name = name
        self.Start()


    def Start(self,noprint=False):
        self.starttime = time.time()
        self.endtime = 0.0
        self.laplist = []
        self.linkcount,self.linkper = 0,0
        if not noprint:
            if self.rewrite:
                print()
                print(self.name+' Starting timer ,' ,str_buffer,'\r', end=' ')
            else:
                print(self.name+' Starting timer ')

    def Stop(self,newl=True):
        self.endtime = time.time()-self.starttime
        self.laplist.append(self.endtime)
        self.linkcount = -1
        self.linkper = 100
        if self.rewrite:
            print(self.name, ' Complete, took: ' ,self.GetTimeStr(self.endtime) , str_buffer)
            if newl: print()
        else:
            print(self.name, ' Complete, took: ' ,self.GetTimeStr(self.endtime))
            if newl: print()
        self.Start(noprint=True)


    def Lap(self,flag=None):
        self.laplist.append(time.time()-self.starttime)
        self.linkcount += 1
        if self.linkcount == self.llsize or self.llsize == 0:
            self.Stop()
        else:
            # if self.llsize == 0.0:
            #     self.linkper = 100
            # else:
            self.linkper = (self.linkcount*100)/self.llsize
            if self.rewrite:
                print(self.GetTimeForm(flag),str_buffer_small,'\r', end=' ')
            else:
                print(self.GetTimeForm(flag))

    def Reset(self):
        self.Start()

    def GetPercent(self):
        return str('{:.1f}'.format(self.linkper)) + '% '

    def GetTimeStr(self,thistime):
        return str(datetime.timedelta(seconds=thistime)) + ' h:m:s '

    def GetTimeLeft(self):
        perdone = self.linkper/100.
        if perdone < 0.001:
            return float(0)
        else:
            return self.laplist[-1]*((1-perdone)/perdone)

    def GetTimeLeftStr(self):
        return self.GetTimeStr(self.GetTimeLeft())

    def GetTimeForm(self,flag=None,screen_size=30):
        if flag is None:
            out_str = self.name+' Current Time: ' + self.GetTimeStr(self.laplist[-1]) + ' at ' + self.GetPercent() + ' Time Left: ' + self.GetTimeLeftStr()
        elif isinstance(flag,(list,tuple,np.ndarray)):
            out_str = self.name+' At '+','.join(map(str,flag))+' Current Time: ' + self.GetTimeStr(self.laplist[-1]) + ' at ' + self.GetPercent() + ' Time Left: ' + self.GetTimeLeftStr()
        else:
            out_str = self.name+' At '+str(flag)+' Current Time: ' + self.GetTimeStr(self.laplist[-1]) + ' at ' + self.GetPercent() + ' Time Left: ' + self.GetTimeLeftStr()
        # print 'DEBUG'
        # print len(out_str),len(out_str) > screen_size
        if len(out_str) > screen_size:
            out_str = out_str.replace(self.name+' ','')
        return out_str
