#2016-04-21 James Ro
#Reads .DAT data files (output of pixie4) and processes data creating various hit pattern combination energy and time spectroscopy data into .csv files

#files automatically produced: Single hit events (anti ch 0,1,2,3) , Double coincidence hit events (6 x 2 types of single energy alignments), Time spectrum data (for double coin events), Triple Coin hit events (4 x 3 single energy alignment combinations), quad coincidence data (1 x 4 combinations)

#additional data produced (parameters set menually)
#: Time ribbon array (3 column: ch a energy, counts, time interval ) - used for plotting 3D ribbon plots (to observe coincident energy spectrum throughout different time intervals between the coincidence events)
#: 2D coincidence plot data (for 2D histogramming)

#to add: individual channel combined total events

#NOTE: ALL DSP ENERGY VALUES ARE TWICE THE value of those displayed in PIXIE4 MCA spectrum (need to divide processed data energies by 2 even before applying calibration factors)
#Current MCA spectrum KeV/Bin ratio: HPGE ~0.294 (ch 0,1), BC408 ~7.558 (ch 2,3)
#divide factor ratios by 2 (HPGE ~0.147 , BC408 = 3.779  )

#1 unit of time stamp/difference value = 13.333 ns

#User Parameters
factor1 = 0.147 #hpge calib/2 factor
factor2 = 3.779 #bc408 plastic scint
CoinWindow = 200 #200*13.33ns = 2400 ns window
filename1 ='background(all-48hr-CW20k-CD7573).dat'
#looks for above file in current directory the script is located in
#######################################################
import string
import glob
import os
import numpy as np
import time

folder = os.getcwd()
os.chdir(folder)

###########################
t1 = time.clock()
filename = filename1
file = open(folder + "\\" + filename, 'r')


print("begin processing, please wait...")
##################################################################################
eventNo = ''
channelNo = ''
energy = ''
trigTime  = ''
eventInt, channelInt, energyInt , trigInt , tdiffInt = 0, 0,0,0,0
start = 'Event No'
###############################################################################
eventTot = []
channelTot =[]
energyTot=[]
trigTot=[]
####################
#Index of all 15 hit patterns in order Follows permutation pattern from below
# 4x4 = 16, however 16th pattern is 0 0 0 0 (not used)
totalIndex= [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]]
####################

#begin .dat read
for line in file:
    if start in line: #continue reading next line until value header
        break
for line in file:#continue reading from last line, begin collecting data
    if line[0].isdigit(): #skip blank lines
        eventNo = line[0:9]
        eventInt = int(eventNo)
        channelNo = line[9:18]
        channelInt = int(channelNo)

        energy = line[21:29]
        energyInt = int(energy)

        trigTime = line[30:38]
        trigInt = int(trigTime)

        eventTot.append(eventInt)
        channelTot.append(channelInt)
        energyTot.append(energyInt)
        trigTot.append(trigInt)

#####################
dataArray = np.column_stack((eventTot, channelTot, energyTot, trigTot ))
maxEvent = max(eventTot)
combArray = np.empty((maxEvent+1,4,2),dtype = int) #xyz 3D array:  #row = events, column = [channel 0,1,2,3], z (height)= [energy, trigger]
lastEvent = 0
CoinIndex = [] #Index of event's hitpatterns
CoinPattern = 0 #Temporary placeholder for hit pattern
################
for row in dataArray:

    i = row[0] #event
    j = row[1] #channel

    currentEvent = i

    combArray[i,j,0] = row[2] #energy
    combArray[i,j,1] = row[3] #trig time

    #hit pattern indexing
    if currentEvent != lastEvent or currentEvent==maxEvent: #fix code to include proper coin pattern for last event
        CoinIndex.append(CoinPattern)

        #Organizing Data into various Hit Patterns (saving index of locations of interest)
        if CoinPattern == 1: #CH 0 Anti (0001)
            totalIndex[0].append(lastEvent)
        elif CoinPattern == 10: #CH1 Anti (0010)
            totalIndex[1].append(lastEvent)
        elif CoinPattern == 100: #CH2 ANTI (0100)
            totalIndex[2].append(lastEvent)
        elif CoinPattern == 1000: #CH3 ANTI (1000)
            totalIndex[3].append(lastEvent)
        elif CoinPattern == 11: #CH 0-1  D-COIN (0011)
            totalIndex[4].append(lastEvent)
        elif CoinPattern == 101: #CH 0-2 D-COIN (0101)
            totalIndex[5].append(lastEvent)
        elif CoinPattern == 1001: #CH 0-3 D-COIN (1001)
            totalIndex[6].append(lastEvent)
        elif CoinPattern == 110: #CH 1-2 D-COIN (0110)
            totalIndex[7].append(lastEvent)
        elif CoinPattern == 1010: #CH 1-3 D-COIN (1010)
            totalIndex[8].append(lastEvent)
        elif CoinPattern == 1100: #CH 2-3 D-COIN (1100)
            totalIndex[9].append(lastEvent)
        elif CoinPattern == 111: #CH 0-1-2 T-COIN (0111)
            totalIndex[10].append(lastEvent)
        elif CoinPattern == 1101: #CH 0-2-3 T-COIN (1101)
            totalIndex[11].append(lastEvent)
        elif CoinPattern == 1011: #CH 0-1-3 T-COIN (1011)
            totalIndex[12].append(lastEvent)
        elif CoinPattern == 1110: #CH 1-2-3 T-COIN (1110)
            totalIndex[13].append(lastEvent)
        elif CoinPattern == 1111: #CH 1-2-3-4 Q-COIN (1111)
            totalIndex[14].append(lastEvent)
        CoinPattern = 0 #reset coinpattern at end of each event sequence
#below is a counter for tracking channels recorded for that particular event
    if j == 0:
        CoinPattern += 1
    if j == 1:
        CoinPattern += 10
    if j == 2:
        CoinPattern += 100
    if j == 3:
        CoinPattern += 1000

    lastEvent = currentEvent

# 170 mb .dat file takes ~30 seconds to process up to here so far
########################################################################
#data processing
#############################
ch0list = []
ch1list = []
ch2list = []
ch1list = []
ch12DCoinList = [[],[],[]] #energy1, energy2, timediff
newdatalist= [[],[],[],[],[],[],[],[],[],[],[],[],[],[],[]] #0 to 14 for the 14 possible hit pattern types, including sublists for coincidence with single energy alignments


#################################################
#Reindexing data for specific hit pattern groups

for i in range(0,15): #iterate through 15 possible hit patterns

    for event in totalIndex[i]:
        ch0e = combArray[event, 0, 0]
        ch0t = combArray[event, 0, 1]
        ch1e = combArray[event, 1, 0]
        ch1t = combArray[event, 1, 1]
        ch2e = combArray[event, 2, 0]
        ch2t = combArray[event, 2, 1]
        ch3e = combArray[event, 3, 0]
        ch3t = combArray[event, 3, 1]

        #single hit patterns (anticoin) list w/ 1 column of energy
        if i== 0:
            newdatalist[0].append(ch0e)
        if i == 1:
            newdatalist[1].append(ch1e)
        if i == 2:
            newdatalist[2].append(ch2e)
        if i == 3:
            newdatalist[3].append(ch3e)

        #Double coincidence hit patterns list with 3 columns (a energy, b energy, time diff)
        if i ==4:
            td = ch0t-ch1t
            newdatalist[4].append([ch0e,ch1e,td])#ch 0-1 coin
        if i ==5:
            td = ch0t-ch2t
            newdatalist[5].append([ch0e,ch2e,td])#ch 0-2 coin
        if i ==6:
            td = ch0t-ch3t
            newdatalist[6].append([ch0e,ch3e,td])#ch 0-3 coin
        if i ==7:
            td = ch1t-ch2t
            newdatalist[7].append([ch1e,ch2e,td]) #Ch12DCoinList
        if i ==8:
            td = ch1t-ch3t
            newdatalist[8].append([ch1e,ch3e,td])#ch 1-3 coin
        if i ==9:
            td = ch2t-ch3t
            newdatalist[9].append([ch2e,ch3e,td])#ch 2-3 coin


        #Triple coincidence hit patterns, list with 3 columns (a energy, b energy, c energy)
        if i==10:
            newdatalist[10].append([ch0e,ch1e,ch2e])#ch 0-1-2 coin
        if i ==11:
            newdatalist[11].append([ch0e,ch2e,ch3e])
        if i ==12:
            newdatalist[12].append([ch0e,ch1e,ch3e])
        if i ==13:
            newdatalist[13].append([ch1e,ch2e,ch3e])

        #Quad coincidence hit pattern, list with 4 column (a, b, c, d energy)
        elif i == 14:
            newdatalist[14].append([ch0e,ch1e,ch2e,ch3e])

#newdatalist will be 15 column, col 0-3 will be single col list, col 4-9 are triple col lists,  col 10-13 will be triple col, col 14 will be quad col


#################################################
name=filename.replace(".dat", "") #Name of original data file to be used for all processed data name scheme
#####################################################################################

singlenames = ["_ch0-anti","_ch1-anti","_ch2-anti","_ch3-anti",]
for i in range(0,4): #histogramming data for single hit Patterns
    try:
        maxe = max(newdatalist[i])

        chdata = np.zeros((maxe+1,2))
        for j in range(0,maxe+1):
            chdata[j,0] = j
        for chenergy in newdatalist[i]:
            chdata[chenergy,1]+=1

        newname = name+singlenames[i]+".csv"
        newdata = np.asarray(chdata)
        np.savetxt(newname, newdata, delimiter=",") #save all single hit anti coincidences into csv file

    except: #exception error handling for empty list (if there were no data avilable for the particular anticoin pattern)
        print(singlenames[i] + " single hit pattern set empty, not produced")
        pass

#histogram adding for double coincidence hit pattern datas
#histogram adding for time spectroscropy, energy spectroscopy (aligned to a), energy spectro (aligned to b)
timenames = ["_ch_0-1-time", "_ch_0-2-time", "_ch_0-3-time", "_ch_1-2-time", "_ch_1-3-time", "_ch_2-3-time"]
doublenames = ["_ch_0-1-align_0", "_ch_0-1-align_1","_ch_0-2-align_0","_ch_0-2-align_2","_ch_0-3-align_0","_ch_0-3-align_3", "_ch_1-2-align_1","_ch_1-2-align_2", "_ch_1-3-align_1","_ch_1-3-align_3", "_ch_2-3-align_2","_ch_2-3-align_3"]

for i in range(4,10): #produce histogram data for double coin hit patterns
    temptdlist=[]#temp lists time diff
    tempealist=[]#temp list energy of a
    tempeblist=[]#temp list energy of b
    for row in newdatalist[i]:
        temptdlist.append(row[2])
        tempealist.append(row[0])
        tempeblist.append(row[1])
    try:
        mintd = min(temptdlist)
        maxtd = max(temptdlist)
        lentd = maxtd-mintd+1
        tdata = np.zeros((lentd,2),dtype=int) #2column list for time spectroscopy datas
        for j in range(0,lentd):
            tdata[j,0] = j+mintd
        for timediff in temptdlist: #histogramming time difference datas
            tdata[timediff-mintd,1] +=1

        maxea = max(tempealist) #max energy of channel a
        maxeb = max(tempeblist) #max energy of channel b

        chaedata = np.zeros((maxea+1,2)) #histrogrammed energy data aligned to channel a
        for j in range(0,maxea + 1):
            chaedata[j,0] = j
        for row in tempealist:
            chaedata[row,1] +=1

        chbedata = np.zeros((maxeb+1,2)) #histrogrammed energy data aligned to channel b
        for j in range(0,maxeb + 1):
            chbedata[j,0] = j
        for row in tempeblist:
            chbedata[row,1] +=1

        tname = name+timenames[i-4]+".csv"
        tdat = np.asarray(tdata,dtype=int)

        ename1 = name+doublenames[(i-4)*2]+".csv"
        edat1 = np.asarray(chaedata)

        ename2 = name+doublenames[((i-4)*2)+1]+".csv"
        edat2 = np.asarray(chbedata)

        np.savetxt(tname ,tdat, delimiter=",") #save double coin time spectroscopy data as csv
        np.savetxt(ename1, edat1, delimiter=",")# save double coin energy spec, align to a
        np.savetxt(ename2, edat2, delimiter=",")# save double coin energy spec, align to b
    except:
        errorname=doublenames[(i-4)*2]
        errorname=errorname[:-8]
        print(errorname + " double coincidence hit pattern set empty, not produced")
        pass

#triple names
triplenames=["_ch_0-1-2-align_0","_ch_0-1-2-align_1","_ch_0-1-2-align_2","_ch_0-2-3-align_0","_ch_0-2-3-align_2","_ch_0-2-3-align_3","_ch_0-1-3-align_0","_ch_0-1-3-align_1","_ch_0-1-3-align_3","_ch_1-2-3-align_1","_ch_1-2-3-align_2","_ch_1-2-3-align_3"]
#produce histogram data for triple coin hit patterns (4x3 single ch energy combinations)
for i in range(10,14):
    tempealist=[]#temp list energy of a
    tempeblist=[]#temp list energy of b
    tempeclist=[]#temp list energy of c
    for row in newdatalist[i]:
        tempealist.append(row[0])
        tempeblist.append(row[1])
        tempeclist.append(row[2])
    try:
        maxea = max(tempealist) #max energy of channel a
        maxeb = max(tempeblist) #max energy of channel b
        maxec = max(tempeclist) #max energy of chanel c

        chaedata = np.zeros((maxea+1,2)) #histrogrammed energy data aligned to channel a
        for j in range(0,maxea + 1):
            chaedata[j,0] = j
        for row in tempealist:
            chaedata[row,1] +=1

        chbedata = np.zeros((maxeb+1,2)) #histrogrammed energy data aligned to channel b
        for j in range(0,maxeb + 1):
            chbedata[j,0] = j
        for row in tempeblist:
            chbedata[row,1] +=1

        chcedata = np.zeros((maxec+1,2)) #histrogrammed energy data aligned to channel c
        for j in range(0,maxec + 1):
            chcedata[j,0] = j
        for row in tempeclist:
            chcedata[row,1] +=1


        ename1 = name+triplenames[(i-10)*3]+".csv"
        edat1 = np.asarray(chaedata)

        ename2 = name+triplenames[((i-10)*3)+1]+".csv"
        edat2 = np.asarray(chbedata)

        ename3 = name+triplenames[((i-10)*3)+2]+".csv"
        edat3 = np.asarray(chcedata)

        np.savetxt(ename1, edat1, delimiter=",")# save double coin energy spec, align to a
        np.savetxt(ename2, edat2, delimiter=",")# save double coin energy spec, align to b
        np.savetxt(ename3, edat3, delimiter=",")# save double coin energy spec, align to c into csv file
    except:
        errorname=triplenames[(i-10)*3]
        errorname=errorname[:-8]
        print(errorname + " triple coincidence hit pattern set empty, not produced")
        pass


#Quadruple Coincidence names
quadnames=["_ch_0-1-2-3-align_0","_ch_0-1-2-3-align_1","_ch_0-1-2-3-align_2","_ch_0-1-2-3-align_3"]
tempealist=[]
tempeblist=[]
tempeclist=[]
tempedlist=[]
#histrogramming for quadruple coincidence hit pattern datas (1x4 single energy alignment combinations)
try:
    for row in newdatalist[14]:
        tempealist.append(row[0])
        tempeblist.append(row[1])
        tempeclist.append(row[2])
        tempedlist.append(row[3])

    maxea = max(tempealist) #max energy of channel a
    maxeb = max(tempeblist) #max energy of channel b
    maxec = max(tempeclist) #max energy of chanel c
    maxed = max(tempedlist) #max energy of chanel d

    chaedata = np.zeros((maxea+1,2)) #histrogrammed energy data aligned to channel a
    for j in range(0,maxea + 1):
        chaedata[j,0] = j
    for row in tempealist:
        chaedata[row,1] +=1

    chbedata = np.zeros((maxeb+1,2)) #histrogrammed energy data aligned to channel b
    for j in range(0,maxeb + 1):
        chbedata[j,0] = j
    for row in tempeblist:
        chbedata[row,1] +=1

    chcedata = np.zeros((maxec+1,2)) #histrogrammed energy data aligned to channel c
    for j in range(0,maxec + 1):
        chcedata[j,0] = j
    for row in tempeclist:
        chcedata[row,1] +=1

    chdedata = np.zeros((maxed+1,2)) #histrogrammed energy data aligned to channel a
    for j in range(0,maxed + 1):
        chdedata[j,0] = j
    for row in tempedlist:
        chdedata[row,1] +=1

    ename1 = name+quadnames[0]+".csv"
    edat1 = np.asarray(chaedata)

    ename2 = name+quadnames[1]+".csv"
    edat2 = np.asarray(chbedata)

    ename3 = name+quadnames[2]+".csv"
    edat3 = np.asarray(chcedata)

    ename4 = name+quadnames[3]+".csv"
    edat4 = np.asarray(chdedata)

    np.savetxt(ename1, edat1, delimiter=",")# save double coin energy spec, align to a
    np.savetxt(ename2, edat2, delimiter=",")# save double coin energy spec, align to b
    np.savetxt(ename3, edat3, delimiter=",")# save double coin energy spec, align to c into csv file
    np.savetxt(ename4, edat4, delimiter=",")# save double coin energy spec, align to c into csv file
except:
    print("ch0-1-2-3 quadrupule coincidence hit pattern set empty, not produced")
    pass





##EXPERIMENTAL SECTIONS BELOW###
##################################################################################################################################

#section below is for creating csv formatted data for 2D coincidence plots
##hit pattern combination datas (for 2D and 3D graph plotting)
ch01dcoin=np.asarray(newdatalist[4])
ch02dcoin=np.asarray(newdatalist[5])
ch03dcoin=np.asarray(newdatalist[6])
ch12dcoin=np.asarray(newdatalist[7])
ch13dcoin=np.asarray(newdatalist[8])
ch23dcoin=np.asarray(newdatalist[9])
#
ch012tcoin=np.asarray(newdatalist[10])
ch023tcoin=np.asarray(newdatalist[11])
ch013tcoin=np.asarray(newdatalist[12])
ch123tcoin=np.asarray(newdatalist[13])
#
ch0123qcoin=np.asarray(newdatalist[14])

try:

    for row in ch13dcoin:
        row[0]*= factor1
        row[1]*= factor2 #apply calibration factor

    for row in ch12dcoin:
        row[0]*= factor1
        row[1]*= factor2 #apply calibration factor (bc408)

    for row in ch123tcoin: #detector 2 with top and bot cosmic plates triple coin (energy alignment to detector 2 and top plate only for double coin data)
        row[0]*= factor1
        row[1]*= factor2 #apply calibration factor (bc408)

    #section for outputting 2D coincidence data (2 column : ch A energy, ch b energy)
    #tempdd = ch01dcoin[:,0:2]
    tempaa = ch13dcoin[:,0:2] #take only the first 2 columns (ch a energy, ch b energy)
    tempbb = ch12dcoin[:,0:2]
    tempcc = ch123tcoin[:,0:2]

    np.savetxt(name+"ch-1-2-doublecoin-data.csv",tempaa,delimiter=",")
    np.savetxt(name+"ch-1-3-doublecoin-data.csv",tempbb,delimiter=",")
    np.savetxt(name+"ch-1-2-3-triplecoin-double-data.csv",tempcc,delimiter=",")
    #np.savetxt(name+"ch-0-1-doublecoin-data.csv",tempdd,delimiter=",")
except:
    print("Indexing error problem with 2D Histogram plotting data, 2D data not produced")
    pass
###############################################################


##3D Ribbon/Stacked time spectrum slices plot
##histogramming 1D double coincidence spectrum energy data with time interval as 3rd axis (for 3D ribbon plot)

try:
    temparray = np.asarray(newdatalist[7])#currently using data set from ch 1-2 coincidence

    maxe = max(temparray[:,0]) + 1

    ribbonlist = np.zeros((maxe,12+1)) #col 0 = energy, col 1= time interval < 0 , col 2 = time interval from 0:10, col3 = time interval 10:20, ...., col12 = >110

    for i in range(0,maxe):
        ribbonlist[i,0] = i
    for row2 in temparray:
        tgroup = 12 #placeholder for time interval range (12 intervals)
        if row2[2] < -20: #under -266 ns
            tgroup = 1
        elif -20 <= row2[2] < -10: #-266 to -133
            tgroup = 2
        elif -10 <= row2[2] < 0: #-133 to 0
            tgroup = 3
        elif 0 <= row2[2] < 10: #0 to 133
            tgroup = 4
        elif 10 <= row2[2] < 20:  # 133 to 266
            tgroup = 5
        elif 20 <= row2[2] < 30: #266 to 399
            tgroup = 6
        elif 30 <= row2[2] < 40: #399 to 533
            tgroup = 7
        elif 40 <= row2[2] < 50: #533 to 666
            tgroup = 8
        elif 50 <= row2[2] < 60 : #666 to 799 ns range
            tgroup = 9
        elif 60 <= row2[2] < 70: #799 to 933
            tgroup = 10
        elif 70 <= row2[2] < 80: #933 to 1066
            tgroup = 11
        elif row2[2] >= 80: #over 1066 ns
            tgroup = 12

        tenergy = row2[0] #extract energy value

        ribbonlist[tenergy, tgroup] += 1 #histogram summing counts

    for row in ribbonlist:
        row[0]*=factor1 #adjust channel energy by calibration factor

    rname = name+"ch-1-2-time-ribbon-calibrated.csv"
    ribbonarray = np.asarray(ribbonlist)
    np.savetxt(rname, ribbonarray, delimiter=",")
except:
    print("Indexing error with ribbon/stacked spectrum plot, csv file not produced")
    pass

###########################################################
# producing 1D energy spectrum data for specific time slices (of coincident events)
# temparray2 = np.asarray(newdatalist[7])#curren for ch 0-2 coincidence
# maxe2 = max(temparray2[:,0]) + 1
#
# slices = np.zeros((maxe2,5))#currently 2 intervals
#
# for i in range(0,maxe2):
#     slices[i,0] = i
# for row in temparray2:
#     if -8 <= row2[2] <= 8: #-106 : +106 ns
#         slices[row[0], 1] += 1 #histogram summing counts
#     if -2 <= row2[2] <= 2: # -26 : 26 ns
#         slices[row[0], 2] += 1 #histogram summing counts
#     if 15 <= row2[2] <= 45: # 200 : 600 ns
#         slices[row[0], 3] += 1 #histogram summing counts
#     if 22 <= row2[2] <= 38: # 293 : 506 ns
#         slices[row[0], 4] += 1 #histogram summing counts
#
# for row in slices:
#     row[0]*=factor1 #adjust channel energy by calibration factor
#
# sname = name+"ch-1-2-time-slices-align_to_ch1.csv"
# slicearray = np.asarray(slices)
# np.savetxt(sname, slicearray, delimiter=",")
#
# #add module for selecting desired hit pattern datas (+ trigger time difference for coincidences)
###################
#takes ~150 seconds to process up to this point for a 170mb .dat file


t2 = time.clock()
processtime = t2 - t1
print('done, processing took {} Seconds' .format(processtime))



##################
#futurework: integrate this script into Gammalyzer program
