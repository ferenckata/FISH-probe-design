
# coding: utf-8

# In[ ]:

# Selection from overlapping oligos to have the possible biggest non-overlapping set


# In[1]:

import numpy
import re, os
import sys


# In[12]:

filtereddir=sys.argv[1]
blasteddir=sys.argv[2]


# In[2]:

#filtereddir='/Users/GG/Desktop/Kata/Quim_Magda_probes/human/missingthree/blastiltered/'
#blasteddir=filtereddir


# In[3]:

os.chdir(filtereddir)
files=os.listdir('.')


# In[4]:

def search(oligos):
    for olitm in oligos:
        yield(olitm)


# In[5]:

# create groups
# set minimal distance between consecutive oligos #
for file in files:
    out=file.split('_')[0] +'_' + file.split('_')[1]
    probesfile = '../final_probes/' + out + '_probes.tsv'
    outfasta = '../final_probes/' + out + '_chosen.fa'
    mind = 2
    with open(file,'r') as fi:
        groups = {}
        # set starting intron number #
        oldi = '1'
        g = 1
        # set start coordinate for the group #
        stc = 0
        for line in fi:
            if line.startswith('>'):
                # get the intron number #
                inm = line.split('_')[1].split(':')[0]
                # get GC content data #
                gc = float(line.split(':')[2])
                # get the start coordinate #
                coord = int(line.split(':')[1])
                # pick consecutive oligos and store in their group #
                # store the GC content data for each oligo #
                olig = {}
                olig[coord] = gc
                if stc + 29 + mind > coord and oldi == inm:
                    oldi = inm
                    if g not in groups:
                        groups[g] = [olig,]
                    else:
                        groups[g].append(olig)
                else:
                    oldi = inm
                    g += 1
                    groups[g] = [olig,]
                # adjust starting coordinate for the next group #
                stc = coord
    # choose non-overlapping oligos
    # and the corresponding GC properties
    sets = {}
    for key in groups:
        sets[key] = {}
        sets[key][1] = []
        # start picking the consecuitve oligos from the frst oligo
        # than shift one and starts from the next, etc.
        for shift in range(len(groups[key])):
            c = 0
            # if an element is already a part of a set of oligos skip it
            #print(sets[key])        
            #print(groups[key][shift])
            for d in sets[key].keys():
                #print(list(sets[key][d]))
                if groups[key][shift] in list(sets[key][d]):
                    c = 1
            if c==1:
                continue
                    #print('aww')
            #if (groups[key][shift] in list(sets[key][d]) for d in sets[key].keys()):
                #print('hurray')
             #   continue
            elzero = list(groups[key][shift].keys())[0]
            # otherwise store it as the starting oligo of the new set
            sets[key][shift+1] = [groups[key][shift],]
            for element in groups[key]:
                el = list(element.keys())[0]
                # store the oligo if it is not present yet and non overlapping with the prevous one
                if el > elzero + 29 + mind:
                    elzero = el
                    sets[key][shift+1].append(element)
#########################################################################
### to calculate the mean GC and the std of the GC for each set and pick the best one from each group
### a preferred mean can be chosen from the statistics before
##########################################################################################################

    probestorage = {}
    with open(probesfile,'a') as pr:
        pline = "group id\tprobe id\tthe mean GC\tnumber of oligos\tthe std of the GC\n"
        pr.write(pline)
        optgc = 50
        for ks in sets.keys():
            probe = {}
            for kks in sets[ks].keys():
                lst = []
                oc = []
                for it in sets[ks][kks]:
                    lst.append(list(it.values())[0])
                    oc.append(list(it.keys())[0])
                # average GC
                mn = numpy.mean(lst)
                # std of the GC
                st = numpy.std(lst)
                # size and spread
                # oc stores the relative coordinates, but to get the real ones (important when there are exons and introns)
                # one would need to include it from the first place
                # Tm should have been included, too
                #line = "%s:%s\t%d\t%.2f\t%.2f\n" % (ks,kks,len(sets[ks][kks]),mn,st)
                #pr.write(line) 
                probe[kks] = [len(sets[ks][kks]),mn,st]
            maxlen = 0
            bestgc = 0.5
            for item in probe:
                if int(probe[item][0]) >= maxlen:
                    maxlen = int(probe[item][0])
                    if abs(float(probe[item][1])-0.5) <= abs(bestgc):
                        bestgc = float(probe[item][1])-0.5
                        stdgc = float(probe[item][2])
                        if (maxlen > 1 or abs(bestgc) < 0.1) and stdgc < 0.1:
                            bestprobe = item
                            probestorage[ks] = bestprobe
            line = '%s\t%s\t%.2f\t%d\t%.2f\n' % (ks,bestprobe,bestgc+0.5,maxlen,stdgc)
            pr.write(line)
    oligos = []
    for gr in probestorage.keys():
        a = sets[gr][probestorage[gr]]
        for i in a:
            oi = ':' + str(list(i.keys())[0]) + ':' + '%.6f' % float(list(i.values())[0])
            oligos.append(oi)
            
# create the fasta file with the chosen oligos
    headr = re.compile('>')
    f=open(file,'r')
    co=open(outfasta,'a')
    ptrn = search(oligos)
    pattern = next(ptrn)
    sr = re.compile(pattern)
    for fline in f:
        if headr.match(fline):
            if sr.search(fline) == None:
                sq = 0
                continue
            else:
                sq = 1
                try:
                    pattern = next(ptrn)
                    sr = re.compile(pattern)
                    co.write(fline)
                except:
                    break
        elif sq == 1:
            co.write(fline)
    f.close()
    co.close()

