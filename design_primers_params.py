#!/usr/bin/env python
##design primers to features in multiple sequences
##usage: python  design_primers.py <fasta-file> <gff file> <file of target IDs> <prod_min_size> <prod_max_size>

##CAUTION will only reliably work with  Primer3 version 1.1.4 or earlier ##


#Copyright 2012 John McCallum & Leshi Chen
#New Zealand Institute for Plant and Food Research
#This program is free software: you can redistribute it and/or modify
#     it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

import os
import StringIO
import re
import tempfile
import subprocess
import copy
import sys
from BCBio import GFF
from BCBio.GFF import GFFExaminer
from Bio import SeqIO
from Bio.Emboss.Applications import Primer3Commandline
from Bio.Emboss import Primer3
import urllib
import urllib2
import xml.etree.ElementTree as ET
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp
from scipy import interpolate
from scipy.interpolate import interp1d
from matplotlib.backends.backend_pdf import PdfPages
import httplib
import socket # For other timeout occurrances urllib2 doesn't catch
import pdb # debugger!

# Query UW melt prediction service for a single sequence, returning default array for helicity, assuming within temperature range of 65-95
def getmelt(input_seq):
  url='http://www.dna.utah.edu/db/services/cgi-bin/udesign.cgi'
	values = {'seq' : input_seq, 'rs':0, 'dmso':0,'cation': 20 ,'mg': 2} # Note the buffer conditions
	data = urllib.urlencode(values)
	# Exception handling, if unable to request data/url etc.
	try:
	    req = urllib2.Request(url, data)
	    response = urllib2.urlopen(req, timeout=5)
	except urllib2.HTTPError, e:
	    print "\n", "_" * 40, "\nRequest not fulfilled by server: \n", e.filename, "\n", e, "\n \n", e.read(), "\n", e.hdrs, "\n", "_" * 40, "\n",
	except urllib2.URLError, e:
	    print "\n", "_" * 40, "\nFailed to reach server: ", e.filename, "\n", e, "\n \n", e.read(), "\n", e.hdrs, "\n", "_" * 40, "\n"
	except socket.timeout:
	    print "\n", "_" * 40, "\nFailed to reach server: timed out \n", "_" * 40, "\n"
	except httplib.HTTPException, e:
	    print "\n", "_" * 40, "\nGeneric HTTP error: \n", e, "\n\n", "_" * 40, "\n"
	except Exception:
	    import traceback
	    print "\n", "_" * 40, "\nGeneric error at 'getmelt(input_seq)' \n", "_" * 40, "\n"
	melt_data = response.read()
	tree = ET.fromstring(melt_data)
	helicity = [amp.find('helicity').text.split() for amp in tree.findall('amplicon')]
	hels = np.array(helicity[0], dtype=np.float32).transpose() # helicity[0] used because the default retreived data is a list of 3 lists (for WT, mu and hets). The 3 lists are identical for our data (no IUPAC) so only [0] is used.
	return hels
	
def getTm(hel_array):
	temps = np.arange(65,95.5,0.5) # Temperature range of 65-95. Step of 0.5
	tck = interpolate.splrep(temps,hel_array,s=0)
	xnew = np.arange(65,95.5,0.05)
	yder = interpolate.splev(xnew,tck,der=1) # der=1, first derivative
	return xnew[yder.argmin()] # Returns the x value corresponding to the minimum (peak) y value -> Tm
	
def graphplot(hel_array):
	hels = hel_array
	temps = np.arange(65,95.5,0.5) # Temperature range of 65-95. Step of 0.5
	tck = interpolate.splrep(temps,hel_array,s=0)
	xnew =np.arange(65,95.5,0.05)
	yder = interpolate.splev(xnew,tck,der=1) # der=1, first derivative
	ypeak = yder.min()
	xpeak = xnew[yder.argmin()] # Tm; x value corresponding to y min
	plt.plot(temps, hels, '--', xnew, yder, xpeak, ypeak, 'o')
	plt.xlabel('Temperature')
	plt.ylabel('Helicity')
	plt.legend(['WT melt curve', 'WT differentiated', 'TM', 'Mutant melt curve', 'Mutant differentiated'])

pdfpage = PdfPages('pdfoutput.pdf')

in_file = sys.argv[1]
gff_file = sys.argv[2]
target_file =  sys.argv[3]
prod_min_size = int(sys.argv[4])
prod_max_size = int(sys.argv[5])
optimum_length = int(sys.argv[6])									## Recommended: 20
max_tm_diff = int(sys.argv[7])                                      ## Recommended: 1
optimum_tm = int(sys.argv[8]) 										## Recommended: 59-61
opt_GC_percent = int(sys.argv[9])                                   ## Recommended: 50
maxpolyx = int(sys.argv[10])                                        ## Recommended: less than 4
gc_clamp = int(sys.argv[11]) 										## Recommended: 2 G/Cs at the very end

## target is specified in start, end format 
productsizerange = str(prod_min_size) + "," + str(prod_max_size)
# open input files
in_seq_handle = open(in_file)
in_gff_handle = open(gff_file)
in_target_handle = open(target_file)
# read  target feature IDs into list
targets = in_target_handle.readlines()
in_target_handle.close()
## and create a hit list of sequences from this
target_seq_id_list = list(set([line.split(":")[0] for line in targets]))
## create iterator returning sequence records
for myrec in SeqIO.parse(in_seq_handle, "fasta"):
    #check if this sequence is included in the target list
    if myrec.id in target_seq_id_list:
        ##create sequence dictionary so we can add in gff annotations
        seq_dict = {myrec.id : myrec}
        ##just limit to gff annotations for this sequence
        limit_info = dict(gff_id = [ myrec.id ])
        ##rewind gff filehandle
        in_gff_handle.seek(0)
        ##read annotations into sequence dictionary for this sequence record only 
        annotations = [r for r in GFF.parse(in_gff_handle, base_dict=seq_dict,limit_info=limit_info)]
        ##if there are any annotations, then  proceed. 
        if annotations:
            rec=annotations[0]
            ##iterate over list of target IDs
            for t in targets:
                target_ID = t.strip('\n')
                target_annotations = [f for f in rec.features if f.id == target_ID ]
                if target_annotations:
                    mytarget = target_annotations[0]
                    #create temporary files
                    tempfastaFile = tempfile.mktemp()
                    tempproutfile = tempfile.mktemp()
                    #just consider slice of sequence in a window of +/- prod_max_size  bp
                    ##from feature UNLESS feature is close to end
                    ##Note that slice is zero-based
                    featLocation = mytarget.location.start.position 
                    if featLocation > prod_max_size:
                        slice_start = featLocation - prod_max_size
                        featPosition = prod_max_size  
                    else:
                        slice_start = 0
                        featPosition = featLocation
                    if (len(rec) - featLocation) < prod_max_size:
                        slice_end = len(rec)
                    else:
                        slice_end = featLocation + prod_max_size
                    ###grab slice of sequence fom this window.
                    targetRec = rec[slice_start:slice_end]
                    matching_feature = [f for f in targetRec.features if f.id == mytarget.id]
                    if matching_feature:
                        target_feat = matching_feature[0]
                        if target_feat.location.start.position == 0:
                            target_feat.location.start.position = 1
                        #we get the mask features by removing the target...all features are masked as just using snp and indels
                        ##a smarter filter could be added 
                        ##note use of list copy to avoid possible side-effects
                        exclude_feat = list(targetRec.features)
                        exclude_feat.remove(target_feat)
                        ##print'targetRec.features',  targetRec.features ##for debug
                        mask_str=map(lambda f: str(f.location.start.position+1) + "," + str(f.location.end.position + 1) ,exclude_feat)
                        #mask_str=map(lambda f: str(f.location).strip('[]'),exclude_feat)
                        p3_exclude_str = str(mask_str).replace('\', \'',':')
                        p3_target = str(target_feat.location.start.position +1)  + "," + str(target_feat.location.end.position +1)
                        #write sequence record into template file as  fasta
                        t_output_handle = open(tempfastaFile, "w")
                        SeqIO.write([targetRec], t_output_handle, "fasta")
                        t_output_handle.close()
                        #create Primer3Commandline() for emboss tool
                        primer_cl = Primer3Commandline()
                        #set the emboss tool to suppress  output as this will make Galaxy  think it is error message although it is a message to state run success
                        primer_cl.set_parameter("-auto",'1')
                        ###pass  sequence file to emboss. FAILS for primer3_core V2 ##
                        primer_cl.set_parameter("-sequence",tempfastaFile)
                        #add target location
                        primer_cl.set_parameter("-target", p3_target)
                        ##mask off other features...dumb masking of everything at present, beware
                        if (p3_exclude_str != ""):
                            primer_cl.set_parameter("-excludedregion", p3_exclude_str)
                        #add temporary output file to get the result
                        primer_cl.set_parameter("-outfile", tempproutfile)
                        #specify maximum different of tm
                        primer_cl.set_parameter("-maxdifftm",max_tm_diff )
                        #other useful parameters
                        primer_cl.set_parameter("-ogcpercent", opt_GC_percent)
                        primer_cl.set_parameter("-opolyxmax", maxpolyx)  
                        primer_cl.set_parameter("-osize", optimum_length)
                        primer_cl.set_parameter("-otm", optimum_tm)
                        primer_cl.set_parameter("-gcclamp", gc_clamp) 
                        #set product size range
                        primer_cl.set_parameter("-prange", productsizerange)
                        #using python subprocess method to run emboss command line programe with the parameters given
                        fnull = open(os.devnull, 'w')
                        result=subprocess.check_call(str(primer_cl),shell=True ,stdout = fnull, stderr = fnull)
                        #read temporary outputfile
                        handle = open(tempproutfile)
                        record = Primer3.read(handle)

                        if len(record.primers) > 0:
                            primer = record.primers[:5] # We want 5 primer results per target
                            for index, p in enumerate(primer):
                            	# Start postition of feature i.e. where the SNP is, to know where to insert mt_base for mt_amplicon
                            	startposition = [f for f in targetRec.features if f.id==target_ID][0].location.start.position 
                            	# The SNP wt/mu base                       
                              	wtbase = [f for f in targetRec.features if f.id==target_ID][0].qualifiers.get('Reference_seq')
                            	mubase = [f for f in targetRec.features if f.id==target_ID][0].qualifiers.get('Variant_seq') 
                            	wt_base = ''.join(c for c in wtbase if c not in '[]') # Removes the [''] from wtbase
                            	mu_base = ''.join(c for c in mubase if c not in '[]') # Removes the [''] from mubase
                            	# Wild-type amplicon
# NOTE: Due to the script's 1 and 0 indexing inconsistencies, the amplicon starts 1 bp later than the primer forward start, hence the "forward_start-1" required
# ADDENDUM: record.primers[0].reverse_start indicates the primer start FROM LEFT TO RIGHT (3' -> 5') hence reverse_start is actually the END position of the reverse primer if reading from 5' -> 3' and therefore not a good basis for determining exactly where to splice for amplicon, hence the following changes:
                            	wt_amplicon = targetRec.seq[record.primers[index].forward_start-1:record.primers[index].forward_start-1+record.primers[index].size]
                            	# Mutant amplicon
                            	mutantSeq = targetRec.seq.tomutable()
                               	mutantSeq[startposition] = mu_base # Replaces wt base with mutant SNP to make mutant sequence
                            	mu_amplicon = mutantSeq[record.primers[index].forward_start-1:record.primers[index].forward_start-1+record.primers[index].size] # Mutant amplicon slice from mutant seq
                            	# Melt data and respective plots for wt and mu amplicons
                            	wt_melt = getmelt(wt_amplicon)
                            	mu_melt = getmelt(mu_amplicon)
                            	wt_Tm = getTm(wt_melt)
                            	mu_Tm = getTm(mu_melt)
                            	Tm_diff = abs(wt_Tm - mu_Tm)
                            	graphplot(wt_melt)
                            	graphplot(mu_melt)
                            	#pdfpage.savefig()
                            	plt.figure()
                            	#Have plt.show(), plt.figure() or pdfpage.savefig() here to compare WT and mutant amplicon pairs, per primerset variant for each target ID

                            	# Output
                            	outputstr = [mytarget.id, p.forward_seq,p.reverse_seq,p.size, wt_amplicon, mu_amplicon, wt_Tm, mu_Tm, Tm_diff]
                            	print ('\t'.join(map(str,outputstr)))

                            	# Function for writing .txt output files
                            	def maketxt (data,filename):
                            		amptxt = open("%s.txt" % filename , "a" )
                            		amptxt.write(str(data) + '\n')
                            		amptxt.close()
                            	#maketxt(wt_amplicon, "wt_amplicon")
                            	#maketxt(mu_amplicon, "mu_amplicon")
                            	#maketxt(('\t'.join(map(str,outputstr))), "outputstr")
                        
                        else:
                            outputstr = [mytarget.id,"NONE","NONE","NONE", "NONE","NONE","NONE","NONE", "NONE"]
                        print('\t'.join(map(str,outputstr)))
                        #Have plt.show(), plt.figure() or pdfpage.savefig() here to compare WT and mutant amplicon pairs, per primerset variant for each target ID

in_gff_handle.close()
in_seq_handle.close()
pdfpage.close() # Very important, otherwise .pdf will not open properly
