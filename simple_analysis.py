import sys 

sys.path.insert(0,"/eos/user/j/jstowell/Software/Analysis/hk-BONSAI/wcte_python_toolkit/python/")
import pyBONSAI

from pyWCTEAnalysisTools import *

import ROOT 
import glob 
import numpy as np 
import pandas as pd 
import array
import pandas as pd
import os
import cppyy

pyBONSAI.set_default_directories(
    WCSIM_BUILD_DIR="/eos/user/j/jstowell/Software/WCSim/WCSim-develop/build4/",
    BONSAIDIR="/eos/user/j/jstowell/Software/Analysis/hk-BONSAI/wcte_python_toolkit/",
    BONSAIPARAM="/eos/user/j/jstowell/Software/Analysis/hk-BONSAI/wcte_python_toolkit/data/fit_param_wcte.dat"
)
bonsai = pyBONSAI.WCTE_hkBONSAI()


       
    
# Load the geometry mapping into a DataFrame
runno = int(sys.argv[1])
tt = get_offline_run_tchain(runno)
    
# Start Filter
ns = 1
count = 0
nhit_cut = 10
mvwindow_start =   10000*ns
mvwindow_end   = 5000000*ns
mvwindow_step  =      10*ns
mvwindow_width =      50*ns

vertex = {
"nhits": [],
"nhitso": [],
"x": [],
"y": [],
"z": [],
"result0": [],
"result1": [],
"result2": [],
"result3": [],
"result4": [],
"result5": [],
"good0":[],
"good1":[],
"good2":[] 
}

hitslist = {
        "x":[],
        "y":[],
        "z":[],
        "q":[],
        "t":[]
}

for e in tt:
    count += 1    
    if count > 10: break
    print(count)
    
    
    # Generate hit_collection from event
    hits = hit_collection(e)
    #     hits = mchit_collection(e)
    
    # Run a moving average 50ns hit filter
    tstart = mvwindow_start 
    tend   = tstart + mvwindow_width
    tmax   = np.max(hits.t)
    
    while tend < tmax:
        
        window = hits.time_slice(tstart, tend)
        tstart += mvwindow_step
        tend   = tstart + mvwindow_width
    
        # Run prehit filter for this window
        if len(window.t) < nhit_cut: continue
        if len(window.t) > 400: continue
            
                
        # Run Bonsai
        bsVertex = array.array('f',3*[0.0])
        bsResult = array.array('f',6*[0.0])
        bsGood = array.array('f',3*[0.0])
        bsNhit = array.array('i',[len(window.cable)])
        bsNsel = array.array('i',[0])

        # Generate hit collection for this triggger
        bsCAB_a = array.array('i', window.cable)
        bsT_a = array.array('f', window.t - np.min(window.t) + 200)
        bsQ_a = array.array('f', window.q)

        # Run Bonsai
        try:
            nhits = bonsai.BonsaiFit(bsVertex, bsResult, bsGood, bsNsel, bsNhit, bsCAB_a, bsT_a, bsQ_a);
        except:
            print("BONSAIFAILED");
            pass
#         print(nhits, bsVertex, bsResult, bsGood, bsNsel, bsNhit)

        vertex["nhits"].append(nhits)
        vertex["nhitso"].append(len(window.t))
        
        vertex["x"].append(bsVertex[0])
        vertex["y"].append(bsVertex[1])
        vertex["z"].append(bsVertex[2])
        vertex["result0"].append(bsResult[0])
        vertex["result1"].append(bsResult[1])
        vertex["result2"].append(bsResult[2])
        vertex["result3"].append(bsResult[3])
        vertex["result4"].append(bsResult[4])
        vertex["result5"].append(bsResult[5])
        vertex["good0"].append(bsGood[0])
        vertex["good1"].append(bsGood[1])
        vertex["good2"].append(bsGood[2])
        
        # Skip to next window if we found a hit
        tstart += mvwindow_width
        tend    = tstart + mvwindow_width
        
        if bsResult[1] > 0.9:
            for j in range(len(window.t)):
                hitslist["x"].append( window.x[j] )
                hitslist["y"].append( window.y[j] )
                hitslist["z"].append( window.z[j] )
                hitslist["q"].append( window.q[j] )
                hitslist["t"].append( window.t[j] )
                
        
for key in vertex:
    print(key, len(vertex[key]))
    
dfhits = pd.DataFrame(hitslist)
dfhits.to_csv(f"hitslow_{runno}.csv")

df = pd.DataFrame(vertex)
df.to_csv(f"vertexlow_{runno}.csv")
