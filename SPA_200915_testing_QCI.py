import sys
#sys.path.append('/Users/yiwenchu/Documents/circuitQED/CAD/SNAIL')
#sys.path.append('/Users/yiwenchu/Documents/python/gdspy-0.6/build/lib.macosx-10.9-x86_64-2.7/gdspy')
import os
import numpy as np
import gdspy
# import waferDefs as wd
# from toolsDefs import *
# import qubitDefs as qd
# import junctionDefs as jd
# import diceDefs as dd
# import markerDefs as md
from SNAILs import *

print('Using gdspy module version ' + gdspy.__version__)

layers = {'chipbounds_layer':0,
		  'dice_layer':12,
		  'coarse_layer':12,
          'extra_coarse_layer':13,
		  'dolan_fine_layer':1,
		  'dolan_undercut_layer':2,
		  'dolan_bridge_layer':3,
		  'substrate_layer':10,
		  'fields_layer':11,};

nm = 0.001;
mm = 1000.0;
um = 1.0;
mil = 25.4;
inches = mil*1000;


# ------- Snail Array-------------------------
alpha = 0.1
l_j = 7.0       # big junction length
w_bridge = 0.500 # dolan bridge width
w_btwn = 2.200 #2.120  # width between bridges
n_j = 3         # number of juncitons in array
loopSize = (5.0, 7.750)
w_leadS = 1.0
underCut = 1.100
Mnew = 10          # number of snails in array
Mold = 20
w_lead = 3.0    # width of lead from edge of snails
overlap = 2.0   # overlap between snail and lumped cap (for EBPG stitching errors)
#l_tot = 70.0
l_tot =  134.5   # Vlad's M=20 array length is 232.0 Lenght of snail array 
oldl_tot = 232.0
#l_tot = 2000.0 # M=200 array length
#l_tot = l_C + g_C + overlap

subx = 34.0*mm #substrate dimensions
suby = 20.0*mm
chipx = 11.0*mm #chip dimensions
chipy = 5.0*mm
dicex = 250.0*um #alignment mark dimensions
dicey = 100.0*um
pump_gap = 50.0*um # gap for pump capacitor
res_width = 250.0*um #width of resonators
SPAshift = -2.8*mm #-2.8

labelshiftx = 1000.0*um #distances in from the upper right corner for labels
labelshifty = 600.0*um

Taperlength = 300.0*um #lengh of the triangular taper 
Res1length = 1000.0*um #length of the resonator segement between the end of the taper and the cap 

Fingergap = 2.0*um
Fingerlength = 70.0*um
Numfingers=7 #Total number of fingers in the coupling capacitor. Must be odd. 

Sig_gap = 4.0*um
StickX = 800.0*um
dice_offset = 200.0*um




#--------Script params---------------------------
design_name = 'SPA_200915_QCI_playing_2' #this is the name of the gds file that gets output
DoseTestBool = 1
TestStructuresBool = 1
SPABool = 1

#-----------Start drawing-----------------------------
# a single SNAIL cell
cellSnail = gdspy.Cell('SNAIL')
snail = make_snail(alpha, n_j, l_j, loopSize, w_bridge, w_btwn, underCut,
                   w_leadS, edgeSnail=(True, True), topLeftLead=(0.0, 0.0),
                   shorts=(False,False), opens=(False,False),RES_LAYER = layers['dolan_fine_layer'], UC_LAYER = layers['dolan_undercut_layer'], BRIDGE_LAYER = layers['dolan_bridge_layer']
                   )

cellSnail.add(snail)
#RES_LAYER = layers['dolan_fine_layer'], UC_LAYER = layers['dolan_undercut_layer'], BRIDGE_LAYER = layers['dolan_bridge_layer']


shorts = [(False,False), (False,False), (False, False), (True,True)]
opens = [(False, False), (False, True), (True, False), (False, False)]
labels = ['', '_test1', '_test2', '_short']



#    width sweep M=20
w_btwn_array = np.array([2.200])
for w in w_btwn_array:
    for kk in range(len(opens)):
        newcell_name_SnailArray_real = 'SPA_a%.2f_v2_M%d_w%d_uc%d_Ej%.1f'%(alpha, Mnew, w*1e3, underCut*1e3, l_j)
        newcell_name_SnailArray = newcell_name_SnailArray_real + labels[kk]
        newcellSnailArray = gdspy.Cell(newcell_name_SnailArray)
        newsnailArray = make_snail_array(Mnew, l_tot, w_lead, alpha, n_j, l_j, loopSize,
                                      w_bridge, w, underCut, w_leadS,
                                      shorts=shorts[kk], opens=opens[kk],
                                      center=(0.0, 0.0), rotation=(np.pi/2, (0,0)),RES_LAYER = layers['dolan_fine_layer'], UC_LAYER = layers['dolan_undercut_layer'], BRIDGE_LAYER = layers['dolan_bridge_layer'])
        newcellSnailArray.add(newsnailArray)    
"""        
    for kk in range(len(opens)):
        oldcell_name_SnailArray_real = 'SPA_a%.2f_v2_M%d_w%d_uc%d_Ej%.1f'%(alpha, Mold, w*1e3, underCut*1e3, l_j)
        oldcell_name_SnailArray = oldcell_name_SnailArray_real + labels[kk]
        oldcellSnailArray = gdspy.Cell(oldcell_name_SnailArray)
        oldsnailArray = make_snail_array(Mold, oldl_tot, w_lead, alpha, n_j, l_j, loopSize,
                                      w_bridge, w, underCut, w_leadS,
                                      shorts=shorts[kk], opens=opens[kk],
                                      center=(0.0, 0.0), rotation=(np.pi/2, (0,0)))
        oldcellSnailArray.add(oldsnailArray)    

"""
############### DRAW Resonator  ############################
resName = res_name = 'Resonatorlength_{0}_nF{1}_Fl{2}_Fg{3}'.format(Res1length, Numfingers, Fingerlength, Fingergap)
resCell = gdspy.Cell(resName) #should change to include useful info length and cc 

resCell.add(gdspy.Rectangle((-chipx/2, -chipy/2), (chipx/2, chipy/2), layer = layers['chipbounds_layer'])) 

#signal bond pad 
res2edge=l_tot/2+Taperlength+Res1length+Fingerlength+Fingergap
resCell.add(gdspy.Rectangle((-chipx/2, -res_width/2), (-res2edge, res_width/2), layer = layers['coarse_layer'])) 

#finger cap array 
Fingerwidth=(res_width-(Numfingers-1)*Fingergap)/Numfingers
finger_right=-res2edge+Fingerlength+Fingergap
for n in range(0,round((Numfingers-1)/2),1):
    print(n)
    resCell.add(gdspy.Rectangle((-res2edge, res_width/2-2*n*(Fingerwidth+Fingergap)), (-res2edge+Fingerlength, res_width/2-2*n*(Fingerwidth+Fingergap)-Fingerwidth), layer = layers['coarse_layer'])) #left finger 
    resCell.add(gdspy.Rectangle((-res2edge+Fingergap, res_width/2-(2*n+1)*(Fingerwidth+Fingergap)), (finger_right, res_width/2-(2*n+1)*(Fingerwidth+Fingergap)-Fingerwidth), layer = layers['coarse_layer'])) #right finger 
resCell.add(gdspy.Rectangle((-res2edge, -res_width/2+Fingerwidth), (-res2edge+Fingerlength, -res_width/2), layer = layers['coarse_layer'])) #last figner  
    
#resonator L
res1edge=l_tot/2+Taperlength+Res1length
resCell.add(gdspy.Rectangle((-res1edge, -res_width/2), (-res1edge+Res1length, res_width/2), layer = layers['coarse_layer'])) 
    
#Left side taper
taperLedge=l_tot/2+Taperlength
resCell.add(gdspy.Polygon([(-taperLedge,res_width/2),(-taperLedge,-res_width/2),(-l_tot/2,0)], layer = layers['coarse_layer'])) 
    
#Right side taper 
resCell.add(gdspy.Polygon([(taperLedge,res_width/2),(taperLedge,-res_width/2),(l_tot/2,0)], layer = layers['coarse_layer']))  
    
#resonator R
resCell.add(gdspy.Rectangle((taperLedge,-res_width/2), (taperLedge+Res1length, res_width/2), layer = layers['coarse_layer'])) 
    
#pump bond pad 

resCell.add(gdspy.Rectangle((taperLedge+Res1length+pump_gap,-res_width/2), (chipx/2, res_width/2), layer = layers['coarse_layer'])) 
 
    

####################################################################################
#gdspy.LayoutViewer()
topCell = gdspy.Cell('topCell')





if SPABool:    
    substrate = gdspy.Rectangle((-subx/2, -suby/2), (subx/2, suby/2), layer = layers['substrate_layer'])
    topCell.add(substrate)
    
    chipCell = gdspy.Cell('chipCell')
    #Chip bounds 
    chipCell.add(gdspy.Rectangle((-chipx/2, -chipy/2), (chipx/2, chipy/2), layer = layers['chipbounds_layer'])) 
    #dicing marks 
    chipCell.add(gdspy.Rectangle((-dicex/2, -chipy/2+dice_offset), (dicex/2, -chipy/2+dice_offset+dicey), layer = layers['dice_layer']))
    chipCell.add(gdspy.Rectangle((-dicex/2, chipy/2-dice_offset), (dicex/2, chipy/2-dice_offset-dicey), layer = layers['dice_layer']))
    
    locx = np.linspace(-chipx, chipx, 3)+0.2*mm #i'm not sure what this is doing 
    locy = np.linspace(-1.5*chipy, 1.5*chipy, 4)
    SPA_locs = np.reshape([[(x, y) for y in locy] for x in locx], (len(locx)*len(locy),2)) 
    #print SPA_locs
    #SPA1_locs = SPA_locs[0::2]#[1:] the commented command removes the first element from the list
    #SPA2_locs = SPA_locs[1::2]
    
    #SPA1_locs = SPA_locs[[0,4]]#[1:] the commented command removes the first element from the list
    #SPA2_locs = SPA_locs[1]
    #SPA3_locs = SPA_locs[2]
    SPA4_locs = SPA_locs[0:11]#SPA_locs[[3,5]]
    # print SPA_locs
    

    

    
    for loc in SPA4_locs:
        print(loc)
        topCell.add(gdspy.CellReference(resCell, loc))
        topCell.add(gdspy.CellReference(oldsnailArray, loc))
        #topCell.add(gdspy.CellReference(oldcell_name_SnailArray_real, np.add(loc, (SPAshift, 0))))
        #print(newcell_name_SnailArray_real)
    
    """
    #draw field edges`
    fieldlength = 600.0*um
    numyfields = np.ceil(suby/fieldlength)   
    numxfields = np.ceil(subx/fieldlength)
    Fieldcell = gdspy.Cell('Fieldcell')
    for xind in np.arange(numxfields):
        fieldx = xind*fieldlength
        for yind in np.arange(numyfields):
            fieldy = yind*fieldlength
            Fieldcell.add(gdspy.Rectangle((fieldx, fieldy), (fieldx + fieldlength, fieldy + fieldlength), layer = layers['fields_layer']))
            
    topCell.add(gdspy.CellReference(Fieldcell, (-subx/2,-suby/2)))
"""    

    

####################### Draw test structures ############################
if TestStructuresBool:
    # test_st_locs = np.array([(2.5*mm, -7.5*mm), (2.5*mm, 7.5*mm)])
    #test_st_locs = np.array([(11.0*mm, -5*mm)])
    M20test_st_locs = np.array([(6.5*mm, 7.5*mm)])
    tsx = 1.5*mm
    tsy = 1.0*mm
    #texts = ['10 real', '10 small', '10 large', '10 short']
    oldtexts = ['20 real', '20 small', '20 large', '20 short']
    
    
    taper_l = 300*um
    extra_l = 85*um
    # array_xs = gdspy.Cell(cell_name_SnailArray_real).get_bounding_box()[:, 0]
    #newpadPts = [(l_tot/2-overlap, w_lead/2),(l_tot/2-overlap+taper_l, res_width/2),(l_tot/2-overlap+taper_l+extra_l, res_width/2),
    #			(l_tot/2-overlap+taper_l+extra_l, -res_width/2),(l_tot/2-overlap+taper_l, -res_width/2),(l_tot/2-overlap, -w_lead/2),]
    oldpadPts = [(oldl_tot/2-overlap, w_lead/2),(oldl_tot/2-overlap+taper_l, res_width/2),(oldl_tot/2-overlap+taper_l+extra_l, res_width/2),
    			(oldl_tot/2-overlap+taper_l+extra_l, -res_width/2),(oldl_tot/2-overlap+taper_l, -res_width/2),(oldl_tot/2-overlap, -w_lead/2),]
    
    '''
    padCell = gdspy.Cell('padCell')
    padCell.add(gdspy.Polygon(newpadPts, layer = layers['coarse_layer']))
    padCell.add(gdspy.Polygon(newpadPts, layer = layers['coarse_layer']).rotate(np.pi))
    '''
    
    oldpadCell = gdspy.Cell('oldpadCell')
    oldpadCell.add(gdspy.Polygon(oldpadPts, layer = layers['coarse_layer']))
    oldpadCell.add(gdspy.Polygon(oldpadPts, layer = layers['coarse_layer']).rotate(np.pi))
    
    testCell = gdspy.Cell('M10testCell')
    oldtestCell = gdspy.Cell('M20testCell')
    # for ind, x in enumerate(np.linspace(-test_space*1.5, test_space*1.5, 4)):
    for ind, x in enumerate(np.array([(0, -3.0*tsy/2.0),(0, -tsy/2),(0, tsy/2),(0,3.0*tsy/2)])):
    	# print cell_name_SnailArray_real + labels[ind]
        '''
    	testCell.add(gdspy.CellReference(newcell_name_SnailArray_real + labels[ind], x))
    	testCell.add(gdspy.CellReference(padCell, x))
    	testCell.add(gdspy.Text(texts[ind], 100*um, position = (x[0]-200*um, x[1]+600*um), layer = layers['coarse_layer']))
        '''
        oldtestCell.add(gdspy.CellReference(oldcell_name_SnailArray_real + labels[ind], x))
        oldtestCell.add(gdspy.CellReference(oldpadCell, x))
        oldtestCell.add(gdspy.Text(oldtexts[ind], 100*um, position = (x[0]-200*um, x[1]+400*um), layer = layers['coarse_layer']))
    
    '''
    for loc in test_st_locs:
    	topCell.add(gdspy.CellReference(testCell, loc))
        '''
        
    for loc in M20test_st_locs:
    	topCell.add(gdspy.CellReference(oldtestCell, loc))
    

####################### Draw dose tests ############################
if DoseTestBool:
    DT_locs = np.array([(11.5*mm, 7.5*mm)])
    
    d1reps = 7
    d2reps = 6
    
    DTx = 1.2*mm
    DTy = 0.7*mm
    
    DTCell = gdspy.Cell('DTCell')
    
    for w in w_btwn_array:
        for dose1 in np.arange(d1reps):
            for dose2 in np.arange(d2reps):
                    curloc = [DTx*(dose1-(np.float(d1reps)-1)/2), DTy*(dose2-(np.float(d2reps)-1)/2)]
                    snailArray = make_snail_array(Mold, oldl_tot, w_lead, alpha, n_j, l_j, loopSize, w_bridge, w, underCut, w_leadS, 
                        center=curloc, rotation=(np.pi/2, curloc))
                    #RES_LAYER = 30+dose2, UC_LAYER = layers['dolan_undercut_layer'], BRIDGE_LAYER = 20+dose1
                    DTCell.add(snailArray)
                    DTCell.add(gdspy.Polygon(np.add(oldpadPts, curloc), layer = layers['coarse_layer']))
                    DTCell.add(gdspy.Polygon(np.add(oldpadPts, curloc), layer  = layers['coarse_layer']).rotate(np.pi, center = curloc))

    for loc in DT_locs:
        topCell.add(gdspy.CellReference(DTCell, loc))
    
    '''
    for w in w_btwn_array:
    	for dose1 in np.arange(d1reps-1):
            curloc = DTy*(dose1-(np.float(d1reps)-1)/2)
            for plgs in capPlgs:
                DTCell.add(gdspy.Polygon(plgs, layer = 20+dose1).translate(1400*um, curloc))
            DTCell.add(gdspy.Rectangle((-300*um+xshift, -res_width/2+curloc), (xshift, res_width/2+curloc), layer = 20+dose1))
            DTCell.add(gdspy.Rectangle((xshift2, -res_width/2+curloc), (300*um+xshift2, res_width/2+curloc), layer = 20+dose1))

            # DTCell.add(gdspy.Rectangle((-300*um-xshift, -res_width/2+curloc), (-xshift, res_width/2+curloc), layer = 20+dose1))
            # DTCell.add(gdspy.Rectangle((2*um-xshift, -res_width/2+curloc), (300*um-xshift, res_width/2+curloc), layer = 20+dose1))
            # DTCell.add(gdspy.Rectangle((-300*um-xshift, -res_width/2+curloc), (-xshift, res_width/2+curloc), layer = 20+dose1))
            # DTCell.add(gdspy.Rectangle((2*um-xshift, -res_width/2+curloc), (300*um-xshift, res_width/2+curloc), layer = 20+dose1))

### generating one test structure for trying out different beamer flow
    curloc = DTy*(d1reps/2)
    for plgs in capPlgs:
        DTCell.add(gdspy.Polygon(plgs, layer = 12).translate(1400*um, curloc))
    DTCell.add(gdspy.Rectangle((-300*um+xshift, -res_width/2+curloc), (xshift, res_width/2+curloc), layer = 12))
    DTCell.add(gdspy.Rectangle((xshift2, -res_width/2+curloc), (300*um+xshift2, res_width/2+curloc), layer = 12))
    
    for loc in DT_locs:
    	topCell.add(gdspy.CellReference(DTCell, loc))
'''



########################## write the gds ##############################
gdspy.write_gds(design_name + '.gds', unit=1.0e-6, precision=1.0e-9)
#gdspy.gds_print(design_name + '.gds', unit=1.0e-6, precision=1.0e-9) #old gdspy version


########################## Output doses ###################################

base_50nA = 500
base_5nA = 150

scale = 1
# 
# device_layers_50nA = [layers['coarse_layer']
# 				]
# device_doses_50nA = np.array([560.0/base_50nA])*scale
# print 'device doses 50 nA: ' + np.array_str(base_50nA*device_doses_50nA, precision = 1)+'\n'

device_layers_50nA = [layers['coarse_layer']]
device_doses_50nA = []
print('device doses 50 nA: '+np.array_str(np.array([750.0]), precision = 1)+'\n')

device_layers_5nA = [layers['dolan_fine_layer'],
				layers['dolan_undercut_layer'],
				layers['dolan_bridge_layer'],
				]
device_doses_5nA = np.array([7.4,
				1.48,
				1.5,
				])*scale

print('device doses 5 nA: ' + np.array_str(base_5nA*device_doses_5nA, precision = 1)+'\n')


layers_5nA = device_layers_5nA
doses_5nA = device_doses_5nA
layers_50nA = device_layers_50nA
doses_50nA = device_doses_50nA

if DoseTestBool:
    #DT_coarse_layers = np.linspace(20, 20+d1reps-2, d1reps-1)
    #DT_coarse_doses = np.linspace(1.0, 1.27, d1reps-1)*scale
    DT_fine_layers = 30+np.arange(d2reps)
    DT_bridge_layers = 20+np.arange(d1reps)
    
    DT_bridge_doses = 1.0/base_5nA*np.linspace(150, 300, d1reps)*scale
    DT_fine_doses = 1.0/base_5nA*np.linspace(900, 1300, d2reps)*scale
    
    print('DT_bridge_doses: ' + np.array_str(base_5nA*DT_bridge_doses, precision = 1)+'\n')
    print('DT_fine_doses: ' + np.array_str(base_5nA*DT_fine_doses, precision = 1)+'\n')
    layers_5nA = np.concatenate((layers_5nA, DT_bridge_layers,DT_fine_layers))
    doses_5nA = np.concatenate((doses_5nA, DT_bridge_doses, DT_fine_doses))
    #print 'DT_coarse_doses: ' + np.array_str(base_50nA*DT_coarse_doses, precision = 1)+'\n'
    
    # print 'hi'
    #layers_50nA = np.concatenate((device_layers_50nA, DT_coarse_layers))
    #layers_50nA = np.concatenate((device_layers_50nA, DT_coarse_layers))
    #doses_50nA = np.concatenate((device_doses_50nA, DT_coarse_doses))
    


np.savetxt('doses50nA.txt', list(zip(layers_50nA, doses_50nA)), fmt='%u(0), %.4f')
np.savetxt('doses5nA.txt', list(zip(layers_5nA, doses_5nA)), fmt='%u(0), %.4f')





