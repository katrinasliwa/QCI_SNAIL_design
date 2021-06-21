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
alpha = 0.09
l_j = 7.0       # big junction length
w_bridge = 0.500 # dolan bridge width
w_btwn = 2.60 #2.120  # width between bridges
n_j = 3         # number of juncitons in array
loopSize = (5.0, 7.750)
w_leadS = 1.0
underCut = 1.100
Mnew = 10          # number of snails in array
Mold = 20
w_lead = 3.0    # width of lead from edge of snails
overlap = 3.5   # overlap between snail and lumped cap (for EBPG stitching errors)

l_tot =  228-4   # Vlad's M=20 array length is 232.0 Lenght of snail array 
oldl_tot=l_tot


#l_tot = 2000.0 # M=200 array length
#l_tot = l_C + g_C + overlap

subx = 34.0*mm #substrate dimensions
suby = 20.0*mm


labelshiftx = 1000.0*um #distances in from the upper right corner for labels
labelshifty = 600.0*um

#rotated chip dimensions 
chipx = 9.0*mm #chip dimensions
chipy = 11.0*mm

dicex = 250.0*um #alignment mark dimensions
dicey = 100.0*um

labelshiftx = 1000.0*um #distances in from the upper right corner for labels
labelshifty = 600.0*um

pump_gap = 50.0*um # gap for pump capacitor
res_width = 250.0*um #width of resonators
SPAshift = -2.8*mm #-2.8

#First Resonator 4 um gap, 9 GHz in middle of band SPA 2
Taperlength = 300.0*um #lengh of the triangular taper 
Res1length = 790*um #length of the resonator segement between the end of the taper and the cap 
Fingergap =4.0*um
Fingerlength = 60.0*um
Numfingers=7 #Total number of fingers in the coupling capacitor. Must be odd. 
Sig_gap = 4.0*um
StickX = 800.0*um
dice_offset = 200.0*um


#Second Resonator 2 um gap, 9 GHz at kerr-free point SPA1
Taperlength_a = 300.0*um #lengh of the triangular taper 
Res1length_a = 650*um #length of the resonator segement between the end of the taper and the cap 
Fingergap_a =2.0*um
Fingerlength_a = 50.0*um
Numfingers_a=5 #Total number of fingers in the coupling capacitor. Must be odd. 
Sig_gap_a = 4.0*um
StickX_a = 800.0*um
dice_offset_a = 200.0*um


#2*Taperlength+Fingerlength+2*Fingergap+2*Res1length+l_tot


#--------Script params---------------------------
design_name = 'SPA_201001' #this is the name of the gds file that gets output
DoseTestBool = 0
TestStructuresBool = 0
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




    
for kk in range(len(opens)):
    #oldcell_name_SnailArray_real = 'SPA_a%.2f_v2_M%d_w%d_uc%d_Ej%.1f'%(alpha, Mold, w*1e3, underCut*1e3, l_j)
    #oldcell_name_SnailArray = oldcell_name_SnailArray_real + labels[kk]
    #oldcellSnailArray = gdspy.Cell(oldcell_name_SnailArray)
    oldsnailArray = make_snail_array(Mold, l_tot, w_lead, alpha, n_j, l_j, loopSize,
                                  w_bridge, w_btwn, underCut, w_leadS,
                                  shorts=shorts[kk], opens=opens[kk],
                                  center=(0.0, 0.0), rotation=(np.pi/2, (0,0)))
    #oldcellSnailArray.add(oldsnailArray)  


kk=0     
Snail_array = make_snail_array(Mold, l_tot, w_lead, alpha, n_j, l_j, loopSize,
                                      w_bridge, w_btwn, underCut, w_leadS,
                                      shorts=shorts[kk], opens=opens[kk],
                                      center=(0.0, 0.0), rotation=(np.pi/2, (0,0)))
kk=1
Snail_array_smallJJ = make_snail_array(Mold, l_tot, w_lead, alpha, n_j, l_j, loopSize,
                                      w_bridge, w_btwn, underCut, w_leadS,
                                      shorts=shorts[kk], opens=opens[kk],
                                      center=(0.0, 0.0), rotation=(np.pi/2, (0,0)))
kk=2
Snail_array_bigJJ = make_snail_array(Mold, l_tot, w_lead, alpha, n_j, l_j, loopSize,
                                      w_bridge, w_btwn, underCut, w_leadS,
                                      shorts=shorts[kk], opens=opens[kk],
                                      center=(0.0, 0.0), rotation=(np.pi/2, (0,0)))
kk=3
Snail_array_short =make_snail_array(Mold, l_tot, w_lead, alpha, n_j, l_j, loopSize,
                                      w_bridge, w_btwn, underCut, w_leadS,
                                      shorts=shorts[kk], opens=opens[kk],
                                      center=(0.0, 0.0), rotation=(np.pi/2, (0,0)))

Snail_array_cell=gdspy.Cell('Snail_array')
Snail_array_cell.add(Snail_array)

Snail_array_cell_sm=gdspy.Cell('Snail_array_smallJJ')
Snail_array_cell_sm.add(Snail_array_smallJJ)

Snail_array_cell_lg=gdspy.Cell('Snail_array_bigJJ')
Snail_array_cell_lg.add(Snail_array_bigJJ)

Snail_array_cell_short=gdspy.Cell('Snail_array_short')
Snail_array_cell_short.add(Snail_array_short)


############### DRAW Resonator 1  ############################
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
resCell.add(gdspy.Polygon([(-taperLedge,res_width/2),(-taperLedge,-res_width/2),(-l_tot/2+overlap,0)], layer = layers['coarse_layer'])) 
    
#Right side taper 
resCell.add(gdspy.Polygon([(taperLedge,res_width/2),(taperLedge,-res_width/2),(l_tot/2-overlap,0)], layer = layers['coarse_layer']))  
    
#resonator R
resCell.add(gdspy.Rectangle((taperLedge,-res_width/2), (taperLedge+Res1length, res_width/2), layer = layers['coarse_layer'])) 
    
#pump bond pad 

resCell.add(gdspy.Rectangle((taperLedge+Res1length+pump_gap,-res_width/2), (chipx/2, res_width/2), layer = layers['coarse_layer'])) 

############### DRAW Resonator 2  ############################
resNamea = res_name = 'Resonatorlength_{0}_nF{1}_Fl{2}_Fg{3}'.format(Res1length_a, Numfingers_a, Fingerlength_a, Fingergap_a)
resCella = gdspy.Cell(resNamea) #should change to include useful info length and cc 

resCella.add(gdspy.Rectangle((-chipx/2, -chipy/2), (chipx/2, chipy/2), layer = layers['chipbounds_layer'])) 

#signal bond pad 
res2edge=l_tot/2+Taperlength_a+Res1length_a+Fingerlength_a+Fingergap_a
resCella.add(gdspy.Rectangle((-chipx/2, -res_width/2), (-res2edge, res_width/2), layer = layers['coarse_layer'])) 

#finger cap array 
Fingerwidth_a=(res_width-(Numfingers_a-1)*Fingergap_a)/Numfingers_a
finger_right=-res2edge+Fingerlength_a+Fingergap_a
for n in range(0,round((Numfingers_a-1)/2),1):
    print(n)
    resCella.add(gdspy.Rectangle((-res2edge, res_width/2-2*n*(Fingerwidth_a+Fingergap_a)), (-res2edge+Fingerlength_a, res_width/2-2*n*(Fingerwidth_a+Fingergap_a)-Fingerwidth_a), layer = layers['coarse_layer'])) #left finger 
    resCella.add(gdspy.Rectangle((-res2edge+Fingergap_a, res_width/2-(2*n+1)*(Fingerwidth_a+Fingergap_a)), (finger_right, res_width/2-(2*n+1)*(Fingerwidth_a+Fingergap_a)-Fingerwidth_a), layer = layers['coarse_layer'])) #right finger 
resCella.add(gdspy.Rectangle((-res2edge, -res_width/2+Fingerwidth_a), (-res2edge+Fingerlength_a, -res_width/2), layer = layers['coarse_layer'])) #last figner  
    
#resonator L
res1edge_a=l_tot/2+Taperlength_a+Res1length_a
resCella.add(gdspy.Rectangle((-res1edge_a, -res_width/2), (-res1edge_a+Res1length_a, res_width/2), layer = layers['coarse_layer'])) 
    
#Left side taper
taperLedge=l_tot/2+Taperlength
resCella.add(gdspy.Polygon([(-taperLedge,res_width/2),(-taperLedge,-res_width/2),(-l_tot/2+overlap,0)], layer = layers['coarse_layer'])) 
    
#Right side taper 
resCella.add(gdspy.Polygon([(taperLedge,res_width/2),(taperLedge,-res_width/2),(l_tot/2-overlap,0)], layer = layers['coarse_layer']))  
    
#resonator R
resCella.add(gdspy.Rectangle((taperLedge,-res_width/2), (taperLedge+Res1length_a, res_width/2), layer = layers['coarse_layer'])) 
    
#pump bond pad 

resCella.add(gdspy.Rectangle((taperLedge+Res1length_a+pump_gap,-res_width/2), (chipx/2, res_width/2), layer = layers['coarse_layer']))
 

###########################Put together the top cell#########################################################
#gdspy.LayoutViewer()

#topCell = gdspy.Cell('topCell')

#make_wafer(wafer_size=75*mm,flat_length_primary=22*mm,flat_length_secondary=0,layer=100,cell=topCell)



#sets up individual chip cells 
if SPABool:    
    substrate = gdspy.Rectangle((-subx/2, -suby/2), (subx/2, suby/2), layer = layers['substrate_layer'])
    #topCell.add(substrate)
    
    chipCell = gdspy.Cell('chipCell')
    #Chip bounds 
    chipCell.add(gdspy.Rectangle((-chipx/2, -chipy/2), (chipx/2, chipy/2), layer = layers['chipbounds_layer'])) 
    #dicing marks 
    chipCell.add(gdspy.Rectangle((-dicex/2, -chipy/2+dice_offset), (dicex/2, -chipy/2+dice_offset+dicey), layer = layers['dice_layer']))
    chipCell.add(gdspy.Rectangle((-dicex/2, chipy/2-dice_offset), (dicex/2, chipy/2-dice_offset-dicey), layer = layers['dice_layer']))
    
    chipCell.add(resCell) 
    chipCell.add(Snail_array_cell)
    
resCell.add(gdspy.Rectangle((-dicex/2, -chipy/2+dice_offset), (dicex/2, -chipy/2+dice_offset+dicey), layer = layers['dice_layer']))
resCell.add(gdspy.Rectangle((-dicex/2, chipy/2-dice_offset), (dicex/2, chipy/2-dice_offset-dicey), layer = layers['dice_layer']))

resCella.add(gdspy.Rectangle((-dicex/2, -chipy/2+dice_offset), (dicex/2, -chipy/2+dice_offset+dicey), layer = layers['dice_layer']))
resCella.add(gdspy.Rectangle((-dicex/2, chipy/2-dice_offset), (dicex/2, chipy/2-dice_offset-dicey), layer = layers['dice_layer']))
    


    

####################### Draw test structures ############################

tsx = 1.5*mm
tsy = 1.0*mm




taper_l = 300*um
extra_l = 85*um
# array_xs = gdspy.Cell(cell_name_SnailArray_real).get_bounding_box()[:, 0]
#newpadPts = [(l_tot/2-overlap, w_lead/2),(l_tot/2-overlap+taper_l, res_width/2),(l_tot/2-overlap+taper_l+extra_l, res_width/2),
#			(l_tot/2-overlap+taper_l+extra_l, -res_width/2),(l_tot/2-overlap+taper_l, -res_width/2),(l_tot/2-overlap, -w_lead/2),]
oldpadPts = [(oldl_tot/2-overlap, w_lead/2),(oldl_tot/2-overlap+taper_l, res_width/2),(oldl_tot/2-overlap+taper_l+extra_l, res_width/2),
    			(oldl_tot/2-overlap+taper_l+extra_l, -res_width/2),(oldl_tot/2-overlap+taper_l, -res_width/2),(oldl_tot/2-overlap, -w_lead/2),]
    
oldpadCell = gdspy.Cell('oldpadCell')
oldpadCell.add(gdspy.Polygon(oldpadPts, layer = layers['coarse_layer']))
oldpadCell.add(gdspy.Polygon(oldpadPts, layer = layers['coarse_layer']).rotate(np.pi))
 



testCell=gdspy.Cell('Snail_array_sm_JJ_cell')
testCell.add(Snail_array_cell_sm)
testCell.add(oldpadCell)



####################### Draw dose tests ############################
if DoseTestBool:
    DT_locs = np.array([(11.5*mm, 7.5*mm)])
    
    d1reps = 10
    d2reps = 6
    
    DTx = 1.2*mm
    DTy = 0.7*mm
    
    DTCell = gdspy.Cell('DTCell')
    
    for w in w_btwn_array:
        for dose1 in np.arange(d1reps):
            for dose2 in np.arange(d2reps):
                    curloc = [DTx*(dose1-(np.float(d1reps)-1)/2), DTy*(dose2-(np.float(d2reps)-1)/2)]
                    snailArray = make_snail_array(Mold, oldl_tot, w_lead, alpha, n_j, l_j, loopSize, w_bridge, w, underCut, w_leadS, 
                        center=curloc, rotation=(np.pi/2, curloc),RES_LAYER = 30+dose2, UC_LAYER = layers['dolan_undercut_layer'],BRIDGE_LAYER = 20+dose1)
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
################################Rotating cells to easily pull them out##################

resCell.copy(name='Resonatorlength_rot_{0}_nF{1}_Fl{2}_Fg{3}'.format(Res1length, Numfingers, Fingerlength, Fingergap),deep_copy=False,translation=None,rotation=np.pi/2)
resCella.copy(name='Resonatorlength_rot_{0}_nF{1}_Fl{2}_Fg{3}'.format(Res1length_a, Numfingers_a, Fingerlength_a, Fingergap_a),deep_copy=False,translation=None,rotation=np.pi/2)
#cellSnail.copy(name='SNAIL_rot',deep_copy=False,translation=None,rotation=np.pi/2)
chipCell.copy(name='Chip_cell_rot',deep_copy=False,translation=None,rotation=np.pi/2)

oldpadCell.copy(name='SNAIL_pads_rot',deep_copy=False,translation=None,rotation=np.pi/2)
#oldcellSnailArray.copy(name='SNAIL_array_rot',deep_copy=False,translation=None,rotation=np.pi/2)

Snail_array_cell.copy(name='SNAIL_rot',deep_copy=False,translation=None,rotation=np.pi/2)
Snail_array_cell_sm.copy(name='SNAIL_rot_test_sm',deep_copy=False,translation=None,rotation=np.pi/2)
Snail_array_cell_lg.copy(name='SNAIL_rot_test_lg',deep_copy=False,translation=None,rotation=np.pi/2)
Snail_array_cell_short.copy(name='SNAIL_rot_test_short',deep_copy=False,translation=None,rotation=np.pi/2)



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
				])*scale #doses for the real SNAIL devices as defined in the top cell 

print('device doses 5 nA: ' + np.array_str(base_5nA*device_doses_5nA, precision = 1)+'\n')


layers_5nA = device_layers_5nA
doses_5nA = device_doses_5nA
layers_50nA = device_layers_50nA
doses_50nA = device_doses_50nA

if DoseTestBool:
   
    #layers and doses for the dose test 
    DT_fine_layers = 30+np.arange(d2reps)
    DT_bridge_layers = 20+np.arange(d1reps)
    
    DT_bridge_doses = 1.0/base_5nA*np.linspace(100, 300, d1reps)*scale
    DT_fine_doses = 1.0/base_5nA*np.linspace(600, 1300, d2reps)*scale
    
    print('DT_bridge_doses: ' + np.array_str(base_5nA*DT_bridge_doses, precision = 1)+'\n')
    print('DT_fine_doses: ' + np.array_str(base_5nA*DT_fine_doses, precision = 1)+'\n')
    layers_5nA = np.concatenate((layers_5nA, DT_bridge_layers,DT_fine_layers))
    doses_5nA = np.concatenate((doses_5nA, DT_bridge_doses, DT_fine_doses))
    #print 'DT_coarse_doses: ' + np.array_str(base_50nA*DT_coarse_doses, precision = 1)+'\n'
    
    # print 'hi'
    #layers_50nA = np.concatenate((device_layers_50nA, DT_coarse_layers))
    #layers_50nA = np.concatenate((device_layers_50nA, DT_coarse_layers))
    #doses_50nA = np.concatenate((device_doses_50nA, DT_coarse_doses))
    


np.savetxt('doses50nA.txt', list(zip(layers_50nA, doses_50nA)))
np.savetxt('doses5nA.txt', list(zip(layers_5nA, doses_5nA)))

#fmt='%u(0), %.4f'




