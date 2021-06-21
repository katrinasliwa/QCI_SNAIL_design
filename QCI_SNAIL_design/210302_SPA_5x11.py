import sys
#sys.path.append('/Users/yiwenchu/Documents/circuitQED/CAD/SNAIL')
#sys.path.append('/Users/yiwenchu/Documents/python/gdspy-0.6/build/lib.macosx-10.9-x86_64-2.7/gdspy')
import os
import numpy as np
import gdspy
from layout import *
# import waferDefs as wd
# from toolsDefs import *
# import qubitDefs as qd
# import junctionDefs as jd
# import diceDefs as dd
# import markerDefs as md
from SNAILs import *

print('Using gdspy module version ' + gdspy.__version__)

clearGDSLibrary=1 #allows you to re-run same file in the same console and re-make the file 


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

Rotate=0 #do you want the Snail resonator along the long chip axis (0) or the short chip axis (1)
OutputDoses=1 #Create dose file 


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
overlap = 4.5   # overlap between snail and lumped cap (for EBPG stitching errors)

l_tot =  228-4   # Vlad's M=20 array length is 232.0 Lenght of snail array 
oldl_tot=l_tot


#l_tot = 2000.0 # M=200 array length
#l_tot = l_C + g_C + overlap

subx = 34.0*mm #substrate dimensions
suby = 20.0*mm
wafer_size=74*mm

labelshiftx = 1000.0*um #distances in from the upper right corner for labels
labelshifty = 600.0*um

#rotated chip dimensions 
chipx = 11.0*mm #chip dimensions, note it assumed this is the axis your resonator will be aligned along 
chipy = 5.0*mm



array_nx=int(np.floor(wafer_size/chipx)) #if instead you want to create a fixed nxm array of chips 
array_ny=int(np.floor(wafer_size/chipy))

dicex = 250.0*um #alignment mark dimensions
dicey = 100.0*um


pump_gap = 50.0*um # gap for pump capacitor
res_width = 250.0*um #width of resonators
SPAshift = -2.8*mm #-2.8


#--------------------Resonator Definitions------------------

Taperlength = 300.0*um #lengh of the triangular taper 
Sig_gap = 4.0*um
StickX = 800.0*um
dice_offset = 200.0*um

num_res=5
Res1Lengths=[790*um,560*um,680*um,440*um,560*um]
Numfingers_arry=[7,7,7,7,7]
Fingergaps=[4.0*um,4.0*um,4.0*um,4.0*um,4.0*um]
Fingerlengths=[60.0*um,60.0*um,60.0*um,60.0*um,80.0*um]
#First Resonator 4 um gap, 9 GHz in middle of band SPA 2

w_btwn_array=[w_btwn] #for the dose test

#--------Script params---------------------------
design_name = 'SPA_21302' #this is the name of the gds file that gets output
DoseTestBool = 1
TestStructuresBool = 1
SPABool = 1


standard_layout=1 #do you want to use a standard 3-inch layout for cells? if so 64 standard resonators 

#number of each of the resonators above. Will be tesselated top to bottom then left to right  
res_layout=[12,13,13,13,13]

if standard_layout and np.sum(res_layout)>64:
    print('Warning, requesting {} resonators, too many for layout'.format(np.sum(res_layout)))
    
if standard_layout and np.sum(res_layout)<64:
    print('Warning, requesting {} resonators, too few for layout'.format(np.sum(res_layout)))
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



# Make test arrays 
    
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

if Rotate:
    Snail_array_cell.copy(name='SNAIL_rot',deep_copy=False,translation=None,rotation=np.pi/2)
    Snail_array_cell_sm.copy(name='SNAIL_rot_test_sm',deep_copy=False,translation=None,rotation=np.pi/2)
    Snail_array_cell_lg.copy(name='SNAIL_rot_test_lg',deep_copy=False,translation=None,rotation=np.pi/2)
    Snail_array_cell_short.copy(name='SNAIL_rot_test_short',deep_copy=False,translation=None,rotation=np.pi/2)
    
    



#--------------------- draw individual resonators---------------------------

i=0


if num_res!=np.amax([len(Res1Lengths),len(Numfingers_arry),len(Fingergaps),len(Fingerlengths)]):
    print("Caution! Expecting a different number of resonator than defined")
else:
    print("Creating {0} resonators".format(num_res))
    
res_cells={}

while i<num_res:
    try:
        Res1length=Res1Lengths[i]
    except:
        print("{0} resonator lenght does not exist. Using previous length".format(i))
        Res1length=Res1Lengths[i-1]
        
    try:
        Numfingers=Numfingers_arry[i]
    except:
        print("{0} num fingers does not exist. Using previous length".format(i))
        Numfingers=Numfingers_arry[i-1]
        
    try:
        Fingergap=Fingergaps[i]
    except:
        print("{0} finger gap does not exist. Using previous length".format(i))
        Fingergap=Fingergaps[i-1]
    try:
        Fingerlength=Fingerlengths[i]
    except:
        print("{0} finger length does not exist. Using previous length".format(i))
        Fingerlength=Fingerlengths[i-1]
        
    resName = res_name = 'Resonatorlength_{0}_nF{1}_Fl{2}_Fg{3}'.format(Res1length, Numfingers, Fingerlength, Fingergap)
    print('wring resontor Resonatorlength_{0}_nF{1}_Fl{2}_Fg{3}'.format(Res1length, Numfingers, Fingerlength, Fingergap))
    
    
    resCell = gdspy.Cell(resName) #should change to include useful info length and cc 
    
    resCell.add(gdspy.Rectangle((-chipx/2, -chipy/2), (chipx/2, chipy/2), layer = layers['chipbounds_layer'])) 
    
    #signal bond pad 
    res2edge=l_tot/2+Taperlength+Res1length+Fingerlength+Fingergap
    resCell.add(gdspy.Rectangle((-chipx/2, -res_width/2), (-res2edge, res_width/2), layer = layers['coarse_layer'])) 
    
    #finger cap array 
    Fingerwidth=(res_width-(Numfingers-1)*Fingergap)/Numfingers
    finger_right=-res2edge+Fingerlength+Fingergap
    for n in range(0,round((Numfingers-1)/2),1):
        #print(n)
        resCell.add(gdspy.Rectangle((-res2edge, res_width/2-2*n*(Fingerwidth+Fingergap)), (-res2edge+Fingerlength, res_width/2-2*n*(Fingerwidth+Fingergap)-Fingerwidth), layer = layers['coarse_layer'])) #left finger 
        resCell.add(gdspy.Rectangle((-res2edge+Fingergap, res_width/2-(2*n+1)*(Fingerwidth+Fingergap)), (finger_right, res_width/2-(2*n+1)*(Fingerwidth+Fingergap)-Fingerwidth), layer = layers['coarse_layer'])) #right finger 
    resCell.add(gdspy.Rectangle((-res2edge, -res_width/2+Fingerwidth), (-res2edge+Fingerlength, -res_width/2), layer = layers['coarse_layer'])) #last figner  
        
    #resonator L
    res1edge=l_tot/2+Taperlength+Res1length
    resCell.add(gdspy.Rectangle((-res1edge, -res_width/2), (-res1edge+Res1length, res_width/2), layer = layers['coarse_layer'])) 
        
    #Left side taper
    taperLedge=l_tot/2+Taperlength
    resCell.add(gdspy.Polygon([(-taperLedge,res_width/2),(-taperLedge,-res_width/2),(-l_tot/2+overlap,-1.5),(-l_tot/2+overlap,1.5)], layer = layers['coarse_layer'])) 
        
    #Right side taper 
    resCell.add(gdspy.Polygon([(taperLedge,res_width/2),(taperLedge,-res_width/2),(l_tot/2-overlap,-1.5),(l_tot/2-overlap,1.5)], layer = layers['coarse_layer']))  
        
    #resonator R
    resCell.add(gdspy.Rectangle((taperLedge,-res_width/2), (taperLedge+Res1length, res_width/2), layer = layers['coarse_layer'])) 
        
    #pump bond pad 
    
    resCell.add(gdspy.Rectangle((taperLedge+Res1length+pump_gap,-res_width/2), (chipx/2, res_width/2), layer = layers['coarse_layer'])) 
    
    if Rotate: 
        resCell.copy(name='Resonatorlength_rot_{0}_nF{1}_Fl{2}_Fg{3}'.format(Res1length, Numfingers, Fingerlength, Fingergap),deep_copy=False,translation=None,rotation=np.pi/2)


    res_cells[i]=resCell
    
    i+=1


###########################Put together the a chip cell#########################################################
    
current_res=res_cells[0]
print(type(current_res.name))

#sets up individual chip cells 
if SPABool:    
    substrate = gdspy.Rectangle((-subx/2, -suby/2), (subx/2, suby/2), layer = layers['substrate_layer'])
    #topCell.add(substrate)
    
    chipCell = gdspy.Cell('chipCell')
    #Chip bounds 
    chipCell.add(gdspy.Rectangle((-chipx/2, -chipy/2), (chipx/2, chipy/2), layer = layers['chipbounds_layer'])) 
    # #dicing marks 
    chipCell.add(gdspy.Rectangle((-dicex/2, -chipy/2+dice_offset), (dicex/2, -chipy/2+dice_offset+dicey), layer = layers['dice_layer']))
    chipCell.add(gdspy.Rectangle((-dicex/2, chipy/2-dice_offset), (dicex/2, chipy/2-dice_offset-dicey), layer = layers['dice_layer']))
    
    # chipCell.add(gdspy.CellReference(current_res, (0,0))) 
    # chipCell.add(Snail_array_cell)
    

###########################Put together the top cell ###############################################

waferCell = gdspy.Cell('waferCell')
make_wafer(wafer_size=wafer_size,flat_length_primary=22*mm,flat_length_secondary=0,layer=100,cell=waferCell)

topCell = gdspy.Cell('topCell')

topCell.add(gdspy.CellReference(waferCell, (0,0)))


offset_x=chipx/2 if array_nx%2==0 else 0
offset_y=chipy/2 if array_ny%2==0 else 0

# print(array_nx,array_ny)
# print(chipx,chipy)

TL_center_x=int(-chipx*np.floor(array_nx/2)+offset_x)
TL_center_y=int(chipy*np.floor(array_ny/2)-offset_y)

text_off_x=0
text_off_y=-1*mm


if standard_layout:
        skip_layers=[(0,0),(0,1),(0,2),(0,3),(0,10),(0,11),(0,12),(0,13),
                     (1,0),(1,13),(4,0),(4,13),(5,0),(5,1),(5,2),(5,3),(5,10),(5,11),(5,12),(5,13)] #elemnts in the array to skip because they hang off the chip
else:
    skip_layers=[]
    
i=0
for x in range(array_nx):
    for y in range (array_ny):
        loc_x=TL_center_x+x*chipx
        loc_y=TL_center_y-y*chipy
        if (x,y) not in skip_layers:
            # print(i)
            for nr in range(num_res):
                if i<np.sum(res_layout[0:(nr+1)]):
                    # print(i,nr,np.sum(res_layout[0:(nr+1)]))
                    current_res=res_cells[nr]
                    # print(current_res)
                    break
            topCell.add(gdspy.CellReference(current_res, (loc_x,loc_y)))
            topCell.add(gdspy.CellReference(Snail_array_cell,(loc_x,loc_y)))
            topCell.add(gdspy.CellReference(chipCell, (loc_x,loc_y)))

            topCell.add(gdspy.Text('{} {} {}'.format(y,x, current_res.name), 2.25,(loc_x+text_off_x-0.5*mm,loc_y+text_off_y),layer=layers['dice_layer']))
            i=i+1          

    
####################### Draw test structure pads ############################

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

if Rotate: 
    oldpadCell.copy(name='SNAIL_pads_rot',deep_copy=False,translation=None,rotation=np.pi/2)

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


########################## write the gds ##############################
gdspy.write_gds(design_name + '.gds', unit=1.0e-6, precision=1.0e-9)
#gdspy.gds_print(design_name + '.gds', unit=1.0e-6, precision=1.0e-9) #old gdspy version

gdspy.current_library = gdspy.GdsLibrary() #resets libarary so you don't need to reset console 


########################## Output doses ###################################
if OutputDoses:
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
    
if clearGDSLibrary:
     gdspy.current_library = gdspy.GdsLibrary() #allows you to re-run same file in the same console and re-make the file 



