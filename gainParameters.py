#############################
# This code uses a tsv file to get a gain spected parameter from the different MCMs and chips inside them.
# Using this information creates a ROOT file to use in combination with histograms to do a calibration process
#
# To use this code you need first an output directory to put the output ROOT files
# You will need a tsv file with expected gains, calculated in previous analysis
# finally you will need to choose a spcecific RunID to get the spected value of gain
#
# The comand line example is somethin like this:
# python3.10 gainParameters.py /path/to/wite_out/root_files_root/ /path/from/file_to_get/gain_info/monitoring_DB_mod.tsv 221
#############################

from ROOT import TFile, TTree
from array import array
import funciontsv  #libreria de carlos moret
import sys


def main(path='/path/to/wite_out/root_files_root/',tsv_path='/path/from/get/gain_info/', RunID=221):
    MCMID=array('f', [-9999])
    Data_OK=array('f', [-9999])
    Gain=array('f', [-9999])
    if path != '/path/to/wite_out/root_files_root/':
        try:
            if isinstance(path, str):
                output_location = path #'/home/oem/datosFits/DarkBeats/data/rootFiles_200kADU/'
            
            
                MCMID_number=["49","19","15","20","18","48","47","21","46","44","42A","45","48B"] #MCMID

                try:
                    if tsv_path != '/path/from/get/gain_info/':
                        #"/home/oem/datosFits/DarkBeats/data/monitoring_DB_mod.tsv"
                        #RunID=211
                        j=0
                        for id in MCMID_number:
                            file_root_name = 'MCM{}_J_gain_parameters.root'.format(id)
                            file = TFile.Open(output_location  + file_root_name, "RECREATE")
                            tree = TTree('tree', 'tree')
                            ##
                            # Aqui va el codigo de Carlos o Stephane
                            ##
                            MCMID_Data, DataOK_data,gainData = funciontsv.listsFromTSV(tsv_path,RunID,id) # "/home/oem/Software/SSocial/Carlosmoret/practicas_icn/monitoring_DB.tsv"
                                
                            ### dummy data
                            # MCMID_Data=['49','49','49','49','49','49','49','49','49','49','49','49','49','49','49','49']
                            # gainData=[274.76, 277.87, 265.73, 278.37, 272.92, 268.7, 290.65, 276.6, 280.04, 275.31, -1, 286.98, 279.15, 274.16, 273.33, 283.78]
                            # DataOK_data=[1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1]
                            ###

                            tree.Branch("MCMID",MCMID,"MCMID/F")

                            tree.Branch("DataOK",Data_OK,"DataOK/F")

                            tree.Branch("Gain",Gain,"Gain/F")

                            



                            for i in range(0, len(MCMID_Data),1):
                                if isinstance(MCMID_Data[i], str): 
                                    MCMID[0] = float(MCMID_Data[i][0:2])
                                else: 
                                    MCMID[0] = MCMID_Data[i]
                                Data_OK[0] = DataOK_data[i]
                                Gain[0] = gainData[i]
                                

                                tree.Fill()

                            tree.Write()
                            file.Close()
                except:
                    print("incorrect path or tsv file path")
        except:
            print("incorrect path or RUNID")

    else:
        print("thank you")

if __name__ == "__main__":
    try:
        #print(str(len(sys.argv)))
        RunID=int(sys.argv.pop())
        tsv_path=sys.argv.pop()
        path=sys.argv.pop()
        print("output directory="+path)
        print("tsv file="+tsv_path)
        print("selected RunID reference="+str(RunID))
    except:
        print("You need to put:\n /path/to/wite_out/root_files_root/ /path/from/get/gain_info/ RunID\n" \
        "RunID must be a integer number")
    
    main(path, tsv_path, RunID)
