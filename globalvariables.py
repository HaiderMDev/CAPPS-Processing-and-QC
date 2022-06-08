from ssl import ALERT_DESCRIPTION_HANDSHAKE_FAILURE

def init_vars():

    #Check to see if packages are installed, if not it will install them. This is a list!
    global pip_packages       
    pip_packages =             "pandas", "datetime", \
                               "tkinter", \
                               "turtle", \
                               "macs2", \
                               "customtkinter", \
                               "argparse", \

    #Dictionary of experimental settings. It is composed of keys (On the left) and values (oN THE ). This is a dictionary
    global EXPERIMENT_SETTINGS
    EXPERIMENT_SETTINGS = {                             
        "Experiment name"       : "",                                         # name of experiment
        "Experiment directory"  : "",                                         # experiment directory
        "Map Quality"           : "",                                         # Map quality
        "Cores"                 : "",                                         # Number of CPU threads
        "picard directory"      : "",                                         # Directory for picard.jar
        "Genome Location"       : "",                                         # Directory for the genome
        "Bioawk Directory"      : ""                                          # Directory where bioawk is installed
    }

    #These lists store the various file names that are generated during filtering. This is a list!
    global bamfileLocation, MapQualityFileLocation, MitochondrialReadsFiltered, Sorted_chrM, peak_folder, DupsRemoved, genomeLocation, broadpeak_files, ShiftReads_files, picard_directory
    bamfileLocation             = " ",                                       
    MapQualityFileLocation      = " ",
    MitochondrialReadsFiltered  = " ",
    Sorted_chrM                 = " ",
    peak_folder                 = " ",
    DupsRemoved                 = " ",
    genomeLocation              = " ",
    broadpeak_files             = " ",
    ShiftReads_files            = " ",
    
    
    #These lists contain file names that are genreated during BigWig generation.This is a list!
    global bedgraphlist, chromosome_sizeDir, clippedBedGraphFiles, sorted_clipped_bedgraph, BigWigFileList
    bedgraphlist                = " ",
    chromosome_sizeDir          = " ",
    clippedBedGraphFiles        = " ",
    sorted_clipped_bedgraph     = " ",
    BigWigFileList              = " ",

    
    #These are the various file names that are given to files generated during filtering and bigwig file conversion. This is a dictionary!
    global File_Flags 
    File_Flags = {
        "Experimental settings":             "_Exp_Settings_.txt",                      
        "mapquality":                        "-MapQuality",                                        
        "fragment size distribution":        "_fragmentsizedistribution",
        "Mitochondrial Reads Filtered":      "-chrM",
        "sortchrM_Reads":                    "-Sorted",
        "Peak_Calls":                        "-PEAKS",
        "Index_BamFile":                     "-index.bai",
        "Remove Duplicates":                 "-Dups_Removed",
        "Bed File":                          "-BroadPeak",
        "Shifted":                           "-Shift",
        "At_chr sizes file":                 "At_chr.sizes",
        "Clipped Bed Graph File":            "-Clipped",
        "Clipped and sorted Bedgraph file":  "-Clipped_sorted",
        "BigWig File":                       "-BigWig" 
    }                                     


    # Names of directories that are generated. This is a dictionary!
    global Directory_Flags
    Directory_Flags = {
        "MapQuality"                 : "MapQuality",
        "Fragment Size Distribution" : "Fragment_Size_Dist_Graphs",
    }
    

    #This contains information for the buttons that are used for analysis. This is a list!
    global progress_bar, current_step, startbutton, PreviousButton                       
    current_step = 0                                                        # stores the current step in analysis
