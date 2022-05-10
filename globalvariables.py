def init_vars():
    #Dictionary of experimental settings. It is composed of keys (On the left) and values (oN THE )
    global EXPERIMENT_SETTINGS
    EXPERIMENT_SETTINGS = {                             
        "Experiment name": "",                                              # name of experiment
        "Experiment directory": "",                                         # experiment directory
        "Map Quality": "",                                                  # Map quality
        "Cores"      : ""
    }

    global bamfileLocation, MapQualityFileLocation, MitochondrialReadsFiltered, Sorted_chrM, peak_folder, DupsRemoved, genomeLocation, broadpeak_files, ShiftReads_files
    
    bamfileLocation             = " ",                                       #Path to the bam file locatio
    MapQualityFileLocation      = " ",
    MitochondrialReadsFiltered  = " ",
    Sorted_chrM                 = " ",
    peak_folder                 = " ",
    DupsRemoved                 = " ",
    genomeLocation              = " ",
    broadpeak_files             = " ",
    ShiftReads_files            = " ",


    global File_Flags 
    File_Flags = {
        "Experimental settings":        "_Exp_Settings_.txt",                      #log file for flags
        "mapquality":                   "-MapQuality",                                        
        "fragment size distribution":   "_fragmentsizedistribution",
        "Mitochondrial Reads Filtered": "-chrM",
        "sortchrM_Reads":               "-Sorted",
        "Peak_Calls":                   "-PEAKS",
        "Index_BamFile":                "-index.bai",
        "Remove Duplicates":            "-Dups_Removed",
        "Bed File":                     "-BroadPeak",
        "Shifted":                      "-Shift"                     
    }                                     

    global DiffBind
    DiffBind = {
        " "
    }

    global Directory_Flags
    Directory_Flags = {
        "MapQuality" : "MapQuality",
        "Fragment Size Distribution" : "Fragment_Size_Dist_Graphs",
    }
   
    global progress_bar, total_steps, current_step, startbutton, PreviousButton                       
    current_step = 0                                                        # stores the current step in analysis