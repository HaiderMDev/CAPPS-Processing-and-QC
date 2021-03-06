# version of the application
#Go to this webpage for color selection: https://encycolorpedia.com/3d6fdb or https://i.stack.imgur.com/nCk6u.jpg
app_version = "2.1"


                         ######### ----------- IMPORTING  PACKAGES ----------- #########


#Import modules and laod global vaiables python file
import os
import subprocess
import sys
import _thread
import tkinter as tk
import customtkinter as ctk
from tkinter import font as tkfont
from tkinter import *
from tkinter import filedialog
from tkinter import messagebox
from tkinter.ttk import Button, Progressbar, Style
import subprocess
from turtle import bgcolor, end_fill, write
import pandas as pd
from datetime import datetime
import argparse
import globalvariables as global_vars

#initialize golabl variabbles and load list of packages
global_vars.init_vars().__init__           



                            ######### ----------- IMPORT COMPLETE ----------- #########






######### ----------- SETUP DISPLAY FOR DOCKER CONTAINER ----------- #########

#function to detect display
if os.environ.get('DISPLAY','') == '':
    print('no display found. Using :0.0')
    os.environ.__setitem__('DISPLAY', '127.0.0.1:0.0')







                                ######### ----------- CAPPS v2.1 CODE ----------- #########


# function to creates a log file and writes messages to it. The Log.txt file is located in the Experiment directory that the user identified:
def write_to_log(message):
    #Identify the file path or change the directory to it: 
    file_path = os.path.join(global_vars.EXPERIMENT_SETTINGS['Experiment directory'] + "/")
    os.chdir(file_path)

    #Print all the messages that are being generated in the "status label" section of the analysis screen
    global_vars.run_log.config(state=NORMAL)
    global_vars.run_log.insert(END,"\n\n %s" % message)
    global_vars.run_log.config(state=DISABLED)

    #Function to write messages to the file called "Log file"
    file_name="Log.txt"
    try:
        file=open(file_name,'a+')
        file.write("\n[%s] %s" %(str(datetime.now()),message))
        file.close()
    except IOError:
        print("\n Warning: Cannot write to the log file! \n Please make sure that file can be edited")


# function to check if the file located at file_path is accessible and can be read
def check_file_access(file_path):
    # returns True only if file can be accessed
    try:
        file=open(file_path, 'r')   #r checks to ensure that the file can be read
        file.close()
    except IOError:
        write_to_log("Could not get access to this file!")
        messagebox.showerror(title="Error",
        message="Your file cannot be opened. Please check permisson for this file")
        return False
    return True

#parsing the color schemes using argparse:
ap =argparse.ArgumentParser()
ap.add_argument("-cbg", "--color_of_background", default= "#d8dcd6", help="color of background")
ap.add_argument("-cbx", "--color_of_boxes", default= "#d8dcd6", help="color of background")
ap.add_argument("-chead", "--color_of_headings", default= "#000000", help="color of background")
args = vars(ap.parse_args())



class ProcessingTool(tk.Tk):

    def __init__(self, *args, **kwargs):
        
        tk.Tk.__init__(self, *args, **kwargs)

        #Main container. This is for the section that everything will be displayed.
        #This sets the color of the frame, and how the items are attached to the grid. 
        ctk.set_appearance_mode("light")  # Modes: system (default), light, dark
        ctk.set_default_color_theme("blue")  # Themes: blue (default), dark-blue, green               
        
        #set up container:
        container = ctk.CTkFrame(self, bg= "#d8dcd6") #8C8C8C: grey; 000000: BLACK; 3d6fdb: BLUE; FFF: white
        container.pack(side="top", fill="both", expand=TRUE)
        container.grid_rowconfigure(0, weight=1)
        container.grid_columnconfigure(0, weight=1)

        # widget styles. This is settings for the widgets that will be displayed in the tkinter windows!
        self.style = Style()                
        self.style.theme_use('clam')     #For windows you can use: vista, clam, alt, classic, default, winnative,xpnative 
        self.frame_relief = RIDGE          #Options are RAISED, FLAT, SUNKEN, GROOVE, RIDGE
        self.frame_borderwidth = 2          #Determines how raised should the widgets be 
       
        # This creates the title at the top of the window
        self.winfo_toplevel().title("CAPPS")

        #set fonts
        self.title_font    = tkfont.Font(family='Arial', size=20, weight="bold")    #Font for the main title, make the weight normal
        self.Mainpage_font  = tkfont.Font(family='Arial',size=16, weight="normal")    #Font for the main title, make the weight normal
        self.label_font    = tkfont.Font(family='Arial', size=14, weight="bold")    #Font for labels
        self.maintext_font = tkfont.Font(family='Arial',size=18, weight="bold")    #Font for MainText
        self.picking_text  = tkfont.Font(family="Arial", size=15, weight="normal")  #font for picking frame                 

        #Define page frames. This will depend on the number of pages you are planning on coding for as classes!
        self.frames = {}
       
        for frame in (MainPage, information,bamfileSelect,SummaryPage,Filtering_PeakCalling):
            Title=frame.__name__
            frame = frame(parent=container, controller=self)
            self.frames[Title] = frame 
            frame.grid(row=0, column=0, sticky="nsew")

        # load main page frame at startup
        self.show_frame("MainPage")
    
    #Brings selected frame to the top
    def show_frame(self, page_name):
        frame = self.frames[page_name]
        frame.tkraise()
    
    # function to select experiment directory
    def browse_folder(self, frame_name, selection_name):
        target_obj = getattr(self.frames[frame_name], selection_name) 
        #The getattr() method returns the value of the named attribute of an object. 
        #provides classes and factory functions for creating file/directory selection windows
        #initialdir - the directory that the dialog starts in
        folder = filedialog.askdirectory(initialdir=os.getcwd(),title='Select Your Folder') 
        target_obj.delete(0, tk.END)  #Ensures that if the user picked a wrong directory, when they reselect it deletes what they previously picked. 
        target_obj.insert(tk.END, folder) 
    
    #define function to select files:
    #I have tested it, and I find that the frame_name needs to be present. It identifies which frame this function is being called in
    #Selection_name refers to the entry or list box in which the directory to *.bam or *.fastq files will be displayed. 
    def browse_file2(self,frame_name,multi_files,selection_name,file_type):
        target_obj=getattr(self.frames[frame_name],selection_name)
        if file_type == "fastq":
            file_options=(("Fastq File","*.fastq"),("All Files","*"))
        elif file_type == "bam":
            file_options=(("Bam File","*.bam"),("All Files","*"))
        elif file_type=="txt":
            file_options=(("Text File","*.txt"),("All Files","*"))
        elif file_type == "fa":
            file_options=(("Fa File","*.fa"),("fna File","*.fna"))
        elif file_type == "jar":
            file_options = (("Java file","*.jar"),("All Files","*"))
        elif file_type == "bioawk":
            file_options = (("bioawk file","*. "),("All Files","*.*"))
        else:
            messagebox.showerror(title="Error",
                message="You must select a file for any downstream work")
        if multi_files == "TRUE":
            files=filedialog.askopenfilenames(initialdir=os.getcwd(),title="Select your Files",filetypes=(file_options))
            for item in files:
                target_obj.insert(END,item)
        else:
            file =filedialog.askopenfilenames(initialdir=os.getcwd(),title="Select your Files",filetypes=(file_options))
            target_obj.delete(0,END)
            target_obj.insert(END,file)

    # function to delete selected file from file list
    def delete_selected(self, frame_name, target_name):
        target_obj = getattr(self.frames[frame_name], target_name)
        # curselection is a predefined function that fetches the value(s) of a selected item or items.
        selected = target_obj.curselection() 
        for i in selected[::-1]:
            target_obj.delete(i)


 
class MainPage(tk.Frame):
    def __init__(self, parent, controller):
   
        tk.Frame.__init__(self, parent) #Initialize the new tkinter frame under the parent

        self.controller = controller    

        #This is controls for setting up the frame: 
        frame = ctk.CTkFrame(self, relief=FLAT, borderwidth=controller.frame_borderwidth, bg= args["color_of_background"]) #grey, FFF is white
        frame.pack(fill=BOTH, padx=5, pady=5, expand=TRUE) 
        #If expand = FALSE, the frame won't take up the whole page. 
        #pady and padx: How many pixels to pad widget, horizontally and vertically, outside v's borders.

        #info for the main TOP WIDGET which contains the name of the tool in mainpage 
        frame1 = tk.Frame(frame, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth, bg = args["color_of_headings"]) #black
        frame1.pack(fill=X, padx=5, pady=5, anchor=N, expand=FALSE)
        label1_info = tk.Label(frame1, text = " ChIP-seq and ATAC-seq Processing & Peak-calling Software (CAPPS) ",
                                wraplength="800", justify="left", font = controller.title_font, bg = args["color_of_headings"], fg= "white")
        label1_info.pack(padx=2, pady=2)

        #Settings for the second widget frame which contains information about the tool in mainpage
        frame2 = ctk.CTkFrame(frame, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth, bg_color= args["color_of_background"]) #blue is 3d6fdb
        frame2.pack(fill=X, side=TOP, padx=5, pady=5, anchor=CENTER, expand=TRUE)
        label2_info = tk.Label(frame2, text=
                                     'INTRODUCTION'
                                    "\n CAPPS is designed for processing and peak calling for ATAC-seq and ChIP-seq raw BAM files."                      
                                    "\n\n OUTPUT"
                                    "\n This tool processes 'BAM' files generated from the Bowtie2 aligner and relies on the MACS2, Bamtools, bedtools, Samtools and Picards for processing."
                                    " The BigWig conversion uses bioawk, bedclip, bedtools and BedgraphtoBigWig packages. These packages can be installed using conda. Data visualization is available separately in R-Shinny as CAAT version 2.0."
                                    "\n \n CAPPS generates indexed, sorted and de-duplicated BAM file(s), as well as, fragment size distribution graph, "
                                    " a peak folder, a BigWig file and bed file for each raw BAM file. There is no limit on the number of BAM files that can be processed simultaneously."
                                    "\n\n IMPORTANT"
                                    "\n CAPPS will generate an error if there are spaces in directories or names of the BAM files. You can fill spaces with underscores."
                                    " For example: './User Directory' will generate an error, but './User_Directory' or './UserDirectory' is acceptable.",
                                wraplength="650", justify="center", font =controller.Mainpage_font, bg= args["color_of_background"], fg="black")
        #Setting for the lable which is contained in the fr_2 frame. This one contains the version of the tool:
        label2_info.pack(side=TOP, padx=5, pady=5)
        label_extra = tk.Label(frame2, text ="\n CAPPS version " + app_version,
                            wraplength="500", justify="center", font =controller.maintext_font, bg = args["color_of_background"], fg="black")
        label_extra.pack(side=TOP, padx=5, pady=5)

        #Button to move to the next frame, which is where the user selects experimental settings. 
        #if you pick "self" as the frame, then it will put the button at the bottom of the tk window. 
        #You can alternatively use "frame" as the option, then it will place the button in the frame you defined. 
        btnNewExp = ctk.CTkButton(self, text="Experiment Information", 
                                  command=lambda: controller.show_frame("information"))
        btnNewExp.pack(side=RIGHT, padx=5, pady=1)
        self.pack(fill=BOTH, expand=TRUE)

        #Button to exit the application
        btnBack = ctk.CTkButton(self, text="Exit", command=lambda: app.destroy())
        btnBack.pack(side=LEFT, padx=5, pady=5)



class information(tk.Frame):
    def __init__(self, parent, controller):
        
        #Load the frame:
        tk.Frame.__init__(self,parent)
        self.controller = controller
        
        #Set up the frame for the body:
        frame_information=tk.Frame(self, relief=FLAT, borderwidth=controller.frame_borderwidth,bg= args["color_of_background"])
        frame_information.pack(fill=BOTH,padx=5,pady=5,expand=TRUE)

        #Heading for the Frame 
        frame1_header=tk.Frame(frame_information,relief=GROOVE,borderwidth=controller.frame_borderwidth,bg= args["color_of_headings"])
        frame1_header.pack(fill=BOTH,pady=5,padx=5,anchor=N)
        frame1_header_label=tk.Label(frame1_header,text="Experiment Setup",bg= args["color_of_headings"], font=controller.title_font, fg="white") #fg is foreground color, and bg is background color
        frame1_header_label.pack(fill=X,pady=5,padx=5,anchor=N,expand=TRUE) #If you set "fill = BOTH", then the label will take the whole page

        #Select the Experiment Name, directory and Phred score:
        frame1=tk.Frame(frame_information, relief=controller.frame_relief,borderwidth=controller.frame_borderwidth,bg= args["color_of_boxes"])
        frame1.pack(fill=X, pady=5,padx=5,expand=TRUE,anchor=N)
        label_frame=tk.Label(frame1,text="Main Experiment settings",width="50", font=controller.maintext_font,bg= args["color_of_boxes"],fg="black")
        label_frame.pack(side=TOP,anchor=N, pady=5,padx=5,expand=TRUE)

        label_frame=tk.Frame(frame1, relief=FLAT,borderwidth=controller.frame_borderwidth,bg= args["color_of_boxes"])
        label_frame.pack(fill=X, pady=5,padx=5,expand=TRUE,anchor=N)
        label_exp_name=tk.Label(label_frame,text="Define Experiment Name", width="40", font=controller.picking_text,bg= args["color_of_boxes"],fg="black")
        label_exp_name.pack(side=LEFT,anchor=E,pady=5,padx=5,expand=FALSE)
        self.txt_exp_name=Entry(label_frame, width="20")
        self.txt_exp_name.pack(fill=X,padx=5,pady=5,expand=TRUE,side=RIGHT)
        self.txt_exp_name.config(state=NORMAL)

        outputdir = tk.Frame(frame1,relief=FLAT, borderwidth=controller.frame_borderwidth,bg= args["color_of_boxes"])
        outputdir.pack(fill=X, pady=5,padx=5,anchor=N,expand=TRUE)
        label_exp_folder=tk.Label(outputdir,text="Experimental Directory",width="40", font=controller.picking_text,bg= args["color_of_boxes"],fg="black")
        label_exp_folder.pack(side=LEFT,anchor=N, pady=5,padx=5,expand=FALSE)
        self.txt_expdir=Entry(outputdir,width="20")
        self.txt_expdir.pack(fill=X,pady=5,padx=5,expand=FALSE,side=LEFT)
        btn_browse = ctk.CTkButton(outputdir,text="BROWSE",command=lambda:self.controller.browse_folder(
            frame_name="information",selection_name="txt_expdir"
        )) 
        btn_browse.pack(fill=X,side=RIGHT,padx=5,pady=5)
        self.txt_expdir.config(state=NORMAL)

        MapQuality=tk.Frame(frame1, relief=FLAT,borderwidth=controller.frame_borderwidth,bg = args["color_of_boxes"])
        MapQuality.pack(fill=X, pady=5,padx=5,expand=TRUE,anchor=N)
        label_map_quality=tk.Label(MapQuality,text="Choose Map Quality Score (Between 20 - 35)", width="40", font=controller.picking_text,bg= args["color_of_boxes"],fg="black")
        label_map_quality.pack(side=LEFT,anchor=N,pady=5,padx=5)
        self.txt_mapquality=Entry(MapQuality, width="20")
        self.txt_mapquality.pack(fill=X,padx=5,pady=5,expand=TRUE)


        #Selecting directory for number of cores, genome, picards and bioawk:
        frame5 = tk.Frame(frame_information,relief=controller.frame_relief, borderwidth=controller.frame_borderwidth,bg= args["color_of_boxes"])
        frame5.pack(fill=X, pady=5,padx=5,anchor=N,expand=TRUE)
        label_exp_folder3=tk.Label(frame5,text="Core Program Settings",width="50", font=controller.maintext_font,bg= args["color_of_boxes"],fg="black")
        label_exp_folder3.pack(side=TOP,anchor=N, pady=5,padx=5,expand=TRUE)

        cores = tk.Frame(frame5, relief=FLAT,borderwidth=controller.frame_borderwidth,bg= args["color_of_boxes"])
        cores.pack(fill=X,pady=5,padx=5,anchor=E,expand=TRUE)
        label_cores=tk.Label(cores,text="Number of CPU threads (Must be > 0)", width="40", font=controller.picking_text,bg= args["color_of_boxes"],fg="black")
        label_cores.pack(side=LEFT,anchor=W, pady=5,padx=5,expand=FALSE)
        self.txt_cores=Entry(cores, width="40")
        self.txt_cores.pack(fill=X,pady=5,padx=5,expand=FALSE,side=LEFT)

        mouse_genome_select = tk.Frame(frame5,relief=FLAT, borderwidth=controller.frame_borderwidth,bg= args["color_of_boxes"])
        mouse_genome_select.pack(fill=X, pady=5,padx=5,anchor=E,expand=TRUE)
        label_exp_folder2=tk.Label(mouse_genome_select,text="Select Genome (eg. mm10.fa)",width="40", font=controller.picking_text,bg= args["color_of_boxes"],fg="black")
        label_exp_folder2.pack(side=LEFT,anchor=W, pady=5,padx=5,expand=FALSE)
        self.txt_expdir2=Entry(mouse_genome_select,width="20")
        self.txt_expdir2.pack(fill=X,pady=5,padx=5,expand=FALSE,side=LEFT)

        frame_picard = tk.Frame(frame5,relief=FLAT, borderwidth=controller.frame_borderwidth,bg= args["color_of_boxes"])
        frame_picard.pack(fill=X, pady=5,padx=5,anchor=E,expand=TRUE)
        label_exp_folder3=tk.Label(frame_picard,text="Location for \"picard.jar\"",width="40", font=controller.picking_text,bg= args["color_of_boxes"],fg="black")
        label_exp_folder3.pack(side=LEFT,anchor=W, pady=5,padx=5,expand=FALSE)

        frame_bioawk = tk.Frame(frame5,relief=FLAT, borderwidth=controller.frame_borderwidth,bg= args["color_of_boxes"])
        frame_bioawk.pack(fill=X, pady=5,padx=5,anchor=W,expand=TRUE)
        label_exp_folder4=tk.Label(frame_bioawk,text="Location for bioawk (Unix Executable File)",width="40", font=controller.picking_text,bg= args["color_of_boxes"],fg="black")
        label_exp_folder4.pack(side=LEFT,anchor=N, pady=5,padx=5,expand=FALSE)

        #picard entry
        self.txt_expdir3=Entry(frame_picard,width="20")
        self.txt_expdir3.pack(fill=X,pady=5,padx=5,expand=FALSE,side=LEFT)
         #bioawk entry
        self.txt_expdir4=Entry(frame_bioawk,width="20")
        self.txt_expdir4.pack(fill=X,pady=5,padx=5,expand=FALSE,side=LEFT)

        #Button for selecting mouse diretory:
        btn_browse2 = ctk.CTkButton(mouse_genome_select,text="BROWSE",command=lambda:self.controller.browse_file2(
                                                                            multi_files="TRUE",
                                                                            frame_name="information",
                                                                            selection_name="txt_expdir2",
                                                                            file_type="fa" 
        )) 
        #Button for selecting picard diretory:
        btn_browse3 = ctk.CTkButton(frame_picard,text="BROWSE",command=lambda:self.controller.browse_file2(
                                                                            multi_files="FALSE",
                                                                            frame_name="information",
                                                                            selection_name="txt_expdir3",
                                                                            file_type="jar" 
        ))
        #Button for selecting bioawk diretory:
        btn_browse4 = ctk.CTkButton(frame_bioawk,text="BROWSE",command=lambda:self.controller.browse_file2(
                                                                            multi_files="FALSE",
                                                                            frame_name="information",
                                                                            selection_name="txt_expdir4",
                                                                            file_type= "bioawk" 
        ))

        #Using browse_folder definition under main class ATACseqAlign. "frame_name" has to be the current class "information". 
        #"target_name" has to be the "txt_expdir" or it won't copy the directory that was picked. 
        btn_browse2.pack(fill=X,side=RIGHT,padx=5,pady=5)
        btn_browse3.pack(fill=X,side=RIGHT,padx=5,pady=5)
        btn_browse4.pack(fill=X,side=RIGHT,padx=5,pady=5)

        #Navigation Buttons
        btnexit=ctk.CTkButton(self, text="Exit",command=lambda: app.destroy())
        btnexit.pack(side=LEFT,padx=5,pady=5)
        #Button for next page
        btnnextpage=ctk.CTkButton(self, text="Next Page", command=lambda:self.store_variables()) #"Store Variables" is defined below!
        #btnnextpage=Button(self, text="Next Page", command=lambda:controller.show_frame("bamfileSelect"))
        btnnextpage.pack(side=RIGHT,pady=5,padx=5)
        #Button to return to the last page
        btnprevious=ctk.CTkButton(self, text="Main Page", command=lambda:controller.show_frame("MainPage"))
        btnprevious.pack(side=RIGHT,pady=5,padx=5)
    

    #Define the functions to store the values entered into global varaibles. These variables will be called later.
    def store_variables(self):
        error_log = ""

        #Check and save user input into global variables:
        #Experiment name
        if self.controller.frames['information'].txt_exp_name.get() != "":
            global_vars.EXPERIMENT_SETTINGS['Experiment name'] = self.controller.frames['information'].txt_exp_name.get()  # experiment name
        else:
            error_log = error_log + "Experiment name cannot be blank"

        #Directory:
        if self.controller.frames['information'].txt_expdir.get() != "":
            global_vars.EXPERIMENT_SETTINGS['Experiment directory'] = self.controller.frames['information'].txt_expdir.get() #Experiment Directory
        else:
            error_log = error_log + "\n A working directory must be selected"

        #check to see if a genome has been stored
        if self.controller.frames['information'].txt_expdir2.get()  != "":
            global_vars.genomeLocation = self.controller.frames['information'].txt_expdir2.get() 
            global_vars.EXPERIMENT_SETTINGS['Genome Location'] = self.controller.frames['information'].txt_expdir2.get() 
        else:
            error_log = error_log + "\n Please select a genome for filtering!"
            
        #Check to make sure that MapQuality are is between 20 and 35 and saves it in global variables:
        if self.controller.frames['information'].txt_mapquality.get() != "":
            if int(self.controller.frames['information'].txt_mapquality.get()) in range(20,36): #checks if user selected score between 20 - 35
                try:
                    global_vars.EXPERIMENT_SETTINGS['Map Quality']= int(self.controller.frames['information'].txt_mapquality.get())
                except ValueError:
                    global_vars.EXPERIMENT_SETTINGS['Map Quality'] = 25
            else:
                error_log = error_log + "\n Map Quality must be a number between 20 - 35!"
            if not self.controller.frames['information'].txt_mapquality.get().isdigit():
                error_log = error_log + "\n Map Quality must be a number!"
        else:
            error_log = error_log + "\n Please select a Maq Quality score!"

        #Save the number of cores in global variables:
        if self.controller.frames['information'].txt_cores.get() != "":
            try:
                global_vars.EXPERIMENT_SETTINGS['Cores'] = int(self.controller.frames['information'].txt_cores.get())
            except ValueError:
                global_vars.EXPERIMENT_SETTINGS['Cores'] = 1
        else:
            error_log = error_log + "\n Please select a number of computing cores"

        #check to make sure that Cores is a number and store its value
        if self.controller.frames['information'].txt_cores.get() != "": 
            if not self.controller.frames['information'].txt_cores.get().isdigit():
                error_log = error_log + "\n CPU threads selection must be identified as a number!"

        #make sure the number of cores entered does not exceed the cores available in the computer:
        if self.controller.frames['information'].txt_cores.get() != "": 
            if int(self.controller.frames['information'].txt_cores.get()) > os.cpu_count():
                cpu_count=os.cpu_count()
                error_log = error_log + "\n Your computer has {} CPU threads. Please choose a number equal to or below that number".format(cpu_count)
       
        #Save picard directory:
        if self.controller.frames['information'].txt_expdir3.get()  != "":
            global_vars.EXPERIMENT_SETTINGS['picard directory'] = self.controller.frames['information'].txt_expdir3.get()
        else:
            error_log = error_log + "\n Please select a directory for picard before continuing!"

          #Save bioawk directory:
        if self.controller.frames['information'].txt_expdir4.get()  != "":
            global_vars.EXPERIMENT_SETTINGS['Bioawk Directory'] = self.controller.frames['information'].txt_expdir4.get()
            
        #Display messages from the error log
        if error_log != "":
            messagebox.showerror(title="Error", message="%s" % error_log)
        else:
            self.controller.show_frame("bamfileSelect")



class bamfileSelect(tk.Frame):
    def __init__ (self, parent, controller):
        #Load the frame:
        tk.Frame.__init__(self,parent)
        self.controller = controller
        
        #Set up the main frame:
        frame_bamfileSelect=tk.Frame(self, relief=FLAT, borderwidth=controller.frame_borderwidth,bg= args["color_of_background"])
        frame_bamfileSelect.pack(fill=BOTH,padx=5,pady=5, expand=TRUE)


        #Heading for the Frame 
        frame1_header=tk.Frame(frame_bamfileSelect,relief=RIDGE,borderwidth=controller.frame_borderwidth,bg= args["color_of_headings"])
        frame1_header.pack(fill=X,pady=5,padx=5,anchor=N)
        frame1_header_label=tk.Label(frame1_header,text="Load BAM Files",bg= args["color_of_headings"], font=controller.title_font, fg="white") #fg is foreground color, and bg is background color
        frame1_header_label.pack(fill=X,pady=5,padx=5,anchor=N,expand=TRUE) 


        #Frame for loading bam files:
        frame4=tk.Frame(frame_bamfileSelect,relief=RIDGE,borderwidth=controller.frame_borderwidth,bg= args["color_of_boxes"])
        frame4.pack(fill=X,pady=5,padx=5,anchor=N)
        frame1_header_label=tk.Label(frame_bamfileSelect,text="\n \n \n \n A BAM file (*. bam) is the compressed \n binary version of a SAM file" \
                                                               "\n that is used to represent aligned  \n sequences up to 128 Mb."\
                                                               " \n A BAM file is generated after \n alignment using the BOWTIE."
                                                               ,bg = args["color_of_boxes"], width="30", font=controller.picking_text, fg="black")
        frame1_header_label.pack(fill=X,pady=5,padx=5,anchor=N,side=LEFT) 

        #Frame for list container and buttons
        #Create a frame which will house the container that is scrollable. 
        listcontainer = tk.Frame(frame_bamfileSelect, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth, bg= args["color_of_boxes"])
        listcontainer.pack(fill=BOTH, padx=5, pady=5, anchor=N, expand=TRUE)
        #Create a frame within the listcontainer_scroll which will contain the scroll options. For this frame use "fill = BOTH" for packing. 
        listcontainer_scroll = tk.Frame(listcontainer, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth, bg= args["color_of_boxes"])
        listcontainer_scroll.pack(fill=BOTH, padx=5, pady=5, side=RIGHT, expand=TRUE) #If you don't select "fill=BOTH" it won't fill the whole listcontainer frame
        #Button for container 
        list_btn_container = tk.Frame(listcontainer_scroll, bg= args["color_of_boxes"])
        list_btn_container.pack(fill=BOTH, padx=0, pady=0, side=RIGHT)

        #Listbox and scrolling 
        self.list_bamfiles = Listbox(listcontainer_scroll,height="22",width="300") #selectmode=EXTENDED)
        #Vertical scroll for final listbox
        yscrl_list_bamfiles= Scrollbar(listcontainer_scroll, orient=VERTICAL)
        yscrl_list_bamfiles.pack(side=RIGHT, fill=Y)
        yscrl_list_bamfiles.config(command=self.list_bamfiles.yview)
        self.list_bamfiles.config(yscrollcommand=yscrl_list_bamfiles.set)
        #Horizontal scroll for final listbox
        xscrl_list_bamfiles= Scrollbar(listcontainer_scroll, orient=HORIZONTAL)
        xscrl_list_bamfiles.pack(side=BOTTOM, fill=X)
        xscrl_list_bamfiles.config(command=self.list_bamfiles.xview)
        self.list_bamfiles.config(xscrollcommand=xscrl_list_bamfiles.set)
        #now pack the lst_final_files
        self.list_bamfiles.pack(side=RIGHT, padx=5, pady=5)

        #Buttons to add bam files
        button_add=ctk.CTkButton(list_btn_container,text="Add Files",
                            command=lambda:self.controller.browse_file2(
                                multi_files="TRUE",
                                frame_name="bamfileSelect",
                                selection_name="list_bamfiles",
                                file_type="bam" 
                            ))
        button_add.pack(side=TOP,padx=0,pady=0)

        #Button to remove the files selected
        button_remove=ctk.CTkButton(list_btn_container,text="Remove Files",command=lambda:self.controller.delete_selected(
            frame_name="bamfileSelect",
            target_name="list_bamfiles"
        ))
        button_remove.pack(side=BOTTOM,padx=0,pady=0)

        #Button to clear the selection:
        button_clear=ctk.CTkButton(list_btn_container,text="Clear Selecton",command=lambda:self.list_bamfiles.delete(0,END))
        button_clear.pack(side=BOTTOM,padx=0,pady=0)

        #Navigation Buttons
        btnexit=ctk.CTkButton(self, text="Exit",command=lambda: app.destroy())
        btnexit.pack(side=LEFT,padx=5,pady=5)
        #Button for next page
        btnnextpage=ctk.CTkButton(self, text="Summary Page", command=lambda:self.bamfile_dir())
        btnnextpage.pack(side=RIGHT,pady=5,padx=5)
        #Button to return to the last page
        btnprevious=ctk.CTkButton(self, text="Previous", command=lambda:controller.show_frame("information"))
        btnprevious.pack(side=RIGHT,pady=0,padx=0)


    #Save the location of bam file to global variables:
    def bamfile_dir(self):
        errorlog = ""
        if self.controller.frames['bamfileSelect'].list_bamfiles.get(0) != "":
            global_vars.bamfileLocation=self.controller.frames['bamfileSelect'].list_bamfiles.get(0,END)
        else:
            errorlog= errorlog + "\n\n Please select a BAM file to continue."

        if errorlog != "":
            messagebox.showerror(title="Error", message= "%s" % errorlog)
        else:
            self.controller.show_frame('SummaryPage')



class SummaryPage(tk.Frame):
    def __init__(self,parent,controller):
        tk.Frame.__init__(self,parent)
        self.controller=controller

        main=tk.Frame(self,relief=FLAT,borderwidth=controller.frame_borderwidth,bg= args["color_of_background"])
        main.pack(fill=BOTH,pady=5,padx=5,expand=TRUE)

        #Heading for the frame
        heading=tk.Frame(main,relief=controller.frame_relief,borderwidth=controller.frame_borderwidth,bg= args["color_of_headings"])
        heading.pack(fill=BOTH,anchor=N,pady=5,padx=5)
        heading_label=tk.Label(heading,text="Summary of Experimental Settings",fg="white",bg= args["color_of_headings"],font=controller.title_font)
        heading_label.pack(fill=BOTH,anchor=N,pady=5,padx=5)

        summaryframe=tk.Frame(main,relief=FLAT,borderwidth=controller.frame_borderwidth,bg= args["color_of_boxes"])
        summaryframe.pack(fill=BOTH,pady=5,padx=5,anchor=S,expand=TRUE)
        self.summary_txt=Text(summaryframe,height="25")
        self.summary_txt.pack(fill=BOTH,pady=5,padx=5,side=LEFT,expand=TRUE)
        #y-scroll
        yscrol_summary_txt=Scrollbar(self.summary_txt,orient=VERTICAL)
        yscrol_summary_txt.pack(side=RIGHT,fill=Y)
        yscrol_summary_txt.config(command=self.summary_txt.yview)
        self.summary_txt.config(yscrollcommand=yscrol_summary_txt.set)

        #Navigation Buttons
        btnexit=ctk.CTkButton(self, text="Exit",command=lambda: app.destroy())
        btnexit.pack(side=LEFT,padx=5,pady=5)
        #Button to print the Experimental Settings and BAM file directories:
        btnPrintFiles=ctk.CTkButton(main,text="Print Summary",command=lambda:self.print_summary()) #have to use self, because the definition for "print_files" is in this frame
        btnPrintFiles.pack(side=RIGHT,padx=5,pady=5)
        #Button for next page
        btnnextpage=ctk.CTkButton(self, text="Next Page", command=lambda:controller.show_frame('Filtering_PeakCalling'))
        btnnextpage.pack(side=RIGHT,pady=5,padx=5)
        #Button to return to the last page
        btnprevious=ctk.CTkButton(self, text="Previous", command=lambda:controller.show_frame("bamfileSelect"))
        btnprevious.pack(side=RIGHT,pady=0,padx=0)
 
    def print_summary(self):
        #global_vars.bamfileLocation = self.controller.frames['bamfileSelect'].list_bamfiles.get(0,END)
        
        #This is defining that the box is the frame for "summary_txt" located in the SummaryPage main frame. 
        Summary_box = self.controller.frames['SummaryPage'].summary_txt
        #These settings resets the summary box
        Summary_box.config(state=NORMAL)
        Summary_box.delete(1.0,END)

        #This will print the keys and values stored in teh EXPERIMENT SETTINGS dictionary:
        for keys, values in global_vars.EXPERIMENT_SETTINGS.items():
            text= "\n\n %s: %s\n" % (keys, values)
            self.controller.frames['SummaryPage'].summary_txt.insert(END,text)

        #This prints the number of items in the list for bam files. 
        self.number_of_items(target=global_vars.bamfileLocation)
        
        #This prints the bam files in the summary box
        self.insertFiles(target=global_vars.bamfileLocation)

        #This makes it so you can't even select the text that has been printed in the summary box
        Summary_box.config(state=DISABLED)
    
    #Function to get the number of items in the list:
    def number_of_items(self,target):
        Summary_box=self.controller.frames['SummaryPage'].summary_txt
        count=0
        for i in target:
            count +=1   
        Summary_box.insert(END, "\n \n You have selected %s file(s) for processing" % count)
        #This line of code will detect if you have selected a duplicate BAM file and give an error. 
        for i in target:
            if target.count(i) > 1:
                messagebox.showerror(title="Error", message= "You have selected the same file more than once.")
            return FALSE
            
    #Function to print the list of bam files stored in the bamfileLocation in GLOBAL VARIABLES.
    def insertFiles(self,target):
        Summary_box=self.controller.frames['SummaryPage'].summary_txt
        for i in target:
            #The index() method returns the position at the first occurrence of the specified value
            #I am add "+1" to the index() count becasue it starts from 0. After adding +1, the first file starts at 1 and not zero. 
            text_template =  "\n\n\n File %s: %s" % (target.index(i)+1,i)
            Summary_box.insert(END,text_template)



class Filtering_PeakCalling(tk.Frame):
    def __init__(self, parent, controller):
        tk.Frame.__init__(self,parent)
        self.controller = controller

        # main frame for page body content
        mainATAC = tk.Frame(self, relief=GROOVE, borderwidth=controller.frame_borderwidth, bg=args["color_of_background"])
        mainATAC.pack(fill=BOTH, padx=5, pady=5, expand=True)

        #Heading
        header = tk.Frame(mainATAC, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth, bg = args["color_of_headings"])
        header.pack(fill=X, padx=5, pady=5, anchor=N)
        lbl_header = tk.Label(header, text="BAM File Filtering and Peak Calling", font=controller.title_font, bg=args["color_of_headings"], fg="white")
        lbl_header.pack(fill=X,padx=5, pady=5, expand=True,anchor=N)

        #Frame for progress bar                            
        progressbar=tk.Frame(mainATAC,relief=FLAT,borderwidth=controller.frame_borderwidth,bg=args["color_of_boxes"])
        progressbar.pack(padx=5,pady=5,expand=TRUE,anchor=N)
  
        s = Style()
        s.theme_use('clam')
        s.configure("red.Horizontal.TProgressbar",background='red', foreground='red',thickness=100)
        global_vars.progress_bar = Progressbar(progressbar, style= "red.Horizontal.TProgressbar",orient="horizontal", mode="determinate", length="650")
        global_vars.progress_bar.pack(fill=BOTH, padx=5, pady=5, expand=True)
        global_vars.progress_bar['maximum'] = 100
        global_vars.progress_bar['value']= 0

        #Frame for status updates. The status updates will be displayed on this frame. 
        status= tk.Frame(mainATAC, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth, bg= args["color_of_boxes"])
        status.pack(fill=BOTH, padx=5, pady=5, side=TOP, anchor=N, expand=True)
        global_vars.run_log = Text(status, height="30")
        global_vars.run_log.pack(fill=BOTH, side=LEFT, padx=5, pady=5, expand=True)
        scrl_txt_summary = Scrollbar(status)
        scrl_txt_summary.pack(side=RIGHT, fill=Y)
        scrl_txt_summary.config(command=global_vars.run_log.yview)
        global_vars.run_log.config(yscrollcommand=scrl_txt_summary.set)

        #Navigation Buttons
        #Exit Button
        btnexit=ctk.CTkButton(self, text="Exit",command=lambda: app.destroy())
        btnexit.pack(side=LEFT,padx=5,pady=5)
        #Back Button:
        global_vars.PreviousButton=ctk.CTkButton(self, text="Previous Page",command=lambda:controller.show_frame("SummaryPage"))
        global_vars.PreviousButton.pack(side=LEFT,pady=0,padx=0)
        #Button to start the filtering process. Saving button this way allows us to control its configuration. 
        global_vars.startbutton=ctk.CTkButton(self,text="START",command = lambda: _thread.start_new_thread(self.analysis_steps, ()))
        global_vars.startbutton.pack(side=RIGHT, padx=5, pady=5)

    
    def update_progress_bar(self):
        global_vars.current_step += 1
        global_vars.progress_bar['value'] = (global_vars.current_step/global_vars.total_steps)*100

    #This will save experimental settings to a file in the experimental directory. 
    #I will be using the "write_to_log" definition constantly here. This was defined at the beginning of the application. The "write_to_log"
    #function just writes messages in a log file titled "logfile" that is saved in the path identified in the "write_to_log" definition
    def save_experimental_settings(self):
        #Disable the start button
        global_vars.PreviousButton.config(state=DISABLED)
        global_vars.startbutton.config(state=DISABLED)

        #Define the file path and the file name
        file_name = global_vars.EXPERIMENT_SETTINGS['Experiment name'] + global_vars.File_Flags['Experimental settings']
        file_path = os.path.join(global_vars.EXPERIMENT_SETTINGS['Experiment directory'] + "/")
        
        #Function to check if the 'Experiment directory' exists. If it doesnt exist, the program will attempt to make it. if it cant make it
        #then it will write an error in the experimental director. 
        if os.path.exists(file_path) is False:
           write_to_log("Experimental directory not found.") #Writes to a log file, which is completely separate from the file this is creating!
           try:
                os.mkdir (global_vars.EXPERIMENT_SETTINGS['Experiment directory'])
           except RuntimeError or IOError:
                write_to_log("Couldn't find Experimental directory. Attempted to create the folder, but was unsucessful") 
                return False

        #Change the working directory to the experimental directory that the user specified
        #Write the experimental settings to file. Everything that was printed on SummaryPage text box will be saved in file
        #with the title specified above as "file_name"
        os.chdir(file_path)
        write_to_log("\n\n The Experimental file %s file is saved in %s \n\n" % (file_name,file_path))
        try:
            write_to_log("Saving experiment settings to file!")
            f=open(file_name, 'w')
            f.write("ATACseq Pre Processing Tool Version: %s" % app_version)
            f.write("\nDate: %s\n\n" % str(datetime.now()))
            f.write("\n Experimental Conditions and BAM files: \n")
            f.write("\n the mouse genome is located in: %s \n\n" % global_vars.genomeLocation)
            f.write(self.controller.frames['SummaryPage'].summary_txt.get(1.0,END))
            f.close()
        except IOError or RuntimeError:
            write_to_log("Cannot write to experimental log file")
            return False
        return True

    def analysis_steps(self):

        #Define the total number of steps for progress bar
        global_vars.total_steps = 10 + len(global_vars.bamfileLocation)

        #Save experimental settings. 
        #If the experimental results cannot be saved, generate an error and save it to logfile. This is general. 
        results = self.save_experimental_settings()

        if results is False:
            write_to_log("Could not save the result!")
            messagebox.showerror(title="Error",message="Error! Could not save the result!")
            return

        #Update progress bar:
        self.update_progress_bar()

        ########################SAVING RESULTS FOR MAPQUALITY FILTERING##########################
        #After defining what happens in case of an error, start the actual process
        #First step is to filter the bam file for mitochondrial reads. Therefore, write to log that filtering mitochondrial reads
        #Then after results are complete write that it was completed.
        #This can be added to every step that will be performed so the user can view the log file and read what was completed. 
        write_to_log("\n\n Step 1: Filter by map quality score %s" % global_vars.EXPERIMENT_SETTINGS['Map Quality'])
        mapquality_results = self.MapQuality(bamfilelist = global_vars.bamfileLocation)
        if mapquality_results is False:
            write_to_log("Error. Could not filter reads by map quality")
            messagebox.showerror(title="Error!",message="Could not filter reads for map quality")
            return 

        #Update progress bar:
        self.update_progress_bar()
        
        ######################## SAVING RESULTS FOR Mitochondrial Reads ##########################
        write_to_log("\n\n STEP 2: Removing Mitochondrial Reads From the BAM Files")
        chrMfiltering_result = self.filteringchrM(bamfilelist=global_vars.MapQualityFileLocation)
        if chrMfiltering_result is False:
            write_to_log("Error. Could not filter the BAM files for mitochondrial reads")
            messagebox.showerror(title="Error!",message="Could not filter the BAM files for mitochondrial reads")
            return 

        #Update progress bar:
        self.update_progress_bar()

        ######################## SAVING RESULTS FOR Sort Mitochondrial Reads ##########################
        write_to_log("\n\n STEP 3: Sorting the BAM file")
        sortchrMreads_result = self.sortchMreads(bamfilelist=global_vars.MitochondrialReadsFiltered)
        if sortchrMreads_result is False:
            write_to_log("Error. Could not sort the BAM files")
            messagebox.showerror(title="Error!",message="Could not sort the bam file")
            return
        
        #Update progress bar:
        self.update_progress_bar()


        ######################## SAVING RESULTS FOR Removing Duplicates ##########################
        write_to_log("\n\n STEP 4: Removing Duplicate Reads From the BAM Files")
        dupsremoved_result = self.RemoveDups(bamfilelist=global_vars.Sorted_chrM)
        if dupsremoved_result  is False:
            write_to_log("Error. Could not remove duplicates from the BAM file")
            messagebox.showerror(title="Error!",message="Could not remove duplicates from the BAM file")
            return 

        #Update progress bar:
        self.update_progress_bar()
        

        ######################## SAVING RESULTS FOR index sorted Reads ##########################
        write_to_log("\n\n STEP 5: indexing sorted the BAM file with dups removed")
        sortchrMreads_result = self.indxSortReads(bamfilelist=global_vars.DupsRemoved)
        if sortchrMreads_result is False:
            write_to_log("Error. Could not index the BAM files")
            messagebox.showerror(title="Error!",message="Could not index the bam file")
            return
        

        #Update progress bar:
        self.update_progress_bar()

        ######################## SAVING RESULTS FOR Fragment Size Distribution ##########################
        write_to_log("\n\n STEP 6: Generating fragment size distribution")
        fragmentsize_result = self.FragmentSizeDist(bamfilelist=global_vars.DupsRemoved)
        if fragmentsize_result is False:
            write_to_log("Error. Could not generate fragment size distribution")
            messagebox.showerror(title="Error!",message="Could not generate fragment size distribution")
            return 
        
        #Update progress bar:
        self.update_progress_bar()


        ######################## SAVING RESULTS FOR MACS2 peak calling ##########################
        write_to_log("\n\n STEP 7: Calling peaks using MACS2!")
        PeakCalls = self.MACS2_peackcalling(bamfilelist=global_vars.DupsRemoved)
        if PeakCalls is False:
            write_to_log("Error. Could not call peaks")
            messagebox.showerror(title="Error!",message="Could not calls peak")
            return 

        #Update progress bar:
        self.update_progress_bar()


        ######################## GENERATING BED FILE FROM BROADPEAK FILE ##########################
        write_to_log("\n\n STEP 8: Generating .BED file!")
        bed_conversion = self.broadpeak_bed_conversion(bamfilelist=global_vars.broadpeak_files)
        if bed_conversion is False:
            write_to_log("Error. Could not generate bed files")
            messagebox.showerror(title="Error!",message="Could not generated bed files")
            return 

        #Update progress bar:
        self.update_progress_bar()


        ######################## REMOVING UNNECESSARY FILES ##########################
        write_to_log("\n\n STEP 9: Removing Unnecessary Files Prior to BigWig Conversion")
        RemoveFiles = self.RemoveFiles(bamfilelist = global_vars.bamfileLocation)
        if RemoveFiles is False:
            write_to_log("Error. Could not remove unnecessary filesfiles")
            messagebox.showerror(title="Error!",message="Could not remvoe unnecessary files")
            return 

         #Update progress bar:
        self.update_progress_bar()

        ######################## Generating At_chr.sizes file ##########################
        write_to_log("\n\n STEP 10: Generating AT_Chr file")
        At_chr = self.At_chr()
        if At_chr is False:
            write_to_log("Error. Could not proceed. Please ensure that bioawk is installed and the executable file is selected. You can find instructions on my GitHub Page!")
            messagebox.showerror(title="Error!",message="Could not proceed")
            return 

        #Update progress bar:
        self.update_progress_bar()

        ######################## Clipping the bedgraph file ##########################
        write_to_log("\n\n STEP 11: Clipping the bedgraph file(s)")
        bedgraphtobigwig = self.bedgraphtobigwig(bedgraphlist = global_vars.bedgraphlist)
        if bedgraphtobigwig is False:
            write_to_log("Error. Could not generate sort the bedgraph. Please ensure you have bedClip installed. This can be done using Anaconda \
                        as follows: Open terminal and type \'conda install -c bioconda/label/cf201901 ucsc-bedclip\'")
            messagebox.showerror(title="Error!",message="Could not proceed with bigwig file conversion")
            return 

        #Update progress bar:
        self.update_progress_bar()

        ######################## Sorting the clipped bedgraph file #########################
        write_to_log("\n\n STEP 12: Sorting the clipped bedgraph files before BigWig conversion")
        sortclippedbedgraph = self.sortclippedbedgraph(bedgraphlist = global_vars.clippedBedGraphFiles)
        if sortclippedbedgraph is False:
            write_to_log("Error. Could not sort the clipped bedgraph file. This is unusual. This uses  Please contact me (hhassan4242@gmail.com) for more instructions.")
            messagebox.showerror(title="Error!",message="Could not proceed with bigwig file conversion")
            return 

        #Update progress bar:
        self.update_progress_bar()

        ######################## Sorting the clipped bedgraph file #########################
        write_to_log("\n\n STEP 13: Now Performing the final BigWig File Conversion")
        BigWigConversion = self.BigWigConversion(bedgraphlist = global_vars.sorted_clipped_bedgraph)
        if BigWigConversion is False:
            write_to_log("Error. Could not make the BigWig File. Ensure you have bedGraphToBigWig installed. This can be done using Anaconda \
                         as follows: Open terminal and type \' conda install -c bioconda/label/cf201901 ucsc-bedgraphtobigwig \'")
            messagebox.showerror(title="Error!",message="Could not proceed with bigwig file conversion")
            return 

        #Update progress bar:
        self.update_progress_bar()

        ######################## Removing Unnecessary Files #########################
        write_to_log("\n\n STEP 14: Removing All Unnecessary Files Generated Using BigWig Conversion! This is the Last Step!")
        finalcleanup = self.finalcleanup(bedgraphlist = global_vars.bedgraphlist)
        if finalcleanup is False:
            write_to_log("Error. Could not remove the files. This is strange. ")
            messagebox.showerror(title="Error!",message="Could not remove unnecessary files")
            return 

        #Update progress bar:
        self.update_progress_bar()
        
        #Update progress bar and write to log:
        self.update_progress_bar()
        write_to_log("All Done!!")
        
        #Now enable the start and next buttons 
        global_vars.startbutton.config(state= NORMAL)
        global_vars.PreviousButton.config(state=NORMAL)



    def MapQuality(self,bamfilelist):
        
        #Change the working direcotory to experimental directory picked by the user
        outputdir= os.path.join(global_vars.EXPERIMENT_SETTINGS['Experiment directory'] + "/")

        #change working directory
        os.chdir(outputdir)
      
        #Save the list of output file paths in "Outputfilesdir"
        Outputfilesdir = []

        #Pipe the individual bamfile for mapquality filtering.
        for i in bamfilelist:

            outputfilename= "File" + str(bamfilelist.index(i)+1) + global_vars.File_Flags['mapquality'] + ".bam"
            readfile = i

            #Saving the resulting files in a list called"outputfilesdir". 
            Output_files_dir = outputdir + outputfilename
            Outputfilesdir.append(Output_files_dir)

            #write to log the file name:
            write_to_log("Filtering for map quality for %s" % readfile)
            
            #Check to ensure that the "readfile", which correspond to bam files can be accessed
            fileaccess = check_file_access(readfile)

            #Now build the BASH pipeline.
            if fileaccess is True:

                #Now build the command pipeline:
                INPUT = "bamtools filter -in %s" % readfile 
                OUTPUT = " -out %s" % outputfilename
                MAPQUALITY = " -mapQuality \">= %i \"" % global_vars.EXPERIMENT_SETTINGS['Map Quality']
                
                #Join the commands
                cmd = INPUT + OUTPUT + MAPQUALITY

                #use the subprocess to run the command in shell 
                cmd_output = subprocess.Popen(cmd,shell=True,stderr=subprocess.PIPE)
               
                # Write the output. The addition of b'' and utf-8 is essential to process the command
                while True:
                    out = cmd_output.stderr.read(1)
                    if out == b'' and cmd_output.poll() != None:
                        break
                    if out != '':
                        sys.stdout.write(out.decode('utf-8'))
                        sys.stdout.flush()
                #If its completed, write to log
                #write_to_log("Successfully filtered %s bam file by map quality!" % outputfilename)
            else:
                messagebox.showerror(title="Error",message="This file could not be accessed")
                write_to_log("Error! Could not filter the reads!")
                return False

        #Saving outputfilesdir list to global_vars.MapQualityFileLocation
        global_vars.MapQualityFileLocation = Outputfilesdir

        #write_to_log("Output files are as follows: %s" % global_vars.MapQualityFileLocation)
        return True
    


    def filteringchrM(self,bamfilelist):
        
        #Define output directory
        outputdir= os.path.join(global_vars.EXPERIMENT_SETTINGS['Experiment directory'] + "/")

        #change working directory
        os.chdir(outputdir)
        #write_to_log("The current working directory is %s" % outputdir)

        #Save the list of output file paths in "Outputfilesdir"
        Outputfilesdir = []
        #Pipe the individual bamfile for mapquality filtering.
        for i in bamfilelist:
            
            #So the file can be read
            separator=''
            readfile = separator.join(i)
            #write_to_log("The mouse genome is located %s" % global_vars.genomeLocation)

            #Saving output file name:
            outputfilename= "File" + str(bamfilelist.index(i)+1) + global_vars.File_Flags['Mitochondrial Reads Filtered'] + ".bam"
            write_to_log("Filtering reads for mitochondrial DNA for %s" % readfile)

            #Saving the resulting files in a list called "outputfilesdir". 
            Output_files_dir = outputdir + outputfilename
            Outputfilesdir.append(Output_files_dir)
            
            #Check to ensure that the "readfile", which correspond to bam files can be accessed
            fileaccess = check_file_access(readfile)

            #Now build the BASH pipeline.
            if fileaccess is True:
                #write_to_log("The %s file is now being processed and will be saved as %s!" % (readfile,outputfilename))
                
                #Download and extract the mouse mm10 chromatin
                #write_to_log("Now Removing Mitochondrial Reads")
                cmd1 = "samtools view -@ %s %s |" % (global_vars.EXPERIMENT_SETTINGS['Cores'], readfile)
                chrm     = " egrep -v chrM |"
                mm10     = " samtools view -@ %s -bT %s" % (global_vars.EXPERIMENT_SETTINGS['Cores'], global_vars.genomeLocation)
                output   = " -o %s" % outputfilename

                cmd = cmd1 + chrm + mm10 + output
        
                #use the subprocess to run the command in shell 
                cmd_output = subprocess.Popen(cmd,shell=True,stderr=subprocess.PIPE)
               
                # Write the output. The addition of b'' and utf-8 is essential to process the command
                while True:
                    out = cmd_output.stderr.read(1)
                    if out == b'' and cmd_output.poll() != None:
                        break
                    if out != '':
                        sys.stdout.write(out.decode('utf-8'))
                        sys.stdout.flush()
                #If its completed, write to log
                #write_to_log("Successfully filtered %s bam file by map quality!" % outputfilename)
            else:
                messagebox.showerror(title="Error",message="This file could not be filtered for mitochondrial reads")
                write_to_log("Error! This file could not be filtered for mitochondrial reads!")
                return False

        #Saving outputfilesdir list to global_vars.MapQualityFileLocation
        global_vars.MitochondrialReadsFiltered = Outputfilesdir

        #write_to_log("Output files are as follows: %s" % global_vars.MitochondrialReadsFiltered)
        return True



    def RemoveDups(self,bamfilelist):
        outputdir= os.path.join(global_vars.EXPERIMENT_SETTINGS['Experiment directory'] + "/")
        os.chdir(outputdir)

        Outputfilesdir = []
        for i in bamfilelist:
            #write_to_log("The current working directory is %s" % outputdir)

            #This is to remove the spaces from the file list that was saved under "gobal_vars.MapQualityFileLocation" in previous def
            #Also write to log
            separator=''
            readfile = separator.join(i)

            write_to_log("Now Removing duplicate reads for %s" % readfile)

            #Define output file names
            outfile_bam = "File" + str(bamfilelist.index(i)+1) + global_vars.File_Flags['Remove Duplicates'] + global_vars.File_Flags['sortchrM_Reads'] + ".bam"
            outfile_metrics = "File" + str(bamfilelist.index(i)+1) + global_vars.File_Flags['Remove Duplicates'] + ".txt"

            #write_to_log("The file being processed is %s. Processed file will be saved as %s" % (readfile,outfile_bam))
            
            fileaccess = check_file_access(readfile) #checking file access

            #Saving the resulting files in a list called "outputfilesdir". 
            Output_files_dir = outputdir + outfile_bam
            Outputfilesdir.append(Output_files_dir)

            if fileaccess is True:
                #Building command pipeline:
                #write_to_log("Initiating Picards")
                
                cmd1 = "java -jar %s MarkDuplicates REMOVE_DUPLICATES=true I=%s" % (global_vars.EXPERIMENT_SETTINGS['picard directory'], readfile)
                cmd2 = " O=%s M=%s" % (outfile_bam,outfile_metrics)

                cmd = cmd1 + cmd2
                
                #use the subprocess to run the command in shell 
                cmd_output = subprocess.Popen(cmd,shell=True,stderr=subprocess.PIPE)
               
                # Write the output. The addition of b'' and utf-8 is essential to process the command
                while True:
                    out = cmd_output.stderr.read(1)
                    if out == b'' and cmd_output.poll() != None:
                        break
                    if out != '':
                        sys.stdout.write(out.decode('utf-8'))
                        sys.stdout.flush()
                #If its completed, write to log
                #write_to_log("Successfully removed duplicates the %s file and saved it as %s bam file!" % (readfile, global_vars.DupsRemoved))
            else:
                messagebox.showerror(title="Error",message="Could not remove duplicates from the bam file. Make sure you have java and JDK installed on your computer!")
                write_to_log("Error! Could not remove duplicates from the bam file")
                return False

        global_vars.DupsRemoved =  Outputfilesdir
        return True


    
    def sortchMreads(self,bamfilelist):
        outputdir= os.path.join(global_vars.EXPERIMENT_SETTINGS['Experiment directory'] + "/")
        os.chdir(outputdir)

        Outputfilesdir = []
        for i in bamfilelist:
            
            #write_to_log("The current working directory is %s" % outputdir)

            #This is to remove the spaces from the file list that was saved under "gobal_vars.MapQualityFileLocation" in previous def
            #Also write to log
            separator=''
            readfile = separator.join(i)

            write_to_log("Sorting reads for  %s" % readfile)

            #Define output file names
            outfile = "File" + str(bamfilelist.index(i)+1) + global_vars.File_Flags['sortchrM_Reads'] + ".bam"
            
            fileaccess = check_file_access(readfile) #checking file access

            #Saving the resulting files in a list called "outputfilesdir". 
            Output_files_dir = outputdir + outfile
            Outputfilesdir.append(Output_files_dir)

            if fileaccess is True:
                #Building command pipeline:
                write_to_log("Initiating Samtools")

                cmd1= "samtools sort -@ %s -O bam -o %s" % (global_vars.EXPERIMENT_SETTINGS['Cores'], outfile)
                cmd2= " -T temp %s" % readfile

                cmd = cmd1 + cmd2
                
                #use the subprocess to run the command in shell 
                cmd_output = subprocess.Popen(cmd,shell=True,stderr=subprocess.PIPE)
               
                # Write the output. The addition of b'' and utf-8 is essential to process the command
                while True:
                    out = cmd_output.stderr.read(1)
                    if out == b'' and cmd_output.poll() != None:
                        break
                    if out != '':
                        sys.stdout.write(out.decode('utf-8'))
                        sys.stdout.flush()
                #If its completed, write to log
                #write_to_log("Successfully sorted the %s file and saved it as %s bam file!" % (readfile, global_vars.MitochondrialReadsFiltered))
            else:
                messagebox.showerror(title="Error",message="Could not sort the bam file")
                write_to_log("Error! Could not sort the bam file")
                return False

        global_vars.Sorted_chrM =  Outputfilesdir
        return True



    def indxSortReads(self,bamfilelist):
        outputdir= os.path.join(global_vars.EXPERIMENT_SETTINGS['Experiment directory'] + "/")
        os.chdir(outputdir)

        for i in bamfilelist:
            
            #write_to_log("The current working directory is %s" % outputdir)

            #This is to remove the spaces from the file list that was saved under "gobal_vars.MapQualityFileLocation" in previous def
            #Also write to log
            separator=''
            readfile = separator.join(i)
            write_to_log("Indexing BAM file for viewing in IGV (This generates a .bai file) for %s" % readfile)

            #Define output file names
            outfile = "File" + str(bamfilelist.index(i)+1) + global_vars.File_Flags['Index_BamFile']
            #write_to_log("The file being processed is %s. Processed file will be saved as %s" % (readfile,outfile))
            
            fileaccess = check_file_access(readfile) #checking file access

            if fileaccess is True:
                #Building command pipeline:
                #write_to_log("Indexing the bam file now")

                cmd = "samtools index -@ %s %s %s" % (global_vars.EXPERIMENT_SETTINGS['Cores'],readfile,outfile)

                #use the subprocess to run the command in shell 
                cmd_output = subprocess.Popen(cmd,shell=True,stderr=subprocess.PIPE)
               
                # Write the output. The addition of b'' and utf-8 is essential to process the command
                while True:
                    out = cmd_output.stderr.read(1)
                    if out == b'' and cmd_output.poll() != None:
                        break
                    if out != '':
                        sys.stdout.write(out.decode('utf-8'))
                        sys.stdout.flush()
                #If its completed, write to log
                #write_to_log("Successfully indexed the %s file!" % readfile)
            else:
                messagebox.showerror(title="Error",message="Could not index the bam file")
                write_to_log("Error! Could not index the bam file")
                return False
        return True



    def FragmentSizeDist(self,bamfilelist):
        outputdir= os.path.join(global_vars.EXPERIMENT_SETTINGS['Experiment directory'] + "/")
        os.chdir(outputdir)

        for i in bamfilelist:
            
            #write_to_log("The current working directory is %s" % outputdir)

            #This is to remove the spaces from the file list that was saved under "gobal_vars.MapQualityFileLocation" in previous def
            #Also write to log
            separator=''
            readfile = separator.join(i)
            write_to_log("Generate fragment size distribution graph using picards for %s" % readfile)
            #write_to_log("The bam file being processed is %s" % readfile)

            #Define output file names
            outfile_pdf = "File" + str(bamfilelist.index(i)+1) + global_vars.File_Flags['fragment size distribution'] + ".pdf"
            outfile_txt = "File" + str(bamfilelist.index(i)+1) + global_vars.File_Flags['fragment size distribution'] + ".txt"
            
            #readfile = i
            fileaccess = check_file_access(readfile) #checking file access

            if fileaccess is True:
                #Building command pipeline:
                #write_to_log("Initiating picard suite")

                INPUT   = "java -jar %s CollectInsertSizeMetrics I=%s" % (global_vars.EXPERIMENT_SETTINGS['picard directory'], readfile) 
                OUTPUT  = " O=%s" % outfile_txt
                HIST    = " H=%s M=0.5" % outfile_pdf

                cmd = INPUT + OUTPUT + HIST 
                
                #use the subprocess to run the command in shell 
                cmd_output = subprocess.Popen(cmd,shell=True,stderr=subprocess.PIPE)
               
                # Write the output. The addition of b'' and utf-8 is essential to process the command
                while True:
                    out = cmd_output.stderr.read(1)
                    if out == b'' and cmd_output.poll() != None:
                        break
                    if out != '':
                        sys.stdout.write(out.decode('utf-8'))
                        sys.stdout.flush()
                
                #If its completed, write to log
                #write_to_log("Successfully graphed fragment size distribution %s bam file by map quality!" % global_vars.MitochondrialReadsFiltered)
            else:
                messagebox.showerror(title="Error",message="This file could not be accessed")
                write_to_log("Error! Could not perform fragment size distribution graph!")
                return False
        return True



    def MACS2_peackcalling(self,bamfilelist):
        outputdir= os.path.join(global_vars.EXPERIMENT_SETTINGS['Experiment directory'] + "/")
        os.chdir(outputdir)

        #list of broadpeak and bedgraph files
        BroadPeakFileName   = []
        Bedgraphfilename    = []

        for i in bamfilelist:
            
            #write_to_log("The current working directory is %s" % outputdir)

            #This is to remove the spaces from the file list that was saved under "gobal_vars.MapQualityFileLocation" in previous def
            #Also write to log
            separator=''
            readfile = separator.join(i)
            write_to_log("Calling peaks using MACS2 for %s" % readfile)

            #Define output file names
            outfile = "File" + str(bamfilelist.index(i)+1) + global_vars.File_Flags['Peak_Calls'] 

            #Broadpeak and bedgraph files append:
            peaks_dir           = "PEAKS_%s" % outfile  #peak directory name
            broad_bedpeaks_dir  = os.path.join(global_vars.EXPERIMENT_SETTINGS['Experiment directory'] + "/" + peaks_dir + "/") #peak directory
            broadpeak_file      = str(broad_bedpeaks_dir) + "splitted_peaks.broadPeak"
            BroadPeakFileName.append(broadpeak_file)
            bedgraph_file       = str(broad_bedpeaks_dir) + "splitted_treat_pileup.bdg"
            Bedgraphfilename.append(bedgraph_file)

            fileaccess = check_file_access(readfile) #checking file access
            if fileaccess is True:
                #Building command pipeline:
                write_to_log("Initiating MACS2")

                cmd1 = "macs2 callpeak -t %s" % readfile
                cmd2 = " -q 0.05 --broad -f BAMPE -n splitted -g mm -B --nomodel --shift -100 --extsize 200 --keep-dup all"
                cmd3 = " --outdir PEAKS_%s" % outfile
                cmd = cmd1 + cmd2 + cmd3

                #use the subprocess to run the command in shell 
                cmd_output = subprocess.Popen(cmd,shell=True,stderr=subprocess.PIPE)
               
                # Write the output. The addition of b'' and utf-8 is essential to process the command
                while True:
                    out = cmd_output.stderr.read(1)
                    if out == b'' and cmd_output.poll() != None:
                        break
                    if out != '':
                        sys.stdout.write(out.decode('utf-8'))
                        sys.stdout.flush()
                #If its completed, write to log
                write_to_log("Successfully called peaks on the %s file and saved it in the MACS2%s folder!" % (readfile, global_vars.MitochondrialReadsFiltered))
            else:
                messagebox.showerror(title="Error",message="Could not call peaks!")
                write_to_log("Error! Could not call peaks!")
                return False

        #Save the broadpeak file locations for conversion to bed files!
        global_vars.broadpeak_files  = BroadPeakFileName
        global_vars.bedgraphlist     = Bedgraphfilename

        return True



    def broadpeak_bed_conversion(self, bamfilelist):
        outputdir= os.path.join(global_vars.EXPERIMENT_SETTINGS['Experiment directory'] + "/")
        os.chdir(outputdir)

        for i in bamfilelist:
            #write_to_log("The current working directory is %s" % outputdir)

            #This is to remove the spaces from the file list that was saved under "gobal_vars.MapQualityFileLocation" in previous def
            #Also write to log
            separator=''
            readfile = separator.join(i)
            write_to_log("Generating BED file from broadpeak file for %s" % readfile)

            #Define output file names
            outfile_name = "File" + str(bamfilelist.index(i)+1) + global_vars.File_Flags['Bed File'] + ".bed"
            
            fileaccess = check_file_access(readfile) #checking file access
            
            if fileaccess is True:
                #Building command pipeline:
                #write_to_log("Generating BED files.....")

                cmd = "cut -f 1-6 %s > %s" % (readfile, outfile_name)

                #use the subprocess to run the command in shell 
                cmd_output = subprocess.Popen(cmd,shell=True,stderr=subprocess.PIPE)
               
                # Write the output. The addition of b'' and utf-8 is essential to process the command
                while True:
                    out = cmd_output.stderr.read(1)
                    if out == b'' and cmd_output.poll() != None:
                        break
                    if out != '':
                        sys.stdout.write(out.decode('utf-8'))
                        sys.stdout.flush()
                #If its completed, write to log
                #write_to_log("Successfully called peaks on the %s file and saved it in the MACS2%s folder!" % (readfile, global_vars.MitochondrialReadsFiltered))
            else:
                messagebox.showerror(title="Error",message="Could not create bed files!")
                write_to_log("Error! Could not create bed files!")
                return False

        return True



    def RemoveFiles(self, bamfilelist):
        #write_to_log("This is function to remove unnecessary files from the working directory")
        
        #Define output directory
        outputdir= os.path.join(global_vars.EXPERIMENT_SETTINGS['Experiment directory'] + "/")

        #change working directory
        os.chdir(outputdir)

        for i in bamfilelist:

            write_to_log("Now removing unnecessary files....")

            #Unnecessary BAM files
            MapQuality = "File" + str(bamfilelist.index(i)+1) + global_vars.File_Flags['mapquality'] + ".bam"
            chrM       = "File" + str(bamfilelist.index(i)+1) + global_vars.File_Flags['Mitochondrial Reads Filtered'] + ".bam"
            sort_chrM  = "File" + str(bamfilelist.index(i)+1) + global_vars.File_Flags['sortchrM_Reads'] + ".bam"

            #Unncessary TXT files:
            dupsremoved_txt = "File" + str(bamfilelist.index(i)+1) + global_vars.File_Flags['Remove Duplicates'] + ".txt"
            fragSize_txt    = "File" + str(bamfilelist.index(i)+1) + global_vars.File_Flags['fragment size distribution'] + ".txt"

            cmd = "rm %s %s %s %s %s" % (MapQuality,chrM,sort_chrM,dupsremoved_txt,fragSize_txt)
            
            cmd_output=subprocess.Popen(cmd,shell=True, stderr=subprocess.PIPE)

            while True:
                out = cmd_output.stderr.read(1)
                if out == b'' and cmd_output.poll() != None:
                    break
                if out != '':
                    sys.stdout.write(out.decode('utf-8'))
                    sys.stdout.flush()
                    #write_to_log("Removed all unnecessary files from the working directory")
                else:
                    messagebox.showerror(title="Error",message="Could not delete files")
                    write_to_log("Could not delete unnecssary files")
                    return False
        return True



    def At_chr(self):
        outputdir = os.path.join(global_vars.EXPERIMENT_SETTINGS['Experiment directory'] + "/")
        os.chdir(outputdir)
        
        write_to_log("Generating the At_Chr file for big wig conversion. This uses Bioawk...")
        #Define output file names
        outfile_name = str(global_vars.File_Flags['At_chr sizes file'])
        #write_to_log("The file being processed is %s. Processed file will be saved as %s" % (readfile,outfile_name))
        genome = global_vars.genomeLocation
        
        cmd = "%s -c fastx \'{print $name, length($seq)}\' %s > %s " % (global_vars.EXPERIMENT_SETTINGS['Bioawk Directory'], global_vars.genomeLocation, outfile_name)
        
        #use the subprocess to run the command in shell 
        cmd_output = subprocess.Popen(cmd,shell=True,stderr=subprocess.PIPE)
        
        # Write the output. The addition of b'' and utf-8 is essential to process the command
        while True:
            out = cmd_output.stderr.read(1)
            if out == b'' and cmd_output.poll() != None:
                break
            if out != '':
                sys.stdout.write(out.decode('utf-8'))
                sys.stdout.flush()
        
        #save directory for chromosome sizes:
        global_vars.chromosome_sizeDir = os.path.join(outputdir + outfile_name)



    def bedgraphtobigwig(self, bedgraphlist):
        outputdir= os.path.join(global_vars.EXPERIMENT_SETTINGS['Experiment directory'] + "/")
        os.chdir(outputdir)

        #list of broadpeak and bedgraph files
        clipped_bdg_file = []

        for i in bedgraphlist:
            
            separator=''
            readfile = separator.join(i)
            write_to_log("Clipping bedgraph File for %s" % readfile)
            
            #output file name:
            outfile_name = "File" + str(bedgraphlist.index(i)+1) + global_vars.File_Flags['Clipped Bed Graph File'] + ".bdg"
            
            #Saving the resulting files in a list called "outputfilesdir". 
            Output_files_dir = outputdir + outfile_name
            clipped_bdg_file.append(Output_files_dir)

            fileaccess = check_file_access(readfile) #checking file access
            if fileaccess is True:
                #Building command pipeline:
                cmd1 = "bedtools slop -i %s -g %s -b 0|" % (readfile, global_vars.chromosome_sizeDir)
                cmd2 = " bedClip stdin %s %s" % (global_vars.chromosome_sizeDir, outfile_name)
                cmd = cmd1 + cmd2
                cmd_output = subprocess.Popen(cmd,shell=True,stderr=subprocess.PIPE)
               
                # Write the output. The addition of b'' and utf-8 is essential to process the command
                while True:
                    out = cmd_output.stderr.read(1)
                    if out == b'' and cmd_output.poll() != None:
                        break
                    if out != '':
                        sys.stdout.write(out.decode('utf-8'))
                        sys.stdout.flush()
            else:
                messagebox.showerror(title="Error",message="Could not clip the bedgraph file!")
                write_to_log("Error! Could not proceed!")
                return False

        #Save the clipped bedgraph files in global_vars
        global_vars.clippedBedGraphFiles = clipped_bdg_file
        
        return True



    def sortclippedbedgraph(self, bedgraphlist):
        outputdir= os.path.join(global_vars.EXPERIMENT_SETTINGS['Experiment directory'] + "/")
        os.chdir(outputdir)

        #list of broadpeak and bedgraph files
        sorted_clipped_bedgraph = []

        for i in bedgraphlist:
            
            separator=''
            readfile = separator.join(i)
            write_to_log("Sorting the clipped bedgraph File for %s" % readfile)

            #output file name:
            outfile_name = "File" + str(bedgraphlist.index(i)+1) + global_vars.File_Flags['Clipped and sorted Bedgraph file'] + ".bdg"
            
            #Saving the resulting files in a list called "outputfilesdir". 
            Output_files_dir = outputdir + outfile_name
            sorted_clipped_bedgraph.append(Output_files_dir)

            fileaccess = check_file_access(readfile) #checking file access
            if fileaccess is True:
                #Building command pipeline:
                cmd = "sort -k1,1 -k2,2n %s > %s" % (readfile,outfile_name)
                cmd_output = subprocess.Popen(cmd,shell=True,stderr=subprocess.PIPE)
               
                # Write the output. The addition of b'' and utf-8 is essential to process the command
                while True:
                    out = cmd_output.stderr.read(1)
                    if out == b'' and cmd_output.poll() != None:
                        break
                    if out != '':
                        sys.stdout.write(out.decode('utf-8'))
                        sys.stdout.flush()
            else:
                messagebox.showerror(title="Error",message="Could not call sort the clipped bedgraph file!")
                write_to_log("Error! Could not proceed!")
                return False

        #Save the clipped bedgraph files in global_vars
        global_vars.sorted_clipped_bedgraph = sorted_clipped_bedgraph
        
        return True



    def BigWigConversion(self, bedgraphlist):
        outputdir= os.path.join(global_vars.EXPERIMENT_SETTINGS['Experiment directory'] + "/")
        os.chdir(outputdir)

        #list of broadpeak and bedgraph files
        BigWigFileList = []

        for i in bedgraphlist:
            separator=''
            readfile = separator.join(i)
            write_to_log("Generating BigWig file from the sorted and clipped bedgraph file for %s" % readfile)
            
            #output file name:
            outfile_name = "File" + str(bedgraphlist.index(i)+1) + global_vars.File_Flags['BigWig File'] + ".bw"
            
            #Saving the resulting files in a list called "outputfilesdir". 
            Output_files_dir = outputdir + outfile_name
            BigWigFileList.append(Output_files_dir)

            fileaccess = check_file_access(readfile) #checking file access
            if fileaccess is True:
                #Building command pipeline:
                cmd = "bedGraphToBigWig %s %s %s" % (readfile,global_vars.chromosome_sizeDir,outfile_name)
                cmd_output = subprocess.Popen(cmd,shell=True,stderr=subprocess.PIPE)
               
                # Write the output. The addition of b'' and utf-8 is essential to process the command
                while True:
                    out = cmd_output.stderr.read(1)
                    if out == b'' and cmd_output.poll() != None:
                        break
                    if out != '':
                        sys.stdout.write(out.decode('utf-8'))
                        sys.stdout.flush()
            else:
                messagebox.showerror(title="Error",message="Could not create bigwig file!")
                write_to_log("Error! Could not proceed!")
                return False

        #Save the clipped bedgraph files in global_vars
        global_vars.BigWigFileList = BigWigFileList
        
        return True



    def finalcleanup(self, bedgraphlist):
        #write_to_log("This is function to remove unnecessary files from the working directory")
        
        #Define output directory
        outputdir= os.path.join(global_vars.EXPERIMENT_SETTINGS['Experiment directory'] + "/")

        #change working directory
        os.chdir(outputdir)

        for i in bedgraphlist:

            write_to_log("Now removing unnecessary files....")

            #Unnecessary BAM files
            chr           =  global_vars.File_Flags['At_chr sizes file']
            clipped       = "File" + str(bedgraphlist.index(i)+1) + global_vars.File_Flags['Clipped Bed Graph File'] + ".bdg"
            sort_clipped  = "File" + str(bedgraphlist.index(i)+1) + global_vars.File_Flags['Clipped and sorted Bedgraph file'] + ".bdg"

            cmd = "rm %s %s %s" % (chr,clipped,sort_clipped)
            
            cmd_output=subprocess.Popen(cmd,shell=True, stderr=subprocess.PIPE)

            while True:
                out = cmd_output.stderr.read(1)
                if out == b'' and cmd_output.poll() != None:
                    break
                if out != '':
                    sys.stdout.write(out.decode('utf-8'))
                    sys.stdout.flush()
                    #write_to_log("Removed all unnecessary files from the working directory")
                else:
                    messagebox.showerror(title="Error",message="Could not delete files")
                    write_to_log("Could not delete unnecssary files")
                    return False
        return True



#Close the mainloop to wrap the software:
if __name__ == "__main__":
    app = ProcessingTool()
    app.geometry("800x570+300+50")
    app.resizable(FALSE, TRUE)
    app.mainloop()
    _thread.start_new_thread(app.mainloop(), ())