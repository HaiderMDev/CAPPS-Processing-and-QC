import tkinter as tk
from tkinter import font as tkfont
from tkinter import *
from tkinter import filedialog
from tkinter import messagebox
from tkinter.ttk import Button, Progressbar, Style
import os
import subprocess
import sys
from turtle import bgcolor, end_fill, write
import pandas as pd
import numpy as np
from scipy import stats
from datetime import datetime
import _thread
import globalvariables as global_vars

#Go to this webpage for color selection: https://encycolorpedia.com/3d6fdb

# version of the application
app_version = "1.8"

#function to detect display
if os.environ.get('DISPLAY','') == '':
    print('no display found. Using :0.0')
    os.environ.__setitem__('DISPLAY', '127.0.0.1:0.0')


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
        #write to the log file
        #with the a option, it opens a file for appending at the end of the file without truncating it, 
        #or it creates a new file if it does not exist.
        #The "+" option opens a file for updating (reading and writing)
        file=open(file_name,'a+')
        #%s specifically is used to perform concatenation of strings together
        #The % symbol tells the program what add where you use "%s" symbol
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


#SELF represents the current object. This is a common first parameter for any method of a class. As you suggested, 
#it's similar to Java's this.

#PARENT represents a widget to act as the parent of the current object. All widgets in tkinter except the root window 
#require a parent (sometimes also called a master)

#controller represents some other object that is designed to act as a common point of interaction for several pages of 
#widgets. It is an attempt to decouple the pages. That is to say, each page doesn't need to know about the other pages. 
#If it wants to interact with another page, such as causing it to be visible, it can ask the controller to make it visible.

#You asked "There is a function already defined called show_frame, but why is the controller being used to call this 
#function?" Notice that show_frame is defined in a separate class, in this case the main program class. It is not defined 
#in the other classes. For the other classes to be able to call it, they must call it on an instance of the main class. 
#That instance is named controller in the context of these other classes.


class ProcessingTool(tk.Tk):
    #Initiate the frame.
    #you need __init__ in order to write that code that needs to get executed 
    # for every new object whenever this object gets initialized / created - not just once when the class is read in.
    #you need __init__ in order to write that code that needs to get executed 
    # for every new object whenever this object gets initialized / created - not just once when the class is read in.
    #All of this is object oriented programing 
    
    def __init__(self, *args, **kwargs):
        #*args and **kwargs allow you to pass an unspecified 
        #number of arguments to a function, so when writing the function definition, you do not need 
        #to know how many arguments will be passed to your function.

        tk.Tk.__init__(self, *args, **kwargs)
        #So with this __init__ function, you are simply setting the shape of the window, stypes, fonts, title of the window
        #page frames and what page will be displayed first. 
        # load global variables. Initiate the global variables. 
       
        global_vars.init_vars().__init__             #Load global variables

        #Main container. This is for the section that everything will be displayed.
        #This sets the color of the frame, and how the items are attached to the grid. 
        container = tk.Frame(self, bg="#8c8c8c")
        container.pack(side="top", fill="both", expand=TRUE)
        container.grid_rowconfigure(0, weight=2)
        container.grid_columnconfigure(0, weight=1)

        # widget styles. This is settings for the widgets that will be displayed in the tkinter windows!
        self.style = Style()                
        self.style.theme_use('clam')     #For windows you can use: vista, clam, alt, classic, default, winnative,xpnative 
        self.frame_relief = GROOVE          #Options are RAISED, FLAT, SUNKEN, GROOVE, RIDGE
        self.frame_borderwidth = 2          #Determines how raised should the widgets be 
       
        # This creates the title at the top of the window
        self.winfo_toplevel().title("CAPPS " + app_version)

        # set fonts
        self.title_font = tkfont.Font(family='Arial', size=20, weight="bold")    #Font for the main title, make the weight normal
        self.label_font = tkfont.Font(family='Arial', size=14, weight="bold")    #Font for labels
        self.maintext_font=tkfont.Font(family='Arial',size=18, weight="bold")    #Font for MainText
        self.picking_text=tkfont.Font(family="Arial", size=15, weight="normal")  #font for picking frame                 

        #Define page frames. This will depend on the number of pages you are planning on coding for as classes!
        self.frames = {}
        #create frames for each page.
        # put all of the pages in the same location;
        # the one on the top of the stacking order will be the one that is visible.
        # draw all frames in same location.
        #nsew, specifying which edges of the cell the widget should be "stuck" to. For example, 
        #a value of n (north) will jam the widget up against the top side, with any extra vertical 
        #space on the bottom; the widget will still be centered horizontally. A value of nw (north-west) 
        #means the widget will be stuck to the top left corner, with extra space on the bottom and right.
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


        #file= filedialog.askopenfilenames(initialdir=os.getcwd(),title="Select your files",filetypes=(file_options))
        #for item in file:
           #target_obj.delete(0, END)
           #target_obj.insert(END,item)


    # function to delete selected file from file list
    def delete_selected(self, frame_name, target_name):
        target_obj = getattr(self.frames[frame_name], target_name)
        # curselection is a predefined function that fetches the value(s) of a selected item or items.
        selected = target_obj.curselection() 
        #For the [::-1] notation:
        #It iterates over any iterable thing, backward. That is an example of 
        # Slice notation. The full notation is [ start: stop: stride ]. When start 
        # is left out it defaults to zero. When stop is left out it defaults to last+1. A negative 
        # stride goes from last backwards to start. -2 would output every other value. -1 outputs every value. 
        # If stride is left out (which is common) it defaults to 1. So [::] just produces a copy of all the 
        # elements in order. [::-1] produces a copy of all the elements in reverse order.
        for i in selected[::-1]:
            target_obj.delete(i)


 
class MainPage(tk.Frame):
    def __init__(self, parent, controller):
        #Regarding self, parent, controller:
        # they're just two other initialisation parameters. 
        # In Tkinter, you generally pass the parent widget within which the new widget 
        # sits to each new widget to define a tree for the whole UI. 
        # Controller is apparently something that can be used to control the overall UI, 
        # rather than make a widget responsible for global changes. 
        #The parent and controller were identified in the ATACseqAlign class under definition of __init__
        
        tk.Frame.__init__(self, parent) #Initialize the new tkinter frame under the parent

        self.controller = controller    #Not necessary

        #This is controls for setting up the frame: 
        frame = tk.Frame(self, relief=GROOVE, borderwidth=15, bg="#8C8C8C")
        frame.pack(fill=BOTH, padx=5, pady=5, expand=TRUE) 
        #If expand = FALSE, the frame won't take up the whole page. 
        #pady and padx: How many pixels to pad widget, horizontally and vertically, outside v's borders.

        #info for the main TOP WIDGET which contains the name of the tool in mainpage 
        frame1 = tk.Frame(frame, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth, bg="#000000")
        frame1.pack(fill=X, padx=5, pady=5, anchor=N, expand=TRUE)
        label1_info = tk.Label(frame1, text="ChIP-seq and ATAC-seq Filtering & Peak calling Software (CAPPS)",
                         wraplength="800", justify="left", font=controller.title_font, bg="#000000", fg="white")
        label1_info.pack(padx=5, pady=5)

        #Settings for the second widget frame which contains information about the tool in mainpage
        frame2 = tk.Frame(frame, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth, bg="#3d6fdb")
        frame2.pack(fill=X, side=TOP, padx=5, pady=5, anchor=NW, expand=True)
        label2_info = tk.Label(frame2, text=
                                    "INTRODUCTION"
                                    "\n CAPPS tool is designed for processing and peak calling for ATAC-seq and ChIP-seq BAM files."
                                    " This tool takes 'BAM' files as input and relies on MACS2, Bamtools, Samtools and Picards for processing."
                                    " Please ensure these packages are installed on your computer."
                                    " Data visualization is available separately in R."
                                    "\n\n OUTPUT"
                                    "\n This tool generates indexed, sorted and de-duplicated BAM file(s), fragment size distribution graph, "
                                    " as well as, called peak(s) and bed file(s). There is no limit on the number of BAM files that can be processed simultaneously."
                                    "\n\n IMPORTANT"
                                    "\n Program will generate an error if there are spaces in directories names."
                                    " For example: './User Directory' will generate an error, but './User_Directory' or './UserDirectory' is acceptable.",
                         wraplength="650", justify="center", font=controller.picking_text, bg="#3d6fdb", fg="white")
        #Setting for the lable which is contained in the fr_2 frame. This one contains the version of the tool:
        label2_info.pack(side=TOP, padx=5, pady=5)
        label_extra = tk.Label(frame2, text="\n BAMing Tool version " + app_version,
                            wraplength="500", justify="center", font=controller.maintext_font, bg="#3d6fdb", fg="white")
        label_extra.pack(side=TOP, padx=5, pady=5)

        #Button to move to the next frame, which is where the user selects experimental settings. 
        #if you pick "self" as the frame, then it will put the button at the bottom of the tk window. 
        #You can alternatively use "frame" as the option, then it will place the button in the frame you defined. 
        btnNewExp = Button(self, text="Experiment Information", 
        command=lambda: controller.show_frame("information"))
        btnNewExp.pack(side=RIGHT, padx=5, pady=1)
        self.pack(fill=BOTH, expand=TRUE)

        #Button to exit the application
        btnBack = Button(self, text="Exit", command=lambda: app.destroy())
        btnBack.pack(side=LEFT, padx=5, pady=5)



class information(tk.Frame):
    def __init__(self, parent, controller):
        
        #Load the frame:
        tk.Frame.__init__(self,parent)
        self.controller = controller
        
        #Set up the frame for the body:
        frame_information=tk.Frame(self, relief=FLAT, borderwidth=controller.frame_borderwidth,bg="#8c8c8c")
        frame_information.pack(fill=BOTH,padx=5,pady=5,expand=TRUE)

        #Heading for the Frame 
        frame1_header=tk.Frame(frame_information,relief=GROOVE,borderwidth=controller.frame_borderwidth,bg="#8c8c8c")
        frame1_header.pack(fill=BOTH,pady=5,padx=5,anchor=N)
        frame1_header_label=tk.Label(frame_information,text="Define Experimental Information",bg="#000000", font=controller.title_font, fg="white") #fg is foreground color, and bg is background color
        frame1_header_label.pack(fill=X,pady=5,padx=5,anchor=N,expand=TRUE) #If you set "fill = BOTH", then the label will take the whole page

        #Select the Experiment Name
        frame1=tk.Frame(frame_information, relief=controller.frame_relief,borderwidth=controller.frame_borderwidth,bg="#8c8c8c")
        frame1.pack(fill=X, pady=5,padx=5,expand=TRUE,anchor=N)
        label_exp_name=tk.Label(frame1,text="Define Experiment Name", width="40", font=controller.maintext_font,bg="#8c8c8c",fg="black")
        label_exp_name.pack(side=LEFT,anchor=E,pady=5,padx=5)
        self.txt_exp_name=Entry(frame1, width="40")
        self.txt_exp_name.pack(fill=X,padx=5,pady=5,expand=TRUE)

        #Select output Directory
        frame2 = tk.Frame(frame_information,relief=controller.frame_relief, borderwidth=controller.frame_borderwidth,bg="#8c8c8c")
        frame2.pack(fill=X, pady=5,padx=5,anchor=N,expand=TRUE)
        label_exp_folder=tk.Label(frame2,text="Experimental Directory",width="40", font=controller.maintext_font,bg="#8c8c8c",fg="black")
        label_exp_folder.pack(side=LEFT,anchor=N, pady=5,padx=5)
        #entry for frame2
        self.txt_expdir=Entry(frame2,width="40")
        self.txt_expdir.pack(fill=X,pady=5,padx=5,expand=TRUE)
        #Button for selecting diretory:
        btn_browse = Button(frame2,text="BROWSE",command=lambda:self.controller.browse_folder(
            frame_name="information",selection_name="txt_expdir"
        )) 
        #Using browse_folder definition under main class ATACseqAlign. "frame_name" has to be the current class "information". 
        #"target_name" has to be the "txt_expdir" or it won't copy the directory that was picked. 
        btn_browse.pack(fill=X,side=RIGHT,padx=5,pady=5)

        #Map quality score to use
        frame3=tk.Frame(frame_information, relief=controller.frame_relief,borderwidth=controller.frame_borderwidth,bg="#8c8c8c")
        frame3.pack(fill=X, pady=5,padx=5,expand=TRUE,anchor=N)
        label_map_quality=tk.Label(frame3,text="Choose Map Quality Score (Between 20 - 35)", width="40", font=controller.maintext_font,bg="#8c8c8c",fg="black")
        label_map_quality.pack(side=LEFT,anchor=N,pady=5,padx=5)
        self.txt_mapquality=Entry(frame3, width="40")
        self.txt_mapquality.pack(fill=X,padx=5,pady=5,expand=TRUE)

        #Cores to use
        frame5=tk.Frame(frame_information, relief=controller.frame_relief,borderwidth=controller.frame_borderwidth,bg="#8c8c8c")
        frame5.pack(fill=X, pady=5,padx=5,expand=TRUE,anchor=N)
        label_cores=tk.Label(frame5,text="Choose the number of CPU cores", width="40", font=controller.maintext_font,bg="#8c8c8c",fg="black")
        label_cores.pack(side=LEFT,anchor=N,pady=5,padx=5)
        self.txt_cores=Entry(frame5, width="40")
        self.txt_cores.pack(fill=X,padx=5,pady=5,expand=TRUE)

        #Selecting directory for mouse/human genome for filtering purposes:
        frame4 = tk.Frame(frame_information,relief=controller.frame_relief, borderwidth=controller.frame_borderwidth,bg="#8c8c8c")
        frame4.pack(fill=X, pady=5,padx=5,anchor=N,expand=TRUE)
        label_exp_folder2=tk.Label(frame4,text="Select Genome (eg. mm10.fa)",width="40", font=controller.maintext_font,bg="#8c8c8c",fg="black")
        label_exp_folder2.pack(side=LEFT,anchor=N, pady=5,padx=5)
        #Create space or entry
        self.txt_expdir2=Entry(frame4,width="40")
        self.txt_expdir2.pack(fill=X,pady=5,padx=5,expand=TRUE)

        #Button for selecting diretory:
        btn_browse = Button(frame4,text="BROWSE",command=lambda:self.controller.browse_file2(
                                                                            multi_files="TRUE",
                                                                            frame_name="information",
                                                                            selection_name="txt_expdir2",
                                                                            file_type="fa" 
        )) 
        #Using browse_folder definition under main class ATACseqAlign. "frame_name" has to be the current class "information". 
        #"target_name" has to be the "txt_expdir" or it won't copy the directory that was picked. 
        btn_browse.pack(fill=X,side=RIGHT,padx=5,pady=5)

        #Navigation Buttons
        btnexit=Button(self, text="Exit",command=lambda: app.destroy())
        btnexit.pack(side=LEFT,padx=5,pady=5)
        #Button for next page
        btnnextpage=Button(self, text="Next Page", command=lambda:self.store_variables()) #"Store Variables" is defined below!
        #btnnextpage=Button(self, text="Next Page", command=lambda:controller.show_frame("bamfileSelect"))
        btnnextpage.pack(side=RIGHT,pady=5,padx=5)
        #Button to return to the last page
        btnprevious=Button(self, text="Main Page", command=lambda:controller.show_frame("MainPage"))
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
        else:
            error_log = error_log + "\n Please select a genome for filtering"
            
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
                error_log = error_log + "\n CPU cores selection must be identified as a number!"

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
        frame_bamfileSelect=tk.Frame(self, relief=FLAT, borderwidth=controller.frame_borderwidth,bg="#8c8c8c")
        frame_bamfileSelect.pack(fill=BOTH,padx=5,pady=5, expand=TRUE)


        #Heading for the Frame 
        frame1_header=tk.Frame(frame_bamfileSelect,relief=GROOVE,borderwidth=controller.frame_borderwidth,bg="#8c8c8c")
        frame1_header.pack(fill=X,pady=5,padx=5,anchor=N)
        frame1_header_label=tk.Label(frame_bamfileSelect,text="Load BAM Files",bg="#000000", font=controller.title_font, fg="white") #fg is foreground color, and bg is background color
        frame1_header_label.pack(fill=X,pady=5,padx=5,anchor=N,expand=TRUE) 

        #Frame for page banner:
        #frame1_banner=tk.Frame(frame_bamfileSelect,relief=controller.frame_relief, borderwidth=controller.frame_borderwidth,bg="#8c8c8c")
        #frame1_banner.pack(fill=X, padx=5, pady=5,anchor=N)
        #image_banner=PhotoImage(file= os.getcwd() + "/dnastrand.png")
        #lable_image=tk.Label(frame1_banner, image=image_banner,bg="#3d6fdb")
        #lable_image.image = image_banner
        #lable_image=tk.Label(frame1_banner,bg="#3d6fdb")
        #lable_image.pack(padx=5, pady=5, expand=True)

        #Frame for loading bam files:
        frame4=tk.Frame(frame_bamfileSelect,relief=GROOVE,borderwidth=controller.frame_borderwidth,bg="#8c8c8c")
        frame4.pack(fill=X,pady=5,padx=5,anchor=N)
        frame1_header_label=tk.Label(frame_bamfileSelect,text="Select Bam Files",bg="#8c8c8c", width="20", font=controller.maintext_font, fg="black")
        frame1_header_label.pack(fill=X,pady=5,padx=5,anchor=N,side=LEFT) 

        #Frame for list container and buttons
        #Create a frame which will house the container that is scrollable. 
        listcontainer = tk.Frame(frame_bamfileSelect, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth, bg="#000000")
        listcontainer.pack(fill=BOTH, padx=5, pady=5, anchor=N, expand=TRUE)
        #Create a frame within the listcontainer_scroll which will contain the scroll options. For this frame use "fill = BOTH" for packing. 
        listcontainer_scroll = tk.Frame(listcontainer, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth, bg="#8c8c8c")
        listcontainer_scroll.pack(fill=BOTH, padx=5, pady=5, side=RIGHT, expand=TRUE) #If you don't select "fill=BOTH" it won't fill the whole listcontainer frame
        #Button for container 
        list_btn_container = tk.Frame(listcontainer_scroll, bg="#8c8c8c")
        list_btn_container.pack(fill=BOTH, padx=0, pady=0, side=RIGHT)

        #Listbox and scrolling 
        self.list_bamfiles = Listbox(listcontainer_scroll,height="16",width="500", selectmode=EXTENDED)
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
        button_add=Button(list_btn_container,text="Add Files",
                            command=lambda:self.controller.browse_file2(
                                multi_files="TRUE",
                                frame_name="bamfileSelect",
                                selection_name="list_bamfiles",
                                file_type="bam" 
                            ))
        button_add.pack(side=TOP,padx=0,pady=0)

        #Button to remove the files selected
        button_remove=Button(list_btn_container,text="Remove Files",command=lambda:self.controller.delete_selected(
            frame_name="bamfileSelect",
            target_name="list_bamfiles"
        ))
        button_remove.pack(side=BOTTOM,padx=0,pady=0)

        #Button to clear the selection:
        button_clear=Button(list_btn_container,text="Clear Selecton",command=lambda:self.list_bamfiles.delete(0,END))
        button_clear.pack(side=BOTTOM,padx=0,pady=0)

        #Navigation Buttons
        btnexit=Button(self, text="Exit",command=lambda: app.destroy())
        btnexit.pack(side=LEFT,padx=5,pady=5)
        #Button for next page
        btnnextpage=Button(self, text="Summary Page", command=lambda:self.bamfile_dir())
        btnnextpage.pack(side=RIGHT,pady=5,padx=5)
        #Button to return to the last page
        btnprevious=Button(self, text="Previous", command=lambda:controller.show_frame("information"))
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

        main=tk.Frame(self,relief=controller.frame_relief,borderwidth=controller.frame_borderwidth,bg="#8c8c8c")
        main.pack(fill=BOTH,pady=5,padx=5,expand=TRUE)

        #Heading for the frame
        heading=tk.Frame(main,relief=controller.frame_relief,borderwidth=controller.frame_borderwidth,bg="#8c8c8c")
        heading.pack(fill=BOTH,anchor=N,pady=5,padx=5)
        heading_label=tk.Label(heading,text="Summary of Experimental Settings",fg="white",bg="#000000",font=controller.title_font)
        heading_label.pack(fill=BOTH,anchor=N,pady=5,padx=5)

        summaryframe=tk.Frame(main,relief=controller.frame_relief,borderwidth=controller.frame_borderwidth,bg="#8c8c8c")
        summaryframe.pack(fill=BOTH,pady=5,padx=5,anchor=N,expand=TRUE)
        self.summary_txt=Text(summaryframe,height="20")
        self.summary_txt.pack(fill=BOTH,pady=5,padx=5,side=LEFT,expand=TRUE)
        #y-scroll
        yscrol_summary_txt=Scrollbar(self.summary_txt,orient=VERTICAL)
        yscrol_summary_txt.pack(side=RIGHT,fill=Y)
        yscrol_summary_txt.config(command=self.summary_txt.yview)
        self.summary_txt.config(yscrollcommand=yscrol_summary_txt.set)

        #Button to print the Experimental Settings and BAM file directories:
        btnPrintFiles=Button(self,text="Print Summary",command=lambda:self.print_summary()) #have to use self, because the definition for "print_files" is in this frame
        btnPrintFiles.pack(side=TOP,padx=5,pady=5)

        #Navigation Buttons
        btnexit=Button(self, text="Exit",command=lambda: app.destroy())
        btnexit.pack(side=LEFT,padx=5,pady=5)
        #Button for next page
        btnnextpage=Button(self, text="Next Page", command=lambda:controller.show_frame('Filtering_PeakCalling'))
        btnnextpage.pack(side=RIGHT,pady=5,padx=5)
        #Button to return to the last page
        btnprevious=Button(self, text="Previous", command=lambda:controller.show_frame("bamfileSelect"))
        btnprevious.pack(side=RIGHT,pady=0,padx=0)
    
    #So, what we are filling in the global variables.EXPERIMENTAL SETINGS is actually a dictionary. 
    #Python dictionary is a collection of key-value pairs where each key is associated with a value.
    #A value in the key-value pair can be a number, a string, a list, a tuple, or even another dictionary. 
    # In fact, you can use a value of any valid type in Python as the value in the key-value pair.
    #A key in the key-value pair must be immutable. In other words, the key cannot be changed, for example, 
    # a number, a string, a tuple, etc.
    #Python uses the curly braces {} to define a dictionary. 
    #Inside the curly braces, you can place zero, one, or many key-value pairs.
    #So the experimental setting is a dictioary with keys (i.e. "Experiment name", "CPU cores" etc), and the entry from the user is
    #value. 
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
        mainATAC = tk.Frame(self, relief=GROOVE, borderwidth=controller.frame_borderwidth, bg="#8C8C8C")
        mainATAC.pack(fill=BOTH, padx=5, pady=5, expand=True)

        #Heading
        header = tk.Frame(mainATAC, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth, bg="#000000")
        header.pack(fill=X, padx=5, pady=5, anchor=N)
        lbl_header = tk.Label(mainATAC, text="BAM File Filtering and Peak Calling", font=controller.title_font, bg="#000000", fg="white")
        lbl_header.pack(fill=X,padx=5, pady=5, expand=True,anchor=N)

        #Frame for progress bar
        progressbar=tk.Frame(mainATAC,relief=controller.frame_relief,borderwidth=controller.frame_borderwidth,bg="#8c8c8c")
        progressbar.pack(padx=5,pady=5,expand=TRUE,anchor=N)
        #The container is the parent component of the progressbar. I called it "progressbar" as above
        # The orient can be either 'horizontal' or 'vertical
        # The length represents the width of a horizontal progress bar or the height of a vertical progressbar
        # The mode can be either 'determinate' or 'indeterminate'.
        s = Style()
        s.theme_use('clam')
        s.configure("red.Horizontal.TProgressbar",background='red', foreground='red')
        global_vars.progress_bar = Progressbar(progressbar, style= "red.Horizontal.TProgressbar",orient="horizontal", mode="determinate", length="600")
        global_vars.progress_bar.pack(fill=BOTH, padx=5, pady=5, expand=True)
        global_vars.progress_bar['maximum'] = 100
        global_vars.progress_bar['value']= 0

        #Frame for status updates. The status updates will be displayed on this frame. 
        status= tk.Frame(mainATAC, relief=controller.frame_relief, borderwidth=controller.frame_borderwidth, bg="#8c8c8c")
        status.pack(fill=BOTH, padx=5, pady=5, side=TOP, anchor=N, expand=True)
        global_vars.run_log = Text(status, height="20")
        global_vars.run_log.pack(fill=BOTH, side=LEFT, padx=5, pady=5, expand=True)
        scrl_txt_summary = Scrollbar(status)
        scrl_txt_summary.pack(side=RIGHT, fill=Y)
        scrl_txt_summary.config(command=global_vars.run_log.yview)
        global_vars.run_log.config(yscrollcommand=scrl_txt_summary.set)

        #Navigation Buttons
        #Exit Button
        btnexit=Button(self, text="Exit",command=lambda: app.destroy())
        btnexit.pack(side=LEFT,padx=5,pady=5)
        #Back Button:
        global_vars.PreviousButton=Button(self, text="Previous Page",command=lambda:controller.show_frame("SummaryPage"))
        global_vars.PreviousButton.pack(side=LEFT,pady=0,padx=0)
        #Button to start the filtering process. Saving button this way allows us to control its configuration. 
        global_vars.startbutton=Button(self,text="START",command = lambda: _thread.start_new_thread(self.analysis_steps, ()))
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

    #This is the definition for running the analysis steps. Basically, every process that will be performed needs to be included here. 
    # For instance, filtering reads by mapquality neesd to be here, and following that filtering by mitochondrial reads etc. 
    # Each of the process that will be performed will be a new definition that willl be called under this definition!


    def analysis_steps(self):

        #Define the total number of steps for progress bar
        global_vars.total_steps = 8 + len(global_vars.bamfileLocation)

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
        write_to_log("\n\n STEP 9: Removing Unnecessary Files")
        RemoveFiles = self.RemoveFiles(bamfilelist = global_vars.bamfileLocation)
        if RemoveFiles is False:
            write_to_log("Error. Could not remove unnecessary filesfiles")
            messagebox.showerror(title="Error!",message="Could not remvoe unnecessary files")
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
            write_to_log("Filtering for map quality...")
            
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
            write_to_log("Filtering reads for mitochondrial DNA...")

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
            write_to_log("Now Removing duplicate reads...")
            #write_to_log("The current working directory is %s" % outputdir)

            #This is to remove the spaces from the file list that was saved under "gobal_vars.MapQualityFileLocation" in previous def
            #Also write to log
            separator=''
            readfile = separator.join(i)

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
                
                cmd1 = "picard MarkDuplicates REMOVE_DUPLICATES=true I=%s" % readfile
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
                messagebox.showerror(title="Error",message="Could not remove duplicates from the bam file")
                write_to_log("Error! Could not remove duplicates from the bam file")
                return False

        global_vars.DupsRemoved =  Outputfilesdir
        return True


    
    def sortchMreads(self,bamfilelist):
        outputdir= os.path.join(global_vars.EXPERIMENT_SETTINGS['Experiment directory'] + "/")
        os.chdir(outputdir)

        Outputfilesdir = []
        for i in bamfilelist:
            write_to_log("Sorting reads...")
            #write_to_log("The current working directory is %s" % outputdir)

            #This is to remove the spaces from the file list that was saved under "gobal_vars.MapQualityFileLocation" in previous def
            #Also write to log
            separator=''
            readfile = separator.join(i)

            #Define output file names
            outfile = "File" + str(bamfilelist.index(i)+1) + global_vars.File_Flags['sortchrM_Reads'] + ".bam"
           #write_to_log("The file being processed is %s. Processed file will be saved as %s" % (readfile,outfile))
            
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
            write_to_log("Indexing BAM file for viewing in IGV (This generates a .bai file)")
            #write_to_log("The current working directory is %s" % outputdir)

            #This is to remove the spaces from the file list that was saved under "gobal_vars.MapQualityFileLocation" in previous def
            #Also write to log
            separator=''
            readfile = separator.join(i)

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
            write_to_log("Generate fragment size distribution graph using picards...")
            #write_to_log("The current working directory is %s" % outputdir)

            #This is to remove the spaces from the file list that was saved under "gobal_vars.MapQualityFileLocation" in previous def
            #Also write to log
            separator=''
            readfile = separator.join(i)
            #write_to_log("The bam file being processed is %s" % readfile)

            #Define output file names
            outfile_pdf = "File" + str(bamfilelist.index(i)+1) + global_vars.File_Flags['fragment size distribution'] + ".pdf"
            outfile_txt = "File" + str(bamfilelist.index(i)+1) + global_vars.File_Flags['fragment size distribution'] + ".txt"
            
            #readfile = i
            fileaccess = check_file_access(readfile) #checking file access

            if fileaccess is True:
                #Building command pipeline:
                #write_to_log("Initiating picard suite")

                INPUT   = "picard CollectInsertSizeMetrics I=%s" % readfile 
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

        #list of broadpeak files
        BroadPeakFileName = []

        for i in bamfilelist:
            write_to_log("Calling peaks using MACS2...")
            #write_to_log("The current working directory is %s" % outputdir)

            #This is to remove the spaces from the file list that was saved under "gobal_vars.MapQualityFileLocation" in previous def
            #Also write to log
            separator=''
            readfile = separator.join(i)

            #Define output file names
            outfile = "File" + str(bamfilelist.index(i)+1) + global_vars.File_Flags['Peak_Calls'] 
            #write_to_log("The file being processed is %s. Processed file will be saved as %s" % (readfile,outfile))

            #Broadpeak files append:
            peaks_dir = "PEAKS_%s" % outfile
            broadpeaks_dir = os.path.join(global_vars.EXPERIMENT_SETTINGS['Experiment directory'] + "/" + peaks_dir + "/")
            broadpeak_file = str(broadpeaks_dir) + "splitted_peaks.broadPeak"
            BroadPeakFileName.append(broadpeak_file)

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
        global_vars.broadpeak_files = BroadPeakFileName

        #write_to_log("The broadpeak file is %s" % global_vars.broadpeak_files)
        #write_to_log("The output directory list is as follows: %s" % global_vars.peak_folder)

        return True



    def broadpeak_bed_conversion(self, bamfilelist):
        outputdir= os.path.join(global_vars.EXPERIMENT_SETTINGS['Experiment directory'] + "/")
        os.chdir(outputdir)

        for i in bamfilelist:
            write_to_log("Generating BED file from broadpeak file...")
            #write_to_log("The current working directory is %s" % outputdir)

            #This is to remove the spaces from the file list that was saved under "gobal_vars.MapQualityFileLocation" in previous def
            #Also write to log
            separator=''
            readfile = separator.join(i)

            #Define output file names
            outfile_name = "File" + str(bamfilelist.index(i)+1) + global_vars.File_Flags['Bed File'] + ".bed"
            #write_to_log("The file being processed is %s. Processed file will be saved as %s" % (readfile,outfile_name))
            
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

            cmd = "rm mm10.fa.fai %s %s %s %s %s" % (MapQuality,chrM,sort_chrM,dupsremoved_txt,fragSize_txt)
            
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

        ###Select experimental bam files###
        #bamfiles= global_vars.bamfileLocation
        #variable=StringVar(main)
        #control1=tk.Frame(main,relief=controller.frame_relief,borderwidth=controller.frame_borderwidth,bg="#8c8c8c")
        #control1.pack(fill=BOTH,pady=0,padx=0,anchor=N,expand=FALSE)
        #self.dropdowncontrol=OptionMenu(control1,variable,*bamfiles)
        #self.dropdowncontrol.pack(fill=BOTH,pady=0,padx=0,side=LEFT,expand=TRUE)
        
        
       
        

#Close the mainloop to wrap the software:
if __name__ == "__main__":
    app = ProcessingTool()
    app.geometry("740x480+300+50")
    app.resizable(True, True)
    app.mainloop()
    _thread.start_new_thread(app.mainloop(), ())
