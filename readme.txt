#####################################################################################
#    Filename : readme.txt
#    Date : Aug 25, 2015
#    What : The readme file for Proj_Resonance_Lorentzian_GitHub_Run_Me.ipynb
#           and def_resonance_lorentzian_for_github.py
#          
#
#####################################################################################

What is here :
================================

"Proj_Resonance_Lorentzian_GitHub_Run_Me.ipynb" 
+
"def_resonance_lorentzian_for_github.py"

=

Together, they will read a file containing a list of numbers, and they will tell if that list of numbers have any resonances in it. Additionally, it will then fit those minimas to a lorentzian curve, and it will tell you what the chi-squared of all those individual minimas are. 

“2015-07-01-starcryo-nevins-initial-sweep.npz” : the file I tested on

“All_the_plots.zip” : all the lorentzian and resonance plots of the test


What you need to put into this :
================================
	 
A) A file, and you would have to specify where it is. As an example, I have left in the  programme where I would have used to locate the file on my computer. You can just cut and paste to get to where the file is on your computer. An example file ".npz" is attached here.
   
B) You need to tell the programme what the frequency and data number are. 
This is usually shortened to 'freq' and 'data' respectively, without the parenthesis.
   
C) Finally, you will get 3 options on how you may want to look at the Lorentzian
.1) This is just all the resonances close up
.2) This is just all the resonances close up and they would be fitted with the Lorentzian (pick me!)
.3) This is just all the resonances close up, followed by all the resonances close up that are fitted with the Lorentzian 
       

A word of caution:
===================
I wouldn't change anything in this programme except these 4 parameters when looking for resonance:
''' Just toggle these 4 parameters'''
flat = 1.2            # flat is the gradient of the points around centre // flat = 1.2
move = 5              # move is the number of points around of centre    // move = 5
noise = 5             # noise is how many points away from 'move'        // noise = 5
tol = 1               # tolerance determine spread of lorentz function   // tol = 1






