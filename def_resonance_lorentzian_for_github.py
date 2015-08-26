#####################################################################################
#    Filename : def_resonance_lorentzian.py
#    Date : Aug 15, 2015
#    What : Main list of definitions for the main file 
#           "Proj_Resonance_Lorentzian_GitHub_Run_Me.ipynb"
#####################################################################################

'''Definitions by AstroTze'''
import numpy as np
import matplotlib.pyplot as plt


'''This gives you two maximas (one to left, one to right) away from minima(Champ)'''
def give_me_2_maximas(champ,y,move,noise,tol):
    away = move + noise
    temp_champ = champ                           
    hl = 0                             
    hr = 0                             
    while y[temp_champ ] <= y[temp_champ - tol*away] and temp_champ !=0:            
        temp_champ = temp_champ - 1                        
        hl = hl +1                   
    temp_champ = champ                          
    while y[temp_champ] <= y[temp_champ + tol*away] and temp_champ != len(y)-tol*away:            
        temp_champ = temp_champ + 1                        
        hr = hr +1
    return hl,hr

'''Gives you two points above 2 Full Width Half Max, counting from Minima to Maxima'''
def give_me_2_half_maxes(champ,y,h_fwhm): 
    j = champ                    
    xl = 0                       
    xr = 0                       
    while y[j] < h_fwhm:         
        j = j-1                  
        xl = xl + 1               
    j = champ                    
    while y[j] < h_fwhm:         
        j = j+1                  
        xr = xr + 1               
    return xl,xr

'''This gives you the all important gamma'''
def give_me_gamma(x,champ,xl,xr):
    gam1 = abs(abs(x[champ - xl]) - abs(x[champ]))
    gam2 = abs(abs(x[champ + xr]) - abs(x[champ]))
    gamma=(gam1 + gam2)/2
    return gamma

'''Returns the difference between an ideal, and a measured y.Used for chi-square'''   
def residuals(p,y_meas, x_ideal):  
    err = y_meas - neg_Lorentz(x_ideal,p)        # y_meas - y_ideal
    return err 

'''This returns the Lorentzian function'''
def neg_Lorentz(x, p):                           #returns the negative lorentzian 
    Numerator = p[2]                             #p[0] = Considered as 'x0': centre of x  
    Denominator = ((x - p[0])**2 + p[2])         #p[1] = (gamma*np.pi) 
    Co_eff = 1/p[1]                              #p[2] = gamma**2 
    Background = p[3]                            #p[3] =  Considered as 'y0': background
    return (-1 * Co_eff * (Numerator/Denominator)) + Background
 
'''This code returns will return the average gradient to the right side of the curve'''   
def right_gradient_average(x,y,center,move,away):
    dx = np. gradient(x[center + move : center + away])
    right_gradient = np.gradient(y[center + move: center + away],dx)
    right_gradient_ave = np.average(right_gradient)
    return right_gradient_ave 

'''This code returns will return the average gradient to the left side of the curve''' 
def left_gradient_average(x,y,center,move,away):
    dx = np. gradient(x[center - away : center - move])
    left_gradient = np.gradient(y[center - away : center - move],dx)
    left_gradient_ave = np.average(left_gradient)
    return left_gradient_ave
    
'''This returns ALL the local minimums'''
def minima(x, y, move, noise,flat):                # move points before and after minimum
    locmins = []
    away = move + noise
    for center in range(away, len(y) - (away)):    # noiseth point after move is checked
        left_min = 0
        right_min = 0 
        for j in range(move):                                                 
            if y[center-j] <= -1* abs(y[center - (j+noise)] + y[center - (j+(noise-1))]) / 2 and y[center] <= y[center-1]:   
                left_min = left_min + 1   
            if y[center+j] <= -1* abs(y[center + (j+noise)] + y[center + (j+(noise-1))]) / 2 and y[center] <= y[center+1]:
                right_min = right_min + 1 
        
        right_gradient_ave = right_gradient_average(x,y,center,move,away)
        left_gradient_ave = left_gradient_average(x,y,center,move,away)
        
        if left_min + right_min == 2*move:
            if left_gradient_ave < 0 and right_gradient_ave > 0:
                if left_gradient_ave < -1*flat and right_gradient_ave > flat:
                    locmins.append([ x[center], y[center] ])
    return locmins                                  # Returns set of all local minimas

'''This removes any local minimums that are the same and next to each other '''
def remove_repeated(locmins):
    for i in range(1,len(locmins)):
        if locmins[i][1] == locmins[i-1][1]:        # The second index [1] is the y-axis!
            locmins[i][0]=0
            locmins[i][1]=0 
    while [0,0] in locmins: locmins.remove([0,0])
    return locmins

'''This code ranks the local minimums for you, according to indexes in y '''
def rank_the_locmins(locmins,y):
    for i in range(len(y)):                         # Sorts locmins into their position
        for j in range (len(locmins)):
            if y[i] == locmins[j][1]:
                if len(locmins[j])==2:              # Prevents picking up repeated points
                    locmins[j].append(i)
    return locmins

'''This will sort out the locmins and remove any minimas which are too close'''
def remove_close_by_mins(locmins,move,noise):
    for i in range(1,len( locmins)):
        after = locmins[i][2]
        before = locmins[i-1][2]
        diff = abs(after - before)
        away = move + noise
        if diff < away :
            y_after = locmins[i]
            y_before = locmins[i-1]
            if y_after < y_before :
                locmins[i-1] = [0,0,0]
            if y_after > y_before :
                locmins[i] = [0,0,0]
    while [0,0,0] in locmins: locmins.remove([0,0,0])
    return locmins

'''Let us find the range of frequencies to look at : if the freq are too closed, the range is the same '''
def ranges(locmins):
    range_of_freq = []
    for i in range(len(locmins)):
        range_of_freq.append(locmins[i][2])
    for i in range(1,len(range_of_freq)):
        if range_of_freq[i] - range_of_freq[i-1] < 300:
            range_of_freq[i-1] = 0
    while 0 in range_of_freq : range_of_freq.remove(0)
    return range_of_freq 



'''This code just allows me to look at ranges of values for the lorentzian '''
def ranges_to_look(ranges):
    list_of_differences = []
    for i in range(1, len(ranges)):
        difference = ranges[i] - ranges[i - 1]
        list_of_differences.append(difference)
    #print list_of_differences
    ave_diff = int(np.mean(list_of_differences))
    #print ave_diff  
    list_of_ranges = []
    for i in range(len(ranges)):
        if i == 0:
            begin = ranges[i] - ave_diff / 2
            end = ranges[i] + (ranges[i + 1]- ranges[i]) / 2
            this_range = [begin, end]
            list_of_ranges.append(this_range)
        if i == len(ranges) - 1:
            begin = ranges[i] - (ranges[i] - ranges[i - 1]) / 2
            end = ranges[i] + ave_diff / 2
            this_range = [begin, end]
            list_of_ranges.append(this_range)
        elif i != 0 and i != len(ranges) -1 :
            begin =  ranges[i] - ( ranges[i] - ranges[i - 1] )/2
            end = ranges[i] + (ranges[i + 1]- ranges[i]) / 2
            this_range = [begin, end]
            list_of_ranges.append(this_range)
    return list_of_ranges

'''This code should tell you all the resonators, according to the Lorentzian definition'''
def resonators(flat,move,noise,tol,x,y):
    order_of_mins = []
    locmins = minima(x, y, move, noise,flat)                 
    remove_repeated(locmins)                              
    locmins = rank_the_locmins(locmins,y)   
    locmins = remove_close_by_mins(locmins,move,noise)
    for i in range(len(locmins)):
        order_of_mins.append([i+1, locmins[i][0], locmins [i][1]])
    return locmins, order_of_mins

'''This code will convert any rankings into frequency'''
def convert_rank_to_freq(ranges_to_look,x):
    freq_to_look= []
    for i in range(len(ranges_to_look)):
        j_begin = ranges_to_look[i][0]
        j_end = ranges_to_look[i][1]
        freq_begin = x[j_begin]
        freq_end = x[j_end]           
        freq_index = [freq_begin, freq_end]
        freq_to_look.append(freq_index)
    return freq_to_look

'''This should the range of values to the min and max of the data'''
def convert_rank_to_min_and_max_data(ranges_to_look,y):
    data_to_look = []
    for i in range(len(ranges_to_look)):
        begin = ranges_to_look[i][0]
        end = ranges_to_look[i][1]
        data_max = max(y[begin:end])
        data_min = min(y[begin:end])
        data_interval = (data_max - data_min)/10
        data_max = data_max + data_interval
        data_min = data_min - data_interval
        data_to_look.append ([data_min,data_max])
    return data_to_look

'''This should give you all the limiting ranges to look at'''
def give_me_all_limits(ranges_to_look,freq_to_look,data_to_look,j):
    ranges_begin,ranges_end = ranges_to_look[j][0], ranges_to_look[j][1]
    xbegin, xend = freq_to_look[j][0] ,freq_to_look[j][1]
    ybegin, yend = data_to_look[j][0], data_to_look[j][1]
    return ranges_begin, ranges_end, xbegin, xend, ybegin, yend

'''This lets you count the number of data points in your range'''
def count_no_of_data_points_per_close_range(order_of_mins,x,ranges_begin,ranges_end,freq_points,index,index_list,i):
    if order_of_mins[i][1] in x[ranges_begin: ranges_end]:
        freq_points.append(order_of_mins[i][1])
        index = index + 1
        index_list.append(index)
    return freq_points, index, index_list

'''This code should give you 2 decimal places'''
def give_me_two_decimal(freq_points):
    freq_points_dec = []
    for j in range(len(freq_points)):
        freq_to_2_dec = "{0:.2f}".format(freq_points[j])
        freq_points_dec.append(freq_to_2_dec)  
    return freq_points_dec

'''This tells you where to look for the resonances'''
def where_do_i_look(locmins,x,y):
    ranging = ranges(locmins)
    where_to_look = ranges_to_look(ranging)
    freq_to_look = convert_rank_to_freq(where_to_look,x)
    data_to_look = convert_rank_to_min_and_max_data(where_to_look,y)
    return where_to_look, freq_to_look, data_to_look


'''This tells you initially where to look, and how to box the overall plot in'''
def box_in_the_plot(ranges_to_look, freq_to_look, data_to_look):
    ystart = 0
    yend = 0
    for j in range(len(ranges_to_look)):
        xstart,xend = freq_to_look[0][0],freq_to_look[j][1]
        ystart = ystart + data_to_look[j][0]
        yend = yend + data_to_look[j][1]
    ystart = ystart / len(data_to_look)
    yend = yend / len(data_to_look)

    x_interval =abs( xend - xstart )/ len(freq_to_look)
    y_interval =abs( yend - ystart )

    xstart = xstart - 2*x_interval
    xend = xend + 2*x_interval
    ystart = ystart - y_interval
    yend = yend + y_interval
    return xstart, xend, ystart, yend

'''This shows you the main plot'''
def main_plot_with_noise(x,y,xstart,xend,ystart,yend):
    plt.plot(x,y)
    plt.title('Main plot with noise')
    plt.xlim(xstart,xend )
    plt.ylim(ystart,yend)  
    
'''This shows you where the resonators are'''
def show_me_resonators(order_of_mins,data,x,y,xstart,xend,ystart,yend):
    print 'length of data',len(data)
    print "There are",len(order_of_mins),"resonators, and they occur at ...'"
    print "Resonator","\t", "Frequency", "\t", "data numbers"
    for i in range(len(order_of_mins)):
        print order_of_mins[i][0],'\t','\t', order_of_mins[i][1],'\t', order_of_mins[i][2]
    plt.plot(x, y)
    plt.title ( 'Plot with the Resonators')   
    plt.xlim(xstart,xend )
    plt.ylim(ystart,yend)                            
    ax = plt.gca()                                              
    for i in range(len(order_of_mins)):
        ax.axvline(order_of_mins[i][1], color = 'red',linewidth = 2, alpha = 0.7)         
        plt.scatter (order_of_mins[i][1], order_of_mins[i][2], color = 'red', s = 40)

        
        
'''This should show you all the resonators in close range'''                   
def show_me_resonators_in_close_range(ranges_to_look,freq_to_look,data_to_look,x,y,order_of_mins):
    index = 0
    index_list=[]
    for j in range(len(ranges_to_look)):
        ranges_begin, ranges_end, xbegin, xend, ybegin, yend \
        = give_me_all_limits(ranges_to_look, freq_to_look, data_to_look,j)
        plt.plot(x,y)
        plt.scatter(x,y, alpha = 0.4)
        plt.xlim (xbegin, xend)
        plt.ylim (ybegin, yend)
        freq_points = []
        ax = plt.gca()
        
        for i in range(len(order_of_mins)):
            ax.axvline(order_of_mins[i][1], color = 'red')         
            plt.scatter (order_of_mins[i][1], order_of_mins[i][2], color = 'red')  
            freq_points, index, index_list = count_no_of_data_points_per_close_range\
            (order_of_mins,x,ranges_begin,ranges_end,freq_points,index,index_list,i)
            
        freq_points_dec = give_me_two_decimal(freq_points)    
        title = " The data point is", index_list,"...and the frequency is", freq_points_dec
        plt.title(title, fontsize = 15)
        plt.xlabel('Frequency') 
        plt.ylabel('Data')
        index_list = []
        plt.figure()        
        
        
        
'''This fits for you the lorentzian, and will tell you which points have chi-squared more 
than 1, and which ones have more than 100 '''
def fit_me_to_lorentzian(champ,x,y,move,noise,tol,flat,count_chi,count_chi_less_than_1):
            
    hl,hr =  give_me_2_maximas(champ,y,move,noise,tol)# hl, hr : height to left, right              
    fwfm = abs (x[champ + hr] - x[champ - hl])        # fwfm   : full width full max
    h_fwfm = abs (y[champ + hr] - y[champ])           # h_fwfm : height fwfm    
    h_fwhm = (abs(y[champ]) - 0.5 * h_fwfm) * -1      # fwhm   : full width half max             
    xl, xr = give_me_2_half_maxes(champ,y,h_fwhm)     # h_fwhm : height fwhm        
    gamma = give_me_gamma(x,champ,xl,xr)              # xl,xr : half max to left,right        
                                                                  
    p=[0.0, 0.0, 0.0, 0.0]                            
    p[0] = x[champ]                                             
    p[1] = (gamma * np.pi)
    p[2] = gamma**2
    p[3] = y[champ + hr]
        #print j+1,p, "The original parameters"  
    
    x_ideal = np.linspace(x[champ - hl], x[champ + hr], hl+hr)# x_fit = x_ideal = x_meas
    y_ideal = neg_Lorentz(x_ideal,p)
    y_meas = y[champ - hl : champ + hr]
    x_meas = x_ideal                                      
    
    from scipy.optimize import leastsq
    plsq = leastsq(residuals, p, args=(y_meas,x_meas))      
        #print j+1,(plsq[0]), "The parameters for leastsq"   
    x_fit = x_ideal                                                 
    y_fit = neg_Lorentz(x_ideal,plsq[0])
        
    chi_squared=0
    for item in range(len(y_meas)):
        element = (y_meas[item]-y_fit[item])**2
        chi_squared = chi_squared + element
    if chi_squared > 100:
        count_chi = count_chi +1
    if chi_squared < 1:
        count_chi_less_than_1=count_chi_less_than_1 +1
            
    return chi_squared, count_chi, count_chi_less_than_1, hl, hr, champ, x_ideal, y_ideal, x_fit, y_fit
        
'''This should show you all the resonators in close range fitted to the lorentzian function'''
def show_me_resonators_in_close_range_with_lorentzian(ranges_to_look,freq_to_look,data_to_look,x,y,\
                                                      move,noise,tol,flat,order_of_mins,locmins):
    index = 0
    index_list=[]
    count_chi=0
    count_chi_less_than_1=0
    chi_squared_list =[]
    chi_squared_total = []
    freq_points_total = []
    for j in range(len(ranges_to_look)):
        ranges_begin, ranges_end, xbegin, xend, ybegin, yend \
        = give_me_all_limits(ranges_to_look, freq_to_look, data_to_look,j) 
        plt.plot(x,y)
        plt.scatter(x,y, alpha = 0.4, label = "Measured")
        plt.xlim (xbegin, xend)
        plt.ylim (ybegin, yend)
        freq_points = []
        ax = plt.gca()
        for i in range(len(order_of_mins)):
            ax.axvline(order_of_mins[i][1], color = 'red')         
            plt.scatter(order_of_mins[i][1], order_of_mins[i][2], color = 'red')  
            freq_points, index, index_list = count_no_of_data_points_per_close_range\
            (order_of_mins,x,ranges_begin,ranges_end,freq_points,index,index_list,i)        
            
        '''This plots the lorentzian per minima'''
        for i in range(len(freq_points)):    
            freq_points_total.append(freq_points[i])
        locmins_list=[]
        chi_squared_list =[]
        for i in range(len(index_list)):
            locmins_list.append(locmins[index_list[i]-1])   
        for j in range(len(locmins_list)):
            champ = locmins_list[j][2]
            chi_squared, count_chi, count_chi_less_than_1, hl, hr, champ, x_ideal, y_ideal, x_fit, y_fit \
            = fit_me_to_lorentzian(champ,x,y,move,noise,tol,flat,count_chi,count_chi_less_than_1) 
            chi_squared_list.append(chi_squared)
            chi_squared_total.append(chi_squared)
            
            plt.scatter(x[champ - hl], y[champ - hl], color = 'orange', s=180 , alpha = 0.8 )
            plt.scatter(x[champ + hr], y[champ + hr], color = 'orange', s=180 , alpha = 0.8)
            plt.plot(x_ideal, y_ideal, color = 'orange', linewidth = 4, linestyle = '--',alpha = 0.8)  
            plt.scatter(x[champ - hl], y[champ - hl], color = 'orange', s=180 , alpha = 0.8 )
            plt.scatter(x[champ + hr], y[champ + hr], color = 'orange', s=180 , alpha = 0.8)
            plt.plot(x_ideal, y_ideal, color = 'orange', linewidth = 4, linestyle = '--',alpha = 0.8) 
            plt.plot(x_fit, y_fit, color = 'green', linewidth = 4 , alpha = 0.5 )
            plt.scatter(x_fit, y_fit, color = 'green', alpha = 0.5 ) 
        #print chi_squared_list
               
#########################################################################################################              
        '''Everything below is just for labelling !!!'''
        if len(locmins_list)!=0:
            plt.scatter(x[champ - hl], y[champ - hl], color = 'orange', s=180 , alpha = 0.8, label = "Maxima" )    
            plt.scatter (locmins[j][0], locmins[j][1], color = 'red', label = "Minima")    
            plt.plot(x_ideal, y_ideal, color = 'orange', linewidth = 4, linestyle = '--', alpha = 0.8 , label = " Ideal" )
            plt.plot(x_fit, y_fit, color = 'green', linewidth = 4 , alpha = 0.5 , label = "Fit" )
            plt.legend(loc=3)
            
        freq_points_dec = give_me_two_decimal(freq_points) 
        chi_squared_dec = give_me_two_decimal(chi_squared_list) 
        title = " Resonator is", index_list,"...and Chi-square is", chi_squared_dec,"...frequency is",freq_points_dec
        plt.title(title, fontsize = 15)
        plt.xlabel('Frequency') 
        plt.ylabel('Data')
        index_list = []
        plt.figure()
        
    return count_chi, count_chi_less_than_1, chi_squared_total,freq_points_total, len(chi_squared_total) 

''' This will tell you all the chi-squared that appears'''
def print_me_those_chi_square(count_chi,count_chi_less_than_1,chi_squared_total,freq_points_total):
    print "The number of fits with Chi_square more than 100 is", count_chi 
    print "The number of fits with Chi_square less than 1 is", count_chi_less_than_1 
    print "Resonator",'\t','Chi_square','\t','at these Frequencies'
    for i in range(len(chi_squared_total)):
        print i+1, '\t','\t',chi_squared_total[i],'\t', freq_points_total[i]
        

    
