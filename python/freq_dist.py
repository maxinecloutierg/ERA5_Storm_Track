import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

"""

    Maxine Cloutier-Gervais

	Created :

    	June 16th, 2023

	Info :

    	This code creates fig 1.e in Chen & al (2022), which gives the frequency distribution (%)
    	of standardised VORS spatiotemporally averaged around a 800 km radius and over the lifetime 
    	in CRCM6 domain


"""

def frequency(df, bins)

    """
   
	Determines the frequency distribution (%) of standardised VORS spatiotemporally averaged 
    around a 800 km radius and over the lifetime in CRCM6 domain

    Parameters :
		df     : name of the dataframe
		bins   : bin interval (numpy array)

    Return     : 
		pct    : Frequency distribution (%) in numpy array form

    """
    
    # Step 1 : Determine average of VORS_av08 over all grid points
    #          that are within CRCM6 domain

    avg_vors = df.groupby('storm')['VORS_av08'].mean()

    # Step 2 : Determine standardized value of the averaged VORS_av08  

    std_vors = (avg_vors - avg_vors.mean()) / avg_vors.std()

    # Step 3 : Calculate occurence of every value for every bin  

    frequency, _ = np.histogram(std_vors, bins=bins)

    # Step 4 Calculate frequency distribution (%)
    
    tot_etc = len(std_vors)
    pct = (frequency / tot_etc) * 100

    return pct


def custom_legend(seasons, colors) : 
	
    """
    
    Creates customized legend

    Parameters  : 
		seasons : list of seasons to plot
		colors  : list of colors for each season

    """
    # Step 1 : Iterate through all seasons and customize accordingly

    for i, color in zip(range(0,len(seasons)), colors):

    legend.get_texts()[i].set_position((-9, 0.5))
    first_text = legend.get_texts()[i]  # Get the first legend text
    first_text.set_fontweight('bold')   # Set font weight to bold
    first_text.set_color(color)
    first_text.set_alpha(0.8)           # transparency


def custom_layout(max_pct, bins) : 

    """
    
    Creates the layout for the barplot

    Parameters  : 
		max_pct : max value of frequency distribution
		bins    : bin values (numpy array)

    """

    # Step 1 : set ticks 
    
	plt.xticks(bins[:-1])
    plt.xticks(fontsize=14)
    plt.yticks(fontsize=14)
    plt.yticks(np.arange(0, max_pct, 3))

    # Step 2 : set labels 
    
    plt.xlabel('Standardized NNA-lifetime and 800 km averaged VORS', fontsize=14)
    plt.ylabel('Frequency (%)', fontsize=14)

    # Step 3 : set grid

    plt.grid(which='major', axis='x', linestyle='dotted')




def set_textbox(mean, sigma)

    """
   
    Create text box containing average VORS_av08 and standard deviation

    Parameters  : 
		mean    : mean VORS_av08
		sigma   : VORS_av08 standard deviation

    Return      : 
		teststr : formatted string ready to be added to the bar plot

    """

    # Step 1 : Change mean and sigma format to keep 2 significant numbers

    mean = '{:0.2e}'.format(mean)
	sigma = '{:0.2e}'.format(sigma)
	
	# Step 2 : Transform first 4 digits into string (ex : 0.93e-05 = "0.93")
	mean = mean[0 : 4]
    sigma = sigma[0 : 4]

    # Step 3 : Create textbox
    
    textstr = '\n'.join((
    r'$\mathrm{mean}:\ %s \times 10^{-5} \ \mathrm{s^{-1}}$' % mean,
    r'$\qquad \ \sigma:\ %s \times 10^{-5} \ \mathrm{s^{-1}}$' % sigma))
    
    return textsrt




def create_fig(file_in , file_out, seasons)

    """
    
    Create the frequency distribution figure and saves it in png format

    Parameters   : 
		file_in  : path of the csv file to read
		file_out : path to save the png image
	seasons  : seasons to be plotted

    """
	# Step 1 : read csv file and affect initial variables

    df24 = pd.read_csv(file_in)
    bins = np.arange(-2, 5, 0.5)
 	
	crcm6 = df24.loc[df24.HU == True]             # keep rows that are within CRCM6 domain
	colors = ['red', 'orange', 'royalblue', 'g']

	i = 0.1                                       # To shift bars position in plot
	max_pct = 0                                   # To format y axis in plot
	
	# Step 2 : Iterate through all seasons and create bar plot
	
	for season, color in zip(seasons, colors) : 
		
		df = df24.loc[(df24.season == season) & (df24.HU == True)]
		pct = frequency(df, bins)
		plt.bar(bins[:-1]+i, pct, width=0.15, alpha=0.7,label=season, color=color)

		# establish the max frequency for y limit axis on the bar plot later

		if max(pct) > max_pct : 
			max_pct = max(pct)
		
		i += 0.1

	# Step 3 : Add custom legend
	
	legend = plt.legend(frameon=False, ncol=len(seasons), fontsize=16, columnspacing=0.5)
	custom_legend(seasons, colors)
	
	# Step 4 : Add x and y axis
	
	custom_layout(max_pct, bins)

	# Step 5 : Calculate mean and standard deviation
	
	mean = vors['VORS_av08'].mean()
	sigma = vors['VORS_av08'].std()

	# Step 6 : Add textbox

	textstr = textstr(mean, sigma)
	plt.text(2.2, 8, textstr, fontsize=16)

	# Step 7 : show plot and save output in file_out
	
	plt.show()
	plt.savefig('/pampa/cloutier/freq_dist.png')



