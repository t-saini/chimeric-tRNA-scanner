#!/usr/bin/env python3
import matplotlib.pyplot as plt
import sys

'''
The purpose of Chimeric tRNAscanner (chi-tRNA) is to visually dispaly tRNA oddities in species.

This program takes in data from tRNAscan-SE as a text document. It parses through
each entry and plots the tRNA type against the predicted Isotype. This is done
by converting the data into numeric vaules which are then plotted on a graph.

Users also have the option have changing the x-axis between the

A "nomral line" is displayed on the graph. This graph is the expected relativity
of the Isotype and tRNA type. For every amino acid encoded by the anticodon on
the X-axis it should line up with the Isotype on the Y-axis.

This program is designed through command line, however if the program is ran
without any parameters, a simple prompt will be executed asking the user for
the following values:

tRNAscan-SE file name
Infernal score cut off
X - axis choice of amino acid or anticodon
Output for supplemental data if desired

From command line
input: chi-tRNA.py tRNAscan-SE 2.0 data file name

ouput: Generate scatter plot showing tRNA types (X) vs isoforms(Y)

For advanced options within command line:
chi-tRNA.py -h

Author: Tanvir Singh Saini
Date: June 5, 2016
'''

class CommandLine() :
    '''
    Handle the command line, usage and help requests.

    CommandLine uses argparse, now standard in 2.7 and beyond.
    it implements a standard command line argument parser with various argument options,
    a standard usage and help, and an error termination mechanism do-usage_and_die.

    attributes:
    all arguments received from the commandline using .add_argument will be
    avalable within the .args attribute of object instantiated from CommandLine.
    For example, if myCommandLine is an object of the class, and requiredbool was
    set as an option using add_argument, then myCommandLine.args.requiredbool will
    name that option.

    '''

    def __init__(self, inOpts=None) :
        '''
        CommandLine constructor.
        Implements a parser to interpret the command line argv string using argparse.
        '''

        import argparse
        self.parser = argparse.ArgumentParser(description = 'This program was designed to parse through tRNAscan-SE\'s output file.'
                                                'This will generate a graph with anticodons on the X-axis with \"tRNA Graph\" as the title of the graph. The data points presented will be greater than or equal to an Infernal score of 90.',
                                             epilog = 'This chimeric tRNA-scanner (chi-tRNA) was designed by Tanvir Singh Saini for Dr.Todd Lowe of UCSC',
                                             add_help = True, #default is True
                                             prefix_chars = '-',
                                             usage = '%(prog)s file.txt -i 90 -o anticodon -s N -t "tRNA Graph" -m 12'
                                             )
        self.parser.add_argument('inFile', action = 'store', help='input file name')
        self.parser.add_argument('-i', '--infscore', type = float ,action = 'store', help='infernal Score cut off for tRNA expression')
        self.parser.add_argument('-o', '--option', action = 'store', help='use aa(amino acid) or anticodon for the X-axis')
        self.parser.add_argument('-s', '--supplemental' , action = 'store', help = 'select y/N to output a text file with chimeric tRNAs')
        self.parser.add_argument('-t', '--title', action = 'store', help = 'title graph; example input "Yeast tRNAs"')
        self.parser.add_argument('-m', '--mod', type = int, action = 'store', help='modify the column entry for the isoType')
        self.parser.add_argument('-f', '--filter', action = 'store', help = 'change score filter from Inf score to Isotype score' )
        self.parser.add_argument('-v', '--version', action='version', version='%(prog)s 1.0')
        if inOpts is None :
            self.args = self.parser.parse_args()
        else :
            self.args = self.parser.parse_args(inOpts)

class Module_Scanner:
    '''
    From command line
    input: ChitRNA-scan.py tRNAscan-SE 2.0 data file name

    ouput: Generate scatter plot showing tRNA types (X) vs isoforms(Y)

    For advanced options within command line:
    ChitRNA-scan.py -h

    Author: Tanvir Singh Saini
    Date: June 5, 2016
    '''

    aa_numeric = {
    'Ala': 1, 'Arg': 2, 'Asn' : 3,'Asp' : 4, 'Cys' : 5,
    'Gln': 6, 'Glu': 7, 'Gly' : 8,'His' : 9, 'Ile' : 10,
    'Leu': 11,'Lys': 12,'Met': 13,'Phe': 14, 'Pro': 15,
    'Ser': 16,'Thr': 17,'Trp': 18,'Tyr': 19, 'Val': 20,
    'SeC': 21, 'Pyl': 22,'iMet':13,'fMet': 13, 'Ile2':23, 'Undet': 24,
    'Sup': 25
    }

    AntiNumeric = {
     #Ala 1-4, Arg 5-10, Asn 11-12, Asp 13-14, Cys 15-16, Gln 17-18, Glu 19-20, Gly 21-24
     #His 25-26, Ile 27-29, leu 30-35, lys 36-37, Met 38, phe, 39-40, Pro 41-44, Ser 45-50
     #Thr 51- 54, Trp 55, Tyr 56-57, Val 58-61, SeC 62, Pyl 63, iMet 38, fMet 38, Stop 64
    'AAA': 39, 'AGA': 45, 'ATA': 56, 'ACA': 15,
    'GAA': 40, 'GGA': 46, 'GTA': 57, 'GCA': 16,
    'TAA': 30, 'TGA': 47, 'TTA': 64, 'TCA': 62,
    'CAA': 31, 'CGA': 48, 'CTA': 63, 'CCA': 55,

    'AAG': 32, 'AGG': 41, 'ATG': 25, 'ACG': 5,
    'GAG': 33, 'GGG': 42, 'GTG': 26, 'GCG': 6,
    'TAG': 34, 'TGG': 43, 'TTG': 17, 'TCG': 7,
    'CAG': 35, 'CGG': 44, 'CTG': 18, 'CCG': 8,

    'AAT': 27, 'AGT': 51, 'ATT': 11, 'ACT': 49,
    'GAT': 28, 'GGT': 52, 'GTT': 12, 'GCT': 50,
    'TAT': 29, 'TGT': 53, 'TTT': 36, 'TCT': 9,
    'CAT': 38, 'CGT': 54, 'CTT': 37, 'CCT': 10,

    'AAC': 58, 'AGC': 1, 'ATC': 13, 'ACC': 21,
    'GAC': 59, 'GGC': 2, 'GTC': 14, 'GCC': 22,
    'TAC': 60, 'TGC': 3, 'TTC': 19, 'TCC': 23,
    'CAC': 61, 'CGC': 4, 'CTC': 20, 'CCC': 24,
    'NNN': 65
    }

    def __init__(self, fname, user_value, option, title, mod, sup, filt):
        '''Initializes objects'''
        self.fname = fname
        self.user_value = user_value
        self.option = option
        self.title = title
        self.mod = mod
        self.sup = sup
        self.filt = filt

    def parser(self):
        '''Take tRNAscan-SE data and append each line to list data'''
        data = []
        with open(self.fname) as f:
            line = f.readline()
            for line in f:
                listit = line.replace('\n','').split('\t') #generates a list at
                #each tab after remove all \n in the tRNAscan-SE output data.
                data.append(listit)#adds lists to data, making data a list of lists.
            data = data[2:] #makes index entry 2 onward the list for data to
            #remove the titles for each column.
            return data

    def tRNA_type(self):
        '''Take list from function parser and extract tRNA type and convert it
        to a numeric value from dictionary aa_numeric'''
        if self.option == 'aa':
            data = self.parser() #calls function parser and stores the output as list data
            transfer = []
            tRNA = [] #stores the tRNA amino acid type
            xtRNA = [] #list that stores the tRNA amino acids as numeric values
            while 0 < len(data): #while the length of list data is greater than 0
                transfer.append(data[0])#take the 0th index and append it to list transfer
                analyze = [entry[4] for entry in transfer]#within in list transfer, look at the list (index 0 that was appended to transfer) and save the 4th index to name analyze
                tRNA.append(analyze)# append analyze to list tRNA
                transfer.pop(0)#pop the 0th index in list transfer to make it empty
                data.pop(0)#pop the 0th index in list data to make index 1 the new index 0
            if 0 == len(data): #once length of list data is zero
                tRNA = [val for sublist in tRNA for val in sublist] #flatten nested list
            for amino_acid in tRNA: #for each amino acid in list tRNA
                xtRNA.append(Module_Scanner.aa_numeric[amino_acid])#find the numeric value for the amino acid in dictionary aa_numeric, and append it to xtRNA
            return xtRNA

    def antiCodon_type(self):
        '''Take list from function parser and extract tRNA type and convert it
            to a numeric value from dictionary aa_numeric'''
        if self.option == 'anticodon':
            data = self.parser() #calls function parser and stores the output as list data
            transfer = []
            antiCodon = [] #stores the Anticodon type
            xtRNA = [] #list that stores the tRNA amino acids as numeric values
            while 0 < len(data): #while the length of list data is greater than 0
                transfer.append(data[0])#take the 0th index and append it to list transfer
                analyze = [entry[5] for entry in transfer]#within in list transfer, look at the list (index 0 that was appended to transfer) and save the 5th index to name analyze
                antiCodon.append(analyze)# append analyze to list tRNA
                transfer.pop(0)#pop the 0th index in list transfer to make it empty
                data.pop(0)#pop the 0th index in list data to make index 1 the new index 0
            if 0 == len(data): #once length of list data is zero
                antiCodon = [val for sublist in antiCodon for val in sublist] #flatten nested list
            for amino_acid in antiCodon: #for each amino acid in list antiCodon
                xtRNA.append(Module_Scanner.AntiNumeric[amino_acid])#find the numeric value for the amino acid in dictionary AntiNumeric, and append it to xtRNA
            return xtRNA

    def isoType(self):
        '''
        Take list from function parser and extract isotype data and convert it
        to a numeric value from dictionary aa_numeric
        '''
        data = self.parser() #calls function parser and stores the output as list data
        transfer = []
        isotype = [] #stores the tRNA isotype
        yisotype = [] #list that stores the tRNA isotype as numeric values
        while 0 < len(data): #while the length of list data is greater than 0
            transfer.append(data[0]) #take the 0th index and append it to list transfer
            analyze = [entry[self.mod] for entry in transfer] #within in list transfer, look at the list (index 0 that was appended to transfer) and save the 4th index to name analyze
            isotype.append(analyze)# append analyze to list isotype
            transfer.pop(0)#pop the 0th index in list transfer to make it empty
            data.pop(0)#pop the 0th index in list data to make index 1 the new index 0
        if 0 == len(data): #once length of list data is zero
            isotype = [val for sublist in isotype for val in sublist]
        for amino_acid in isotype: #for each amino acid in list isotype
            yisotype.append(Module_Scanner.aa_numeric[amino_acid]) #find the numeric value for the amino acid in dictionary aa_numeric, and append it to yisotype
        return yisotype

    def analysis(self):
        '''
        Compare the tRNA type and isotype and output them as a smaller string
        if they are not the same
        '''
        if self.sup == 'y':
            data = self.parser() #take output list from parser function and set it to name data
            transfer = [] #used as temporary container for the list within the larger list of data
            analysis = ''
            while 0 < len(data): #as long as the length of data is larger than 0
                scan = (data[0])#grab the 0th index from the list
                if scan[4] != scan[self.mod]: #if indexes do not match
                    if scan[self.mod] == 'iMet':#and making sure that index value is not iMet
                        pass
                    elif scan[self.mod] == 'fMet':#and that it is not fMet
                        pass
                    else:
                        analysis += scan[0] +'\t'+ scan[4] +'\t\t'+ scan[self.mod] +'\n' #add the indecies to analysis
                        scan = ''#empties scan
                        data.pop(0)#pops 0th index making 1st index the new 0th
            return analysis#returns the strings

    def grab_numbers(self):
        '''
        Create list of Scores
        '''
        if self.filt == 'inf':
            data = self.parser()#set function parser output to data
            transfer = []#empty lists used for transfering contents
            inf_score = []
            inf_scores = []
            while 0 < len(data): #while lenght of list data is larger than 0
                transfer.append(data[0]) #take the 0th item and append it to transfer generating a list of lists
                analyze = [entry[8] for entry in transfer] #for the list within list transfer grab the 8th index and save it to name analyze
                inf_score.append(analyze) #append analyze to list inf_score
                transfer.pop(0)# pop 0th index in lists data and transfer
                data.pop(0)
            if 0 == len(data): #once the length of data is equal to 0
                inf_score = [val for sublist in inf_score for val in sublist] #flatten the list of lists
                inf_scores = [float(i) for i in inf_score] #for each entry in list inf_score convert the string values to intgers.
            return inf_scores

        elif self.filt == 'iso':
            data = self.parser()#set function parser output to data
            transfer = []#empty lists used for transfering contents
            inf_score = []
            inf_scores = []
            while 0 < len(data): #while lenght of list data is larger than 0
                transfer.append(data[0]) #take the 0th item and append it to transfer generating a list of lists
                analyze = [entry[-1] for entry in transfer] #for the list within list transfer grab the 8th index and save it to name analyze
                inf_score.append(analyze) #append analyze to list inf_score
                transfer.pop(0)# pop 0th index in lists data and transfer
                data.pop(0)
            if 0 == len(data): #once the length of data is equal to 0
                inf_score = [val for sublist in inf_score for val in sublist] #flatten the list of lists
                inf_scores = [float(i) for i in inf_score] #for each entry in list inf_score convert the string values to intgers.
            return inf_scores

    def restrictions(self):
        '''Fetch values from functions grab_numbers, isoType, and tRNA_type and zip them'''
        inf_scores = self.grab_numbers() #set output from grab numbers to inf_scores
        isotype = self.isoType() #set output from isoType to name isotype
        if self.option == 'aa':
            tRNA = self.tRNA_type()#set ouput from tRNA_type to name tRNA
        elif self.option == 'anticodon':
            tRNA = self.antiCodon_type()#set ouput from antiCodon_type to name tRNA
        flint = zip (tRNA, isotype, inf_scores) #zips the lists
        return list(flint)#generates list of tuples

    def plot_filter(self):
        '''
        Fetch list of tuples and compare the third entry of each tuple
        (the Inf score from input file) and determine if the value is equal to
        or greater than the user's cut off value.
        '''
        data = self.restrictions() #set output from restrictions to data
        transfer = [] #empty container used for transfering data
        plot = [] #will become a list of lists containing points that will be plotted
        while 0 < len(data):#while the length of data is greater than 0
            transfer.append(data[0])#append the 0th element to transfer
            analyze = [entry[2] for entry in transfer] #take the 2nd element in the list within list transfer
            if self.user_value <= analyze[0]:# if the vaule is greater than or equal to the users value
                plot.append(data[0])#append it to plot
            transfer.pop(0)# pop the 0th index in transfer and data.
            data.pop(0)
        return plot #returns

    def xtRNA(self):
        '''
        Grab the tRNA type numeric vaule from the list of tuples that made the
        cut off done in plot_filter.
        '''
        data = self.plot_filter()
        transfer = []
        xtRNA = []
        while 0 < len(data):
            transfer.append(data[0])
            analyze = [entry[0] for entry in transfer]
            xtRNA.append(analyze)
            transfer.pop(0)
            data.pop(0)
        if 0 == len(data): #once length of list data is zero
            xtRNA = [val for sublist in xtRNA for val in sublist] #flatten nested list
        return xtRNA

    def yisoType(self):
        '''
        Grab the yisotype numeric value from the list of tuples that made the
        cut off done in plot_filter
        '''
        data = self.plot_filter()
        transfer = []
        yisotype = []
        while 0 < len(data):
            transfer.append(data[0])
            analyze = [entry[1] for entry in transfer]
            yisotype.append(analyze)
            transfer.pop(0)
            data.pop(0)
        if 0 == len(data):
            yisotype = [val for sublist in yisotype for val in sublist]
        return yisotype

    def plot(self):
        '''
        Fetch numeric data from functions tRNA_type and isoType. Plot points
        on a 2D plane.
        '''
        plt.rc('font', size = 14)
        x = self.xtRNA() #points that will be plotted for the x axis
        y = self.yisoType() #points that will be plotted for the y axis
        if self.option == 'aa':
            x2 = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22]
            y2 = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22]
            #x2 and y2 are used for the "normal line"
            x3 = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25]
            y3 = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23]
            #x3 and y3 are used for labeling the X and Y axis as amino acids in
            #conjunciton with labelx and labely
            labelx = ['Ala','Arg','Asn','Asp','Cys',
                'Gln', 'Glu', 'Gly','His', 'Ile',
                'Leu','Lys','Met','Phe', 'Pro',
                'Ser','Thr','Trp','Tyr', 'Val',
                'SeC', 'Pyl','Ile2', 'Undet',
                'Sup']
            labely = ['Ala','Arg','Asn','Asp','Cys',
                'Gln', 'Glu', 'Gly','His', 'Ile',
                'Leu','Lys','Met','Phe', 'Pro',
                'Ser','Thr','Trp','Tyr', 'Val',
                'SeC', 'Pyl', 'Ile2']
            plt.plot(x2,y2)#generatees "normal line"
        elif self.option =='anticodon':
             #Ala 1-4, Arg 5-10, Asn 11-12, Asp 13-14, Cys 15-16, Gln 17-18, Glu 19-20, Gly 21-24
             #His 25-26, Ile 27-29, leu 30-35, lys 36-37, Met 38, phe, 39-40, Pro 41-44, Ser 45-50
             #Thr 51- 54, Trp 55, Tyr 56-57, Val 58-61, SeC 62, Pyl 63,
            x3 = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,
                29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47
                ,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,64,65]
            y3 = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23]

            #x3 and y3 are used for labeling the X and Y axis as amino acids in
            #conjunciton with labelx and labely

            x4 = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,
                29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47
                ,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63]

            y4 = [1,1,1,1,2,2,2,2,2,2,3,3,4,4,5,5,6,6,7,7,8,8,8,8,9,9,10,10,10,
            11,11,11,11,11,11,12,12,13,14,14,15,15,15,15,16,16,16,16,16,16,17,17,17,17,
            18,19,19,20,20,20,20,21,22]

            #x4 and y4 are used to generate a "normal line" in the form of a
            #step-wise function showing the expected isotype in regards to the
            #anticodon Sequence

            labelx = ['AGC', 'GGC', 'UGC', 'CGC', 'ACG', 'GCG', 'UCG',
                    'CCG', 'UCU', 'CCU', 'AUU', 'GUU', 'AUC', 'GUC',
                    'ACA', 'GCA', 'UUG', 'CUG', 'UUC', 'CUC', 'ACC', 'GCC',
                    'UCC', 'CCC', 'AUG', 'GUG', 'AAU', 'GAU', 'UAU', 'UAA',
                    'CAA', 'AAG', 'GAG', 'UAG', 'CAG', 'UUU', 'CUU', 'CAU',
                    'AAA', 'GAA', 'AGG', 'GGG', 'UGG', 'CGG', 'AGA', 'GGA',
                    'UGA', 'CGA', 'ACU', 'GCU', 'AGU', 'GGU', 'UGU', 'CGU',
                    'CCA', 'AUA', 'GUA', 'AAC', 'GAC', 'UAC', 'CAC', 'UCA',
                    'CUA', 'UUA', 'NNN']

            labely = ['Ala','Arg','Asn','Asp','Cys',
            'Gln', 'Glu', 'Gly','His', 'Ile',
            'Leu','Lys','Met','Phe', 'Pro',
            'Ser','Thr','Trp','Tyr', 'Val',
            'SeC', 'Pyl', 'Ile2']
            plt.xticks(x, labelx, rotation='vertical')
            #names used for x and y axis
            plt.step(x4,y4, lw = 2)
        plt.xlabel('Anticodon')
        plt.ylabel('Isotype')
        plt.title(self.title)#title of graph is global name fname
        plt.xticks(x3, labelx)
        plt.yticks(y3, labely)
        #labels the x and y axis' as amino acids.
        plt.scatter(x,y, color = 'green') #generates data points using values from names x and y
        #accoidates for low frequency of apperances.
        plt.scatter (x,y, alpha=.1, s = 400, color = 'm' ) #generates data points
        #using data from names x and y. Each point is transparent and becomes
        #more visible as points overlap. Selected color is magenta.
        plt.grid(b = True, which = 'major', linestyle = '-' ) #creates darker
        #grid lines
        plt.show()

def main (myCommandLine = None):
    if myCommandLine is None: #default parameters
        fname = input('State tRNAs to Analyze: ')
        xfile = fname.replace('.out','').replace('.txt','')
        select = input('Select Amino Acid or Anticodon for X-Axis: ')
        option = select.lower()
        user_value = float(input('State Minimum Inf Score or select \"0\": '))
        title = xfile
        mod = 10
        computing = Module_Scanner(xfile, user_value, option, title, mod)
        suplemental = input('Output chimeric tRNAs?[y/N]: ').lower()
        if suplemental == 'y' or 'yes':
            saveFile = open(xfile+'-suplemental data.txt', 'w')#standardize tabs by using length
            saveFile.write('Sequence Name\t\t\t\t\ttRNA Type\tIsotype\n')
            saveFile.write('-------------\t\t\t\t\t---------\t-------\n')
            saveFile.write(computing.analysis())
            saveFile.close()


        computing.plot()

    else: #additional options for users if they choose to utilize them
        myCommandLine = CommandLine(myCommandLine)
        fname = myCommandLine.args.inFile
        xfile = fname.replace('.txt','').replace('.out','')
        select = myCommandLine.args.option
        if select == None:
            select = 'aa'
        option = select.lower()
        user_value = myCommandLine.args.infscore
        if user_value == None:
            user_value = 0

        title = xfile
        if myCommandLine.args.title != None:
            title = myCommandLine.args.title

        mod = 10
        if myCommandLine.args.mod != None:
            mod = myCommandLine.args.mod - 1

        supp = myCommandLine.args.supplemental
        if supp == 'y':
            saveFile = open(fname+'-suplemental data.txt', 'w')#standardize tabs by using length
            saveFile.write('Sequence Name\t\t\t\t\ttRNA Type\tIsotype\n')
            saveFile.write('-------------\t\t\t\t\t---------\t-------\n')
            saveFile.write(computing.analysis())
            saveFile.close()

        filt = myCommandLine.args.filter
        if filt == None:
            filt = 'inf'

        computing = Module_Scanner(fname, user_value, option, title, mod, supp, filt)

        computing.plot()

if __name__ == "__main__":
    if (sys.argv[1:]):
        main(sys.argv[1:])
    else :
        main()
