#I created the functions stats, find and dna2protein based on the assumption that my inputs will only be upper case
#the code for the case when we assume that they are not necessarily upper case is written in comments
#the functions were runned both as case sensitive and case insensitive and the last one, which includes transforming all in upper case showed a significant increase in the run time


def load(filename):

    """Takes the filename of a .fna file and returns the DNA sequence and description line"""
    
    description='file was not loaded'
    sequence=[]    #create an empty list
    
    if filename.endswith('.fna'): #we want to change description and sequence only if the file is a .fna file
        
        try:  #to avoid throwing an error in case the file cannot be open;
            #description and sequence will be changed just if filename is a .fna file that can be open
        
            fin=open(filename)
            fin.read(1) #because we don't want to include the first character '>' in the description
            description=fin.readline().strip() #takes the first line of the .fna file, without the '>' character, and remove leading and trailing whitespaces
            
            line =fin.readline().strip().upper() #need to define line before the loop in order to test the loop condition
            while line:
                for i in line:     
                    sequence.append(i)    #each character in the line will be appended to the sequence list
                line =fin.readline().strip().upper() #takes each line, removes the leading and trailing whitespaces (strip used especially to remove the new line character) and transforms all in upper case
            
            return sequence, description  #updated sequence and description

        except:
            return sequence, description  #sequence and description will be as initialised

    else:
        return sequence, description  #sequence and description will be as initialised
    fin.close()




def stats(sequence):

    """Takes a sequence and returns a table that includes the number of times a nucleic acid code occurs"""
    table={'A':0,'C':0,'G':0,'T':0,'N':0,'U':0,'K':0,'S':0,'Y':0,'M':0,'W':0,'R':0,'B':0,'D':0,'H':0,'V':0,'-':0} #create a dictionary and set all the values for each key equal to 0
    other=0 #we will count the characters that are not one of the 17 nucleic acid codes separately
    
    for m in sequence:
        if m in table:   #to make the function case insensitive, this is replaced by: if m.upper() in table:
            table[m]+=1 #increment value of the corresponding key    #to make the function case insensitive, this is replaced by: table[m.upper()]+=1
        else:
            other+=1 #count the number of times we find a character that is not key of the table
            
    table['other']=other #add the 'other' key with its corresponding value in the dictionary
    
    return table






def format_sequence(sequence,first_index,last_index):

    """Takes a sequence along with two indices and returns the subsequence with a particular format """

    (a,b)=divmod(last_index-first_index+1,80)  #a is the number of strings of 80 characters each; if b is greater than zero, we obtain one more string with b characters; this holds for last_index<=length of sequence
    formatted_sequences=[] 
    m=first_index
    group=''.join(sequence[m:m+80])  #we form the first group of characters; note that we don't necessarily need 80 characters for this
    
    for n in range (0,a):  #a=number of possible 80 characters groups; more efficient than going through each element
        if group!='':  #we may have last_index greater than the length of our sequence, so that the variable a will take a value greater than the actual number of groups of 80 characters, so we need to prevent appending empty strings
            formatted_sequences.append(group)  
            m=m+80 #increase m to be able to take the next group of elements
            group=''.join(sequence[m:m+80])  # create the next group of 80 characters
        else:
            break  #if we find an empty string, it means that we no longer have groups of 80 characters within our sequence, so we can break the loop
        
    if b!=0:  #if we have groups of less than 80 characters
        group=''.join(sequence[m:m+b])  #create the group containing the remaining characters
        if group!='': #check again if the group is an empty string, for example, we may have last_index greater than the length of the sequence, so all the possible groups of characters will be appended in the previous for loop
            formatted_sequences.append(group)  #add the string formed by concatenating the  b elements that remain after creating groups of 80 characters each

    return formatted_sequences





def write(filename, description,sequence,first_index,last_index):

    """Takes a description, sequence, and sequence range, and writes to a .fna file"""

    fout=open(filename,'w') #this will overwrite the file if the file exists
    fout.write(f'>{description}\n') #this will write the description line with the required format
    fout.write(''.join(sequence[first_index:last_index+1])) #write the subsequence of nucleic acids, use (last_index+1) to include sequence[last_index]
    fout.close()





def find(sequence,sequence_to_find):

    """Finds a sequence within another sequence and records the indices where they occured"""

    matches=[] #create an empty list, where we will append the indices where matches were found


    #to make the function case insensitive, add the following 2 lines:
    #sequence=[element.upper() for element in sequence]
    #sequence_to_find=[element.upper() for element in sequence_to_find]
    
    for n in range(0,len(sequence)-len(sequence_to_find)+1):   # go up to (len(sequence)-len(sequence_to_find)+1) to include index (len(sequence)-len(sequence_to_find)), where the last possible match may occur
                                                             #because index n increases by 1 with every step, matches will include overlapping occurances
        if sequence[n:n+len(sequence_to_find)]==sequence_to_find:   #check every possible group of consecutive characters with the same length as sequence_to find;
            matches.append(n)  #add the index where the match was found to the list

    return matches






def add(sequence,sequence_to_add,index):

    """Adds a sequence into an existing sequence at a specified index"""

    new_sequence=sequence[0:index]+sequence_to_add+sequence[index:] #create the new_sequence; also works for adding the sequence_to_add at the end of sequence in case index is beyond sequence

    #check:
    #if len(sequence)+len(sequence_to_add)==len(new_sequence):
        #print('the length of new_sequence is correct')
    #else:
        #print('incorrect')

    return new_sequence




        

def delete(sequence,index,number_of_codes):

    """Deletes a subsequence from a sequence as specified by a starting index and the number of codes to delete"""

    if index<len(sequence): #because the last element of sequence is sequence[len(sequence)-1]
        new_sequence=sequence[0:index]+sequence[index+number_of_codes:] #includes the case when the index and number_of_codes exceed the length of the sequence and we delete nucleic acid codes until the end of the sequence
    else:
        new_sequence=sequence #if index is out of bounds, the sequence remains unchanged, so the new_sequence takes the value of sequence

    return new_sequence



def replace(sequence,sequence_to_add,index,number_of_codes):

    """Replaces a section of sequence with a new subsequence"""

    incomplete_sequence=delete(sequence,index,number_of_codes) #delete number_of_codes elements from sequence, starting from sequence[index]
    new_sequence=add(incomplete_sequence,sequence_to_add,index) #add sequence_to_add to the incomplete_sequence before element on position index
    
    return new_sequence 




def dna2protein(dna_sequence):

    """Converts the DNA sequence to its corresponding protein sequence"""

    table={}
    protein_sequence=[]
    
    fin=open('dna2protein.csv')
    dna2protein_list=list(fin) #create a list with elements of type str corresponding to each line of the file fin
    
    for line in dna2protein_list:
        line=line.split(',') #split each line into its components (elements separated by ',') and transform the str line in a list containing these elements
        table[line[0]+line[1]+line[2]]=line[4] #associate each 3-letter nucleic acid code(entry) with its specific aminoacid(value) in a dictionary
    table['???']='?' #add the entry '???', which maps to the character '?' 

    for n in range (0,len(dna_sequence)//3): #use floor division to split the dna_sequence in disjoint pairs of 3 nucleic acid codes; this will not consider any remaining nucleic acid codes
        aminoacid=''.join(dna_sequence[n*3:(n+1)*3])  #we build the aminoacid corresponding to the next 3 letters in the dna_sequence; #to make the function case insensitive, this line is replaced by: aminoacid=''.join(dna_sequence[n*3:(n+1)*3]).upper()

        if aminoacid in table:
            protein_sequence.append(table[aminoacid]) #if the aminoacid is key for our table, we add to the list protein_sequence the value of our aminoacid, ie the 1-letter notation from the .csv file; this will include the entry '???'
        else:
            protein_sequence.append(table['???']) #if not found as entry, the 3-letter pair will be associated with '?', which is the value for '???'
        
    fin.close()
    return protein_sequence,table








