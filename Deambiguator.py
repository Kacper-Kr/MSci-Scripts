import os
import sys
import random
print(sys.executable)

from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction

def MakeDatabase(): #Makes the database with all of the sequences.
    # List of strain names
    sequences = ["wAlbB-HN2016", "wAlbB-Uju", "wAnD", "wAnM", "wBm", "wCauB", "wCfeT", "wHa", "wKueWO", "wMel", "wMelCS112", "wPip", "wRi","WO-1","WO-2","WO-3","WO-4","wmk-1","wmk-2","wmk-3","wmk-4","cifA","cifB"]
    # Dictionary to store SeqRecord objects with strain names as keys
    global Wolbachia
    Wolbachia = {}
    # Loop through each strain and read its sequence
    for strain in sequences:
        print(f"\nReading {strain} Fasta")
        # Parse the sequence for the current strain
        for record in SeqIO.parse(f'Fasta Sequences/{strain}.fasta', "fasta"):
            # Stores the SeqRecord in the dictionary
            Wolbachia[strain] = record

# To access a sequence information, use Wolbachia[strain].x, where x can be id (accession number), description or seq (for sequence).
# E.g. Wolbachia["wMel"].seq gives the wMel sequence
# With a sequence, translate by ammending .translate(to_stop={bool}), .complement(), .count(subsequence), .find(subsequence)
# Certain functions word with the sequence; len(sequence), gc_fraction(sequence), 


def TriTranslate(sequence): # Translates the nucleotide sequence into proteins based on the three possible reading, returning an array containing all three possibilities
    translations = []
    for i in range(0,3): # Repeats three times, for each potential set of codons
        while not len(sequence)%3 ==0: # Checks if the sequence is divisible by 3. If not, adds an ambiguous character to the end of the sequence to shift the frame.
            sequence+= "N"
        translations.append(sequence.translate(to_stop=False))
        sequence = sequence[1:-1] #Trims one base from the beginning of the sequence to account for alternative reading frames
    return translations

def contains_ambiguous(sequence): # Checks if a sequence contains any ambiguous bases
    # Define ambiguous bases
    ambiguous_bases = {"N", "V", "H", "D", "B", "M", "K", "W", "S", "Y", "R"}
    # Check if any character in the sequence is ambiguous
    for nucleotide in sequence:
        if nucleotide in ambiguous_bases:
            return True  # Found an ambiguous base
    return False  # No ambiguous base found

def Disambiguate(sequence):  #Replaces any ambiguous bases in the sequence with a base at random (using weights based on the A,T,G,C contents of the genome.
    # Count base frequencies in the sequence
    A = sequence.count("A")
    C = sequence.count("C")
    G = sequence.count("G")
    T = sequence.count("T")
    log.append(strain +" genome consists of " +str(len(sequence)) +" bases: A ("+str(round((A/(A + C + G + T))*100,2))+"%), C ("+str(round((C/(A + C + G + T))*100,2)) +"%), G (" +str(round((G/(A + C + G + T))*100,2)) +"%), T ("+str(round((T/(A + C + G + T))*100,2)) +"%).")
    try: #Checks for replicates; If replicates are being made this is logged, so you can tell what happened at each replicate.
        if i >= 0:
            log.append("Replicate " +str(i+2))
    except:
        "No Replicates Needed"
    # Dictionary for ambiguous base mappings
    ambiguous_bases = {
        "N": (["A", "C", "G", "T"], [A, C, G, T]),
        "V": (["A", "C", "G"], [A, C, G]),
        "H": (["A", "C", "T"], [A, C, T]),
        "D": (["A", "G", "T"], [A, G, T]),
        "B": (["C", "G", "T"], [C, G, T]),
        "M": (["A", "C"], [A, C]),
        "K": (["G", "T"], [G, T]),
        "W": (["A", "T"], [A, T]),
        "S": (["C", "G"], [C, G]),
        "Y": (["C", "T"], [C, T]),
        "R": (["A", "G"], [A, G]),
    }

    #Contiguate long stretches of Ns by replacing ALL Ns in that stretch into a single $ (which will be used to define a contig break in another function)
    
    # Resulting sequence
    sequencelist = []
    # Loop through the sequence
    counter = -1
    for nucleotide in sequence:
        ambiguous = False
        contiguate = False
        counter = counter + 1
        #Checks for many Ns, and if at least 10 are present in a row they are all replaced with a single $.
        if nucleotide == "N":
            if len(sequencelist) > 0: #Ensures (sequencelist) isn't empty, first (Prevents error)
                if sequencelist[-1] == "$": #Checks if contiguation just took place (this means Ns are part of this unknown segment) and updates log on end of contig break.
                    log[-1] = (strain +": Long ambiguous region with a length of " +str(counter - contigstart) +" bases, ends at " +str(counter+1) +".")
                    contiguate = True
                elif len(sequence)-counter > 49: #Makes sure there is at least 49 bases left at the end of the sequence (Prevents error).
                    if len(set(sequence[counter+1:counter+49])) == 1 and sequence[counter+1] == "N":
                        sequencelist.append("$") #A $ is added, which denotes contiguation when the sequence is exported
                        print("Long ambiguous region detected spanning at least 50 bases, which usually denotes unknown length; splitting into contigs.")
                        log.append(strain +": Long ambiguous region starts at " +str(counter+1) +"; Contiguating.") #Start of cotig break is logged
                        log.append(strain +": Long ambiguous region with a length of 50 bases, ends at " +str(counter+1) +".") #End of contig break is logged
                        contigstart = counter
                        contiguate = True
            elif len(sequence)-counter > 49: #Makes sure there is at least 10 bases left at the end of the sequence (Prevents error).
                if len(set(sequence[counter+1:counter+49])) == 1 and sequence[counter+1] == "N":
                        sequencelist.append("$") #A $ is added, which denotes contiguation when the sequence is exported
                        print("Long ambiguous region detected spanning at least 50 bases, which usually denotes unknown length; splitting into contigs.")
                        log.append(strain +": Long ambiguous region starts at " +str(counter+1) +"; Contiguating.") #Start of contig break is logged
                        log.append(strain +": Long ambiguous region with a length of 50 bases, ends at " +str(counter+1) +".") #End of contig break is logged
                        contiguate = True
        if nucleotide in ambiguous_bases and contiguate == False:  # Skip if not one of the valid ambiguous nucleotides (or if contiguation is required which means N nucleotides are skipped)
            # Retrieve the base options and weights for the ambiguous nucleotide
            bases, weights = ambiguous_bases.get(nucleotide, ([], []))
            old_nucleotide = nucleotide
            ambiguous = True
            global nonreplacing
            nonreplacing = False #Changes non-replacing to false, since a base had to be replaced randomly
            # If the sum of weights is zero, assign equal weights
            if sum(weights) == 0:
                weights = [1] * len(bases)
            # Choose a base based on weights
            nucleotide = random.choices(bases, weights=weights, k=1)[0]
        if contiguate == False: #Only adds bases if contiguation didn't occur (Prevents adding Ns when contiguation occurred)
            # Append the disambiguated nucleotide to the list
            sequencelist.append(nucleotide) #Nucleotide is added to the list (even if this didnt involve randomization)
            if ambiguous == True: #HOWEVER: If randomization occurred, this is logged.
                print(f"Replacing ambiguous base: {old_nucleotide} with {nucleotide}")
                log.append(strain +": Replacing ambiguous base "+old_nucleotide+" at " +str(counter+1) +" with " +nucleotide +".")
    # Return the disambiguated sequence as a Bio.Seq object
    log.append("")
    return Seq("".join(sequencelist))

def FastaExport(exportee):
    """
    Export sequences stored in the exportee dictionary to individual FASTA files.
    If a sequence contains a '$' (denoting a contig break), split the sequence
    into multiple contigs and write each contig as a separate entry.
    Each header will include accession number, description, and contig info if applicable.

    Args:
        exportee (dict): Dictionary where keys are strain names and values are SeqRecord objects.
        output_directory (str): Directory where the output FASTA files will be saved.
    """
    

    for strain, record in exportee.items():
        output_directory = r"C:\Users\btb20169\OneDrive - University of Strathclyde\MSci\Fasta Sequences\Exported"
        # Ensure output directory exists
        os.makedirs(output_directory, exist_ok=True)
        if record.repbool == True: #Checks if there are replicates, If there are:
            output_directory += f"\\{strain}" #Replicates are put in a folder with the strain name.
            # If folder doesn't exist, create it.
            if not os.path.exists(output_directory):
                os.makedirs(output_directory)
            counter = 0
            for replicate in record.rep: #Iterate through each replicate in the array, to export a separate file each.
                counter += 1
                # Generate output file path for this strain with "_Exported_#" appended
                output_file = os.path.join(output_directory, f"{strain}_Exported_" +str(counter) +".fasta")
        
                records = []  # List to hold SeqRecord objects for this strain
                # Split the sequence at '$' to get individual contigs
                contigs = str(replicate).split("$")
                for i, contig in enumerate(contigs):
                    if contig.strip():  # Skip empty contigs resulting from consecutive '$'
                        # Generate the header based on whether it's a single contig or multiple
                        contig_id = record.id
                        # Description format: "Accession Number Description, contig X"
                        description = f"{record.description}, contig {i+1}/"+str(len(contigs)) if len(contigs) > 1 else f"{record.description}"
                
                        # Create SeqRecord for the contig
                        contig_record = SeqRecord(Seq(contig), id=contig_id, description=description)
                        records.append(contig_record)

                # Write all records to the output FASTA file for this strain
                with open(output_file, "w") as fasta_out:
                    SeqIO.write(records, fasta_out, "fasta")
        
                print(f"\nSequences for {strain} replicate " +str(counter) +f" have been exported to {output_file}")
        else:
            # Generate output file path for this strain with "_Exported" appended
            output_file = os.path.join(output_directory, f"{strain}_Exported.fasta")
        
            records = []  # List to hold SeqRecord objects for this strain
            # Split the sequence at '$' to get individual contigs
            contigs = str(record.seq).split("$")
            for i, contig in enumerate(contigs):
                if contig.strip():  # Skip empty contigs resulting from consecutive '$'
                    # Generate the header based on whether it's a single contig or multiple
                    contig_id = record.id
                    # Description format: "Accession Number Description, contig X"
                    description = f"{record.description}, contig {i+1}/"+str(len(contigs)) if len(contigs) > 1 else f"{record.description}"
                
                    # Create SeqRecord for the contig
                    contig_record = SeqRecord(Seq(contig), id=contig_id, description=description)
                    records.append(contig_record)

            # Write all records to the output FASTA file for this strain
            with open(output_file, "w") as fasta_out:
                SeqIO.write(records, fasta_out, "fasta")
        
            print(f"\nSequences for {strain} have been exported to {output_file}")

log = []        
MakeDatabase()
#Automatically disambiguates all sequences that need it, running 10 replicates of them.
for strain in Wolbachia:
    Wolbachia[strain].repbool = False #Creates a variable to store if replicates are needed.
    print(f"\nAnalyzing {strain}")
    if contains_ambiguous(strain) == True: #If the strain contains ambiguous bases...
        nonreplacing = True
        replicates = [Disambiguate(Wolbachia[strain].seq)] #Creates initial replicate, which is checked for presence of randomization
        if nonreplacing == False: #Checks if randomization isn't required (e.g. ambiguous bases caused by long sequences of N, which are contiguated). Otherwise:
            for i in range(0,9): #Creates 9 more replicates
                print("\nCreating Replicate "+str(i+2) +f" for {strain}")
                replicates.append(Disambiguate(Wolbachia[strain].seq))
            i = -1
            Wolbachia[strain].rep = replicates
            Wolbachia[strain].repbool = True #Boolean is changed, since replicates are present and Wolbachia[strain].seq will now be a different format (Array)
        else:
            Wolbachia[strain].seq = replicates[0] #Replicates not required so initial replicate is set to the final value          
        
    else:
        Wolbachia[strain].seq = Disambiguate(Wolbachia[strain]) #If disambiguation isn't required, all strains are run through the function anyway (since it calculates ATGC content)
FastaExport(Wolbachia) #Strains are exported


output_directory = r"C:\Users\btb20169\OneDrive - University of Strathclyde\MSci\Fasta Sequences\Exported"
# Ensure output directory exists
os.makedirs(output_directory, exist_ok=True)
    
# Creates a log file, which tells the user what the program did
with open(output_directory+"\log.txt", 'w') as file:
    for item in log:
        file.write(item + '\n')  # Write each item followed by a newline
print(f"\nLog file has been exported to {output_directory}\log.txt")




