import collections
import random
from structures import *
class GeneToolKit:
    """ Gene Tool Kit """
    """ Default sequence = ATCG """
    """ Default sequenceType = DNA """
    
    def __init__(self,sequence="ATCG",sequenceType="DNA"):
        self.sequenceType = sequenceType.upper()
        self.sequence = sequence.upper()
        self.isvalid = self.validateSequence()
        assert self.isvalid, f"Provided sequence dosen't seem to be correct {self.type} type"
    
    def validateSequence(self):
        """ check if Sequence is Valid DNA or not """
        return set(NUCLEOTIDE_BASE[self.sequenceType]).issuperset(self.sequence)
    
    def info(self):
        """ Print general info about the sequence """
        print(f"\n[Sequence] : {self.sequence} \n[Type] : {self.type} \n[Length] : {len(self.sequence)} \n")


    def genRandomSequence(self,length=10,seqtype="DNA"):
        """ Generate a random DNA/RNA sequence with default length = 10 and default seqtype = DNA """
        self.sequence =  ''.join([random.choice(NUCLEOTIDE_BASE[seqtype]) for i in range(self.length)])


    def nucleotideFrequency(self):
        """ Count nucleotide in a given sequence, return dictionary """
        return dict(collections.Counter(self.sequence))
 
    def transcript(self):
        """ DNA => RNA Transcription Replacing Thymine with Uracil """
        if self.sequenceType == "DNA":
            return self.sequence.replace("T","U")
        else:
            print("Not a DNA sequence")

    def complement(self,reverse = False):
        """ Swapping adenin with thymine and guanine with cytosine 
            with default reverse=False """
        if self.sequenceType == "DNA":
            if reverse:
                return ''.join([DNA_ReverseComplement[nuc] for nuc in self.sequence])[::-1]
            else:
                return ''.join([DNA_ReverseComplement[nuc] for nuc in self.sequence])
        else:           
            if reverse:
                return ''.join([RNA_ReverseComplement[nuc] for nuc in self.sequence])[::-1]
            else:
                return ''.join([RNA_ReverseComplement[nuc] for nuc in self.sequence])
            
    def Color(self,seq):
        """ Used to give colors to DNA/RNA Sequences """
        bcolors = {
            'A' : '\033[92m',
            'C' : '\033[94m',
            'G' : '\033[93m',
            'T' : '\033[91m',
            'U' : '\033[91m',
            'reset' : '\033[0;0m'
        }
        tempstr = ""
        
        for nuc in seq:
            if nuc in bcolors:
                tempstr += bcolors[nuc] + nuc
            else:
                tempstr += bcolors['reset'] + nuc
        return tempstr + '\033[0;0m'
    
    def gcContent(self):
        """ GC Content in a DNA/RNA sequence """
        return round((self.sequence.count('C') + self.sequence.count('G')) / len(self.sequence) * 100)

    def gcContentSubsec(self,k = 20):
        """ GC Content in a DNA/RNA sub-sequence length k , k=20 by default """
        result = []
        for i in range(0,len(self.sequence)-k+1,k):
            subseq = self.sequence[i:i+k]
            result.append(round((subseq.count('C') + subseq.count('G')) / len(subseq) * 100))
        return result

    def translate_seq(self,init_pos=0):
        """ Translate DNA sequence into an aminoacids sequence """
        if self.sequenceType == "DNA":
            return [DNA_Codons[self.sequence[pos:pos+3]] for pos in range(init_pos,len(self.sequence)-2,3) ]
        elif self.sequenceType == "RNA":
            return [RNA_Codons[self.sequence[pos:pos+3]] for pos in range(init_pos,len(self.sequence)-2,3) ]

    def codon_usage(self,aminoacid):
        """ Provides the Frequency of each codon encoding a given aminoacid in a DNA sequence """
        templist = []
        if self.sequenceType == "DNA":
            for i in range(0,len(self.sequence)-2,3):
                if DNA_Codons[self.sequence[i:i+3]] == aminoacid:
                    templist.append(self.sequence[i:i+3])
        elif self.sequenceType == "RNA":
            for i in range(0,len(self.sequence)-2,3):
                if RNA_Codons[self.sequence[i:i+3]] == aminoacid:
                    templist.append(self.sequence[i:i+3])
                    
        freqdict = dict(collections.Counter(templist))
        totalWeight = sum(freqdict.values())
        for seq in freqdict:
            freqdict[seq] = round(freqdict[seq] / totalWeight ,2)
        return freqdict

    def genReadingFrames(self):
        """ Generating the six reading frames of a DNA sequence, including the reverse complement """
        frames = []
        frames.append(self.translate_seq(0))
        frames.append(self.translate_seq(1))
        frames.append(self.translate_seq(2))
        tempseq = GeneToolKit(self.complement(reverse = True),self.sequenceType)
        frames.append(tempseq.translate_seq(0))
        frames.append(tempseq.translate_seq(1))
        frames.append(tempseq.translate_seq(2))
        del tempseq
        return frames

    def protienFromReadingFrames(self,AASeq):
        """ Compute all possible protiens in an aminoacid sequence and return a list of all possible protiens """
        currentProtien = []
        protiens = []
        
        for aa in AASeq:
            if aa == '_':
                # Stop accumulating amino acids if stop codon was found
                if currentProtien:
                    for p in currentProtien:
                        protiens.append(p)
                    currentProtien = []
            else:
                # Start accumulating amino acids if Start codon was found
                if aa == 'M':
                    currentProtien.append("")
                for i in range(len(currentProtien)):
                    currentProtien[i] += aa
        return protiens

    def allProtiensFromOpenReadingFrames(self, startReadPos=0, endReadPos=0, ordered=False):
        """ Compute all possible protiens for all open reading frames """
        """ Protien search DB: https://www.ncbi.nlm.nih.gov/nuccore/NM_001185097.2 """
        """ API can be used to pull protien info """
        if endReadPos > startReadPos:
            tempseq = GeneToolKit(self.sequence[startReadPos:endReadPos],self.sequenceType)
            rfs = tempseq.genReadingFrames()
        else:
            rfs = self.genReadingFrames()
        
        result = []
        for rf in rfs:
            protiens = self.protienFromReadingFrames(rf)
            for p in protiens:
                result.append(p)
        if ordered:
            return sorted(result,key=len,reverse=True)
        return result