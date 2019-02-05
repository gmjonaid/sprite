#!/usr/bin/env python

import time
import os
import sys




config = {}
reversedDpmDict = {}
with open("config_1.txt", "r") as c:
	for line in c:
		_, tag, seq, _ = line.rstrip("\n").split("\t")
		config[tag] = seq
		if "DPM" in tag:
			reversedDpmDict[seq] = tag

print("CONFIG READING..")
#dpmDict = { k:v for k, v in config.items() if "DPM" in k}
yDict = {k:v for k, v in config.items() if "Ybot" in k}
oddDict = {k:v for k, v in config.items() if "Odd" in k}
evenDict = {k:v for k, v in config.items() if "Even" in k}
#reversedDpmDict = { v:k for k, v in dpmDict.items()}
reversedYDict = { v:k for k, v in yDict.items()}
reversedOddDict = { v:k for k, v in oddDict.items()}
reversedEvenDict = { v:k for k, v in evenDict.items()}
maY = max([len(i) for i in reversedYDict])
miY = min([len(i) for i in reversedYDict])	
dif_y = maY - miY +1
lo = max([len(i) for i in reversedOddDict]) #length of Odd sequences; same for all
le = max([len(i) for i in reversedEvenDict])

print("CONFIG PROCESSED..")


def FastqGeneralIterator(handle):
    handle_readline = handle.readline
    line = handle_readline()
    if not line:
        return  
    if isinstance(line[0], int):
        raise ValueError("Is this handle in binary mode not text mode?")

    while line:
        if line[0] != "@":
            raise ValueError("Records in Fastq files should start with '@' character")
        title_line = line[1:].rstrip()
        seq_string = handle_readline().rstrip()
        while True:
            line = handle_readline()
            if not line:
                raise ValueError("End of file without quality information.")
            if line[0] == "+":
                second_title = line[1:].rstrip()
                if second_title and second_title != title_line:
                    raise ValueError("Sequence and quality captions differ.")
                break
            seq_string += line.rstrip()  

        if " " in seq_string or "\t" in seq_string:
            raise ValueError("Whitespace is not allowed in the sequence.")
        seq_len = len(seq_string)
        quality_string = handle_readline().rstrip()
        while True:
            line = handle_readline()
            if not line:
                break  
            if line[0] == "@":
                if len(quality_string) >= seq_len: 
                    break
                
            quality_string += line.rstrip()

        if seq_len != len(quality_string):
            raise ValueError("Lengths of sequence and quality values differs "
                             " for %s (%i and %i)."
                             % (title_line, seq_len, len(quality_string)))


        yield (title_line, seq_string, quality_string)





def dpmFinder(sequence, readName1):
	
	seq = sequence[:8]
	if seq in reversedDpmDict:
			readName1 = readName1 + "::" +"[" + reversedDpmDict[seq] + "]"
	else:
		readName1 = readName1 + "::" + "[NOT_FOUND]"
	return readName1

def seqWalk(Seq, s, m, l):
	seqs = [Seq[s:i+m] for i in range(l)]
	return seqs

def laxSeqWalk( Seq, sl,l): #l for laxity
	seqs = [Seq[i:i+sl] for i in range(l)]
	return seqs



def yFinder(sequence, readName2):
	seq = sequence[:maY]
	for j in seqWalk(seq, 0, miY, dif_y ):
		if j in reversedYDict:
			readName2 = readName2 + "::" + "[" + reversedYDict[j] + "]"
			break 
	else:
		readName2 = readName2 + "::" + "[NOT_FOUND]"
	return readName2 

# sp = spacer, l = laxity 

def oddFinder(sequence, readName2, sp, l):
	seq = sequence[miY+sp-1 :miY+sp+ lo + l+dif_y]
	for j in laxSeqWalk(seq, lo, l+dif_y):
		if j in reversedOddDict:
			readName2 = readName2 + "[" + reversedOddDict[j] + "]"
			break
	else: 
		readName2 = readName2 + "[NOT_FOUND]"
	return readName2

def evenFinder(sequence, readName2, sp, l):
	
	seq = sequence[miY+sp+lo+sp :miY+sp+lo+sp+ le + l + dif_y]
	for j in laxSeqWalk(seq, le, l+dif_y):
		if j in reversedEvenDict:
			readName2 = readName2 + "[" + reversedEvenDict[j] + "]"
			break
	else: 
		readName2 = readName2 + "[NOT_FOUND]"
	return readName2

def oddFinder2(sequence, readName2, sp, l):
	seq = sequence[miY+sp+lo+sp+le+sp : miY+sp+lo+sp+le+sp+lo + l+dif_y]
	for j in laxSeqWalk(seq, lo, l+dif_y):
		if j in reversedOddDict:
			readName2 = readName2 + "[" + reversedOddDict[j] + "]"
			break
	else: 
		readName2 = readName2 + "[NOT_FOUND]"
	return readName2

"""
if __name__ == "__main__":
	c=0
	fw = open(sys.argv[3], "w")
	print("OPENING FILE1 AND FILE2..")
	with open(sys.argv[1]) as f1, open(sys.argv[2]) as f2:
		for x, y in zip(FastqGeneralIterator(f1), FastqGeneralIterator(f2)):
			name1, seq1, qual1 = x 
			name2, seq2, qual2 = y
			read1 = name1.split(" ")[0] 
			read2 = name2.split(" ")[0]
			dpmRead = dpmFinder(seq1, read1)
			yRead = yFinder(seq2, read2)
			oRead = oddFinder(seq2, yRead, 6,9)
			eRead = evenFinder(seq2, oRead, 6,9)
			o2Read = oddFinder2(seq2, eRead, 6,9)
			finalRead = dpmRead + o2Read.split("::")[1]
			fw.write("@"+finalRead+"\n")
			fw.write(seq1+"\n")
			fw.write("+"+"\n")
			fw.write(qual1+"\n")
	fw.close()
	print("DONE!")
			
	

"""

		
def main():
	fw = open(sys.argv[3], "w")
	with open(sys.argv[1]) as f1, open(sys.argv[2]) as f2:
		while True:
			name1 = f1.readline().rstrip("\n")
			seq1 = f1.readline().rstrip("\n")
			f1l3 = f1.readline().rstrip("\n")
			qual1 = f1.readline().rstrip("\n")
			if not qual1:
				break 
			name2 = f2.readline().rstrip("\n")
			seq2 = f2.readline().rstrip("\n")
			f2l3 = f2.readline().rstrip("\n")
			qual2 = f2.readline().rstrip("\n")
			if not qual2:
				break 
			read1 = name1.split(" ")[0] 
			read2 = name2.split(" ")[0]
			dpmRead = dpmFinder(seq1, read1)
			yRead = yFinder(seq2, read2)
			oRead = oddFinder(seq2, yRead, 6,9)
			eRead = evenFinder(seq2, oRead, 6,9)
			o2Read = oddFinder2(seq2, eRead, 6,9)
			finalRead = dpmRead + o2Read.split("::")[1]
			fw.write(finalRead+"\n")
			fw.write(seq1+"\n")
			fw.write("+"+"\n")
			fw.write(qual1+"\n")
	fw.close()

if __name__ == "__main__":
	main()

		
