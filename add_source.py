# This script will take a list of the filenames (i.e., filenames.txt) to be appended with the SOURCE_ID. In this case, the SOURCE_ID appended to each protein sequence corresponds to the name of each fasta file minus the ".fasta". The script will loop through each protein sequence (identified by ">") in each fasta file (listed in "filenames.txt") to add the SOURCE_ID. The appended output will be saved at a new fasta file that starts with "new_".
def main():
    names = open("filenames.txt",'r')
    for file_name in names:
        file_name=file_name.strip("\n")
        file_out=open("new_"+file_name,'w')
        with open(file_name, 'r') as file:
            for line in file:
                if ">" in line:
                    file_name2=file_name.replace(".fasta","")
                    lineNew=line.rstrip() + " /SOURCE_ID="+file_name2
                    print(lineNew,file=file_out)
                else:
                    print(line,file=file_out)

if __name__ == "__main__":
    main()
