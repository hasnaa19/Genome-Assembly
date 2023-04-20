"""
Function:
    getReads
Description:
    Generic for (paired and single)
    Read First line get K and G (Case Paired)
    Time Complexity: O(1)
Input:
    FileName: String carries the text fileName
    Paired: boolean carries 0 for single reads and 1 for paired reads
Returns:
    Either K-mers, gap and Reads for Paired Reads or K-mers and Reads for Single Reads
Invoked/Called:
    Called at the beginning in Main
"""
def getReads(fileName, Paired):
    Reads = open(fileName, "r")

    line = Reads.readline()
    line = line.rstrip().split(" ")

    k = int(line[0])

    # For paired read the gap and return it for assembly
    if Paired:
        g = int(line[1])
        return k, g, Reads
    # Single Reads
    else:
        return k, Reads

'''
Function:  
    prefixExists
Description: 
    Generic for (paired and single)
    Checks if the prefix (aka key of the dictionary used to store reads) exists before i.e. prefix has more than one suffix.
    thus append the other suffix (where suffix is a list). if prefix doesn't exist then add it & it's first suffix element in 
    the list to the dictionary.
    Time Complexity: O(1)
Input: 
    Dictionary that is filled with { prefix : suffix } called "Graph"
    string carriers prefix used as the key for dict called "prefix"
    list of suffix carries suffix called "suffix"
Returns: 
    Graph
Invokes/Called: 
    Called in DebruijinGraph Function
'''
# Creating Graph checks if prefix has 2 suffix
def prefixExists(Graph, prefix, suffix):
    if prefix in Graph:
        Graph[prefix].append(suffix)
    else:
        Graph[prefix] = []
        Graph[prefix].append(suffix)
    return Graph

'''
Function: 
    Debruijn Graph
Description: 
    Generic Function for both single and paired
    Creates Graph using Splitting and Slicing Techniques
    Data Structure for graph -> used dictionary where value is a list
    Time Complexity -> O(N), where n is # of reads i.e. lines in file - 1
Input:
    Paired boolean
    Kmers and Reads 
Returns:
    Table/ Graph in a dictionary
Invokes/Called:
    Invokes: prefixExists
    Called: Main
'''
def DebruijnGraph(Paired, Kmers, Reads):
    dGraph = {}
    for Read in Reads:
        Read = Read.rstrip()
        # Paired Reads
        if Paired:
            Read = Read.split("|")
            prefix = "|".join([Read[0][0:Kmers - 1], Read[1][0: Kmers - 1]])
            suffix = "|".join([Read[0][1:Kmers], Read[1][1: Kmers]])
            dGraph = prefixExists(dGraph, prefix, suffix)
        # Single Reads
        else:
            prefix = Read[0: Kmers - 1]
            suffix = Read[1:Kmers]
            dGraph = prefixExists(dGraph, prefix, suffix)

    return dGraph

'''
Function: 
    eulerian_Path
Description: 
    Generic Function (paired and single) for finding path 
    Time Compexity -> O(2N) , N = dictionary length (keys)
Input:
    Graph dictionary
Returns:
    list carries the path of the reads called "walk"
Invokes/Called:
    Called in Main
'''
# Finding path
def eulerian_Path(graph):
    for key in graph.keys():
        is_unique = True
        for val in graph.values():
            if key in val:
                is_unique = False
                break
        if is_unique:
            start = key
            break
    walk = [start]
    while walk[-1] in graph.keys():
        if len(graph[str(walk[-1])]) == 1:
            # converted to string because lists are unhashable (i.e you can't use lists as a key to access dictionary
            val = graph.pop(walk[-1]).pop()
            # first pop is to remove corresponding key,value pair from dictionary while returning the value only without the key, the second pop is used to get the string from list
        else:
            if walk[-1] in graph[walk[-1]]:
                val = walk[-1]
                graph[walk[-1]].remove(walk[-1])
            else:
                val = graph[walk[-1]].pop()
        walk.append(val)
    return walk

'''
Function:  
    Paired Assemble
Description: 
    Get Assembly for paired using Slicing Technique
    Time Complexity -> O(M); M = Path length
Input:
    Path (list)
Returns:
    3 strings carries prefix, suffix and Genome
Invokes/Called:
    called in Main
'''
def Paired_Assemble(path,gap):
    Prefix = []
    Suffix = []
    for node in path:
        if node != path[len(path) - 1]:
            Prefix.append(node[0])
            Suffix.append(node[k])
        else:
            Prefix.append(node[0:k - 1])
            Suffix.append(node[k:])
    Prefix = ''.join(Prefix)
    Suffix = ''.join(Suffix)
    Genome = Prefix + Suffix[len(Suffix) - (k + gap):]
    return Prefix, Suffix, Genome

'''
Function: 
    Single_Assemble
Description: 
    Time Complexity -> O(M); M = lenght(path) 
Input: 
    list of Path
Returns: 
    String of genome
Invokes/Called: 
    Called in Main
'''
def Single_Assemble(path):
    Genome = ""
    Genome += str(path[0])
    for i in range(1, len(path)):
        temp = path[i]
        Genome += temp[-1]
    return Genome


# ------------------------------------------------- MAIN  -------------------------------------------------

paired = int(input("For Single Reads type 0 \t \t \t For Paired Reads type 1\n"))

if paired:
    FileName = "ReadPairsInputV1.txt"
    k, g, Reads = getReads(FileName, paired)

    print("\nDebruijn Graph : ")
    Graph = DebruijnGraph(paired, k, Reads)
    [print(key, '->', ",".join(value)) for key, value in Graph.items()]
    # print(dict(list(Graph.items())[0:1]))  # Checking Pairs

    print("\nEulerian path : ")
    path = eulerian_Path(Graph)
    [print(' -> ' + node, end='') for node in path]

    print("\n\nAssembly: ")
    Prefix, Suffix, Genome = Paired_Assemble(path, g)
    print("Genome :\n", Genome)
    print("Prefix :\n", Prefix)
    print("Suffix :\n", Suffix)

else:
    FileName = "SingleReadsInputV1.txt"
    k, Reads = getReads(FileName, paired)

    print("\nDebruijn Graph : ")
    Graph = DebruijnGraph(paired, k, Reads)
    [print(key, '->', ",".join(value)) for key, value in Graph.items()]

    print("\nEulerian path : ")
    path = eulerian_Path(Graph)
    [print(' -> ' + node, end='') for node in path]

    print("\n\nAssembly: ")
    Genome = Single_Assemble(path)
    print("Genome :\n", Genome)

