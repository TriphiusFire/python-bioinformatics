# import the random package here
import random

# Input:  Positive integers k and t, followed by a list of strings Dna
# Output: RandomizedMotifSearch(Dna, k, t)
def RandomizedMotifSearch(Dna, k, t):
    # insert your code here
    M = RandomMotifs(Dna, k, t)
    BestMotifs = M
    while True:
        Profile = ProfileWithPseudocounts(M)
        M = Motifs(Profile, Dna)
        if Score(M) < Score(BestMotifs):
            BestMotifs = M
        else:
            return BestMotifs

# Input:  A set of kmers Motifs
# Output: A consensus string of Motifs.
def Consensus(Motifs):
    # insert your code here
    k = len(Motifs[0])
    count = Count(Motifs)
    consensus = ""
    for j in range(k):
        m = 0
        frequentSymbol = ""
        for symbol in "ACGT":
            if count[symbol][j] > m:
                m = count[symbol][j]
                frequentSymbol = symbol
        consensus += frequentSymbol
    return consensus

# Input:  A set of kmers Motifs
# Output: Count(Motifs)
def Count(Motifs):
    count = {} # initializing the count dictionary
    # your code here
    k = len(Motifs[0])
    for symbol in "ACGT":
        count[symbol] = []
        for j in range(k):
            count[symbol].append(0)
    t = len(Motifs)
    for i in range(t):
        for j in range(k):
            symbol = Motifs[i][j]
            count[symbol][j] += 1
    return count

# Input:  A set of kmers Motifs
# Output: CountWithPseudocounts(Motifs)
def CountWithPseudocounts(Motifs):
    count = {} # initializing the count dictionary
    # your code here
    k = len(Motifs[0])          #how long a motif row is
    for symbol in "ACGT":       #for each nucleotide
        count[symbol] = []      #create a list that stores counts for each col
        for j in range(k):      #for each column
            count[symbol].append(1) #initialize to zero
    t = len(Motifs)             #how many rows in motif matrix
    for i in range(t):          #in each row
        for j in range(k):      #go throw each colum
            symbol = Motifs[i][j]   #Get the nucleotide
            count[symbol][j] += 1   #and increment its count for that column
    return count

# Input:  A list of kmers Dna, and integers k and t (where t is the number of kmers in Dna)
# Output: GreedyMotifSearch(Dna, k, t)
def GreedyMotifSearch(Dna, k, t):
    # type your GreedyMotifSearch code here.
    BestMotifs = []
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])
    n = len(Dna[0])
    for i in range(n-k+1):
        Motifs = []
        Motifs.append(Dna[0][i:i+k])
        for j in range(1, t):
            P = Profile(Motifs[0:j])
            Motifs.append(ProfileMostProbablePattern(Dna[j], k, P))
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs

# Input:  A list of kmers Dna, and integers k and t (where t is the number of kmers in Dna)
# Output: GreedyMotifSearch(Dna, k, t)
def GreedyMotifSearchWithPseudocounts(Dna, k, t):
    # type your GreedyMotifSearch code here.
    BestMotifs = []
    for i in range(0, t):
        BestMotifs.append(Dna[i][0:k])
    n = len(Dna[0])
    for i in range(n-k+1):
        Motifs = []
        Motifs.append(Dna[0][i:i+k])
        for j in range(1, t):
            P = ProfileWithPseudocounts(Motifs[0:j])
            Motifs.append(ProfileMostProbablePattern(Dna[j], k, P))
        if Score(Motifs) < Score(BestMotifs):
            BestMotifs = Motifs
    return BestMotifs

# Input:  A profile matrix Profile and a list of strings Dna
# Output: Motifs(Profile, Dna)
def Motifs(Profile, Dna):
    motifs = []
    for s in Dna:
        motifs.append(ProfileMostProbablePattern(s, len(Profile['A']), Profile))
    return motifs

# Input:  String Text and profile matrix Profile
# Output: Pr(Text, Profile)
def Pr(Text, Profile):
    # insert your code here
    p = 1
    for i in range(len(Text)):
        p = p * Profile[Text[i]][i]
    return p

# Input: A dictionary Probabilities, where keys are k-mers and values are the probabilities of these k-mers (which do not necessarily sum up to 1)
# Output: A normalized dictionary where the probability of each k-mer was divided by the sum of all k-mers' probabilities
def Normalize(Probabilities):
    sump = 0.0
    for p in Probabilities:
        sump += Probabilities[p]
    for p in Probabilities:
        Probabilities[p] = Probabilities[p]/sump
    return Probabilities

# Input:  A list of kmers Motifs
# Output: the profile matrix of Motifs, as a dictionary of lists.
def Profile(Motifs):
    profile = Count(Motifs)
    k = len(profile['A'])
    for c in profile:
        for j in range(k):
            profile[c][j] = round(float(profile[c][j]) / float(len(Motifs)),1)            
    return profile

# Input:  A string Text, a profile matrix Profile, and an integer k
# Output: ProfileGeneratedString(Text, profile, k)
def ProfileGeneratedString(Text, profile, k):
    n = len(Text)
    probabilities = {} 
    for i in range(0,n-k+1):
        probabilities[Text[i:i+k]] = Pr(Text[i:i+k], profile)
    probabilities = Normalize(probabilities)
    return WeightedDie(probabilities)

# Input:  A set of kmers Motifs
# Output: ProfileWithPseudocounts(Motifs)
def ProfileWithPseudocounts(Motifs):
    profile = CountWithPseudocounts(Motifs)
    t = len(profile)
    k = len(profile['A'])
    for c in profile:
        for j in range(k):
            profile[c][j] = profile[c][j] / (len(Motifs)+4)
    return profile

# Input:  String Text, an integer k, and profile matrix Profile
# Output: ProfileMostProbablePattern(Text, k, Profile)
def ProfileMostProbablePattern(Text, k, Profile):
    # insert your code here. Make sure to use Pr(Text, Profile) as a subroutine!
    mostProbable = ""
    maxPr = float(-1);
    for i in range(len(Text)-k+1):
        p = Pr(Text[i:i+k],Profile)
        if p > maxPr:
            mostProbable = Text[i:i+k]
            maxPr = p
    return mostProbable

# Input:  A list of strings Dna, and integers k and t
# Output: RandomMotifs(Dna, k, t)
# HINT:   You might not actually need to use t since t = len(Dna), but you may find it convenient
def RandomMotifs(Dna, k, t):
    # place your code here.
    r_list = []
    for s in Dna:
        index = random.randint(0,len(Dna[0])-k)
        r_list.append(s[index:index+k])
    return r_list

# Input:  A set of k-mers Motifs
# Output: The score of these k-mers.
def Score(Motifs):
    countMotif = Count(Motifs)               
    consensusMotif = Consensus(Motifs)   
    score = 0
    for i in range(len(consensusMotif)):
        for c in countMotif:
            if consensusMotif[i] != c:
                score += countMotif[c][i]
    return score

#input Dna, k, t
#output best motifs list --WITH N instead of 1000 see next function
# def RepeatedRandomizedMotifSearch(Dna, k, t):
#     BestScore = float('inf')
#     BestMotifs = []
#     for i in range(1000):
#         Motifs = RandomizedMotifSearch(Dna, k, t)
#         CurrScore = Score(Motifs)
#         if CurrScore < BestScore:
#             BestScore = CurrScore
#             BestMotifs = Motifs
#     return BestMotifs

#input Dna, k, t
#output best motifs list
def RepeatedRandomizedMotifSearch(Dna, k, t, N):
    BestScore = float('inf')
    BestMotifs = []
    for i in range(N):
        Motifs = RandomizedMotifSearch(Dna, k, t)
        CurrScore = Score(Motifs)
        if CurrScore < BestScore:
            BestScore = CurrScore
            BestMotifs = Motifs
    return BestMotifs

# Input:  A dictionary Probabilities whose keys are k-mers and whose values are the probabilities of these kmers
# Output: A randomly chosen k-mer with respect to the values in Probabilities
def WeightedDie(Probabilities):
    kmer = ''     # output variable
    r = random.uniform(0,1)
    s = 0.0
    for kmer in Probabilities:
        if s + Probabilities[kmer] > r:
            return kmer
        s += Probabilities[kmer]
        
        
        
        
        