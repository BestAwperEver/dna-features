import numpy as np
import collections
import os
from scipy import stats

s1 = "ACTTATTTACCAAGCATTGGAGGAATATCGTAGGTAAAAATGCCTATTGGATCCAAAGAG\
AGGCCAACATTTTTTGAAATTTTTAAGACACGCTGCAACAAAGCAGGTGGCGCGAGCTTC\
TGAAACTAGGCGGCAGAGGCGGAGCCGCTGTGGCACTGCTGCGCCT\
CTGCTGCGCCTCGGGTGTCTTTTGCGGCGGTGGGTCGCCGCCGGGAGAAGCGTGAGGGGA\
CAGATTTGTGACCGGCGCGGTTTTTGTCAGCTTACTCCGGCCAAAAAAGAACTGCACCTC\
TGGAGCGGATTTAGGACCAATAAGTCTTAATTGGTTTGAAGAACTTTCTTCAGAAGCTCCACCCTATA\
ATTCTGAACCTGCAGAAGAATCTGAACATAAAAACAACAATTACGAACCAAACCTATTTA\
AAACTCCACAAAGGAAACCATCTTATAATCAGCTGGCTTCAACTCCAATAATATTCAAAG\
AGCAAGGGCTGACTCTGCCGCTGTACCAATCTCCTGTAAAAGAATTAGATAAATTCAAAT\
TAGACTTAGGAAGGAATGTTCCCAATAGTAGACATAAAAGTCTTCGCACAGTGAAAACTAAAATGGATC\
AAGCAGATGATGTTTCCTGTCCACTTCTAAATTCTTGTCTTAGTGAAAG\
TCCTGTTGTTCTACAATGTACACATGTAACACCACAAAGAGATAAGTCAG\
TGGTATGTGGGAGTTTGTTTCATACACCAAAGTTTGTGAAG\
GGTCGTCAGACACCAAAACATATTTCTGAAAGTCTAGGAGCTGAGGTGGATCCTGATATG\
TCTTGGTCAAGTTCTTTAGCTACACCACCCACCCTTAGTTCTACTGTGCTCATAG\
TCAGAAATGAAGAAGCATCTGAAACTGTATTTCCTCATGATACTACTGCT\
AATGTGAAAAGCTATTTTTCCAATCATGATGAAAGTCTGAAGAAAAATGATAGATTTATC\
GCTTCTGTGACAGACAGTGAAAACACAAATCAAAGAGAAGCTGCAAGTCATG\
GATTTGGAAAAACATCAGGGAATTCATTTAAAGTAAATAGCTGCAAAGACCACATTGGAA\
AGTCAATGCCAAATGTCCTAGAAGATGAAGTATATGAAACAGTTGTAGATACCTCTGAAG\
AAGATAGTTTTTCATTATGTTTTTCTAAATGTAGAACAAAAAATCTACAAAAAGTAAGAA\
CTAGCAAGACTAGGAAAAAAATTTTCCATGAAGCAAACGCTGATGAATGTGAAAAATCTA\
AAAACCAAGTGAAAGAAAAATACTCATTTGTATCTGAAGTGGAACCAAATGATACTGATC\
CATTAGATTCAAATGTAGCAAATCAGAAGCCCTTTGAGAGTGGAAGTGACAAAATCTCCA\
AGGAAGTTGTACCGTCTTTGGCCTGTGAATGGTCTCAACTAACCCTTTCAGGTCTAAATG\
GAGCCCAGATGGAGAAAATACCCCTATTGCATATTTCTTCATGTGACCAAAATATTTCAG\
AAAAAGACCTATTAGACACAGAGAACAAAAGAAAGAAAGATTTTCTTACTTCAGAGAATT\
CTTTGCCACGTATTTCTAGCCTACCAAAATCAGAGAAGCCATTAAATGAGGAAACAGTGG\
TAAATAAGAGAGATGAAGAGCAGCATCTTGAATCTCATACAGACTGCATTCTTGCAGTAA\
AGCAGGCAATATCTGGAACTTCTCCAGTGGCTTCTTCATTTCAGGGTATCAAAAAGTCTA\
TATTCAGAATAAGAGAATCACCTAAAGAGACTTTCAATGCAAGTTTTTCAGGTCATATGA\
CTGATCCAAACTTTAAAAAAGAAACTGAAGCCTCTGAAAGTGGACTGGAAATACATACTG\
TTTGCTCACAGAAGGAGGACTCCTTATGTCCAAATTTAATTGATAATGGAAGCTGGCCAG\
CCACCACCACACAGAATTCTGTAGCTTTGAAGAATGCAGGTTTAATATCCACTTTGAAAA\
AGAAAACAAATAAGTTTATTTATGCTATACATGATGAAACATCTTATAAAGGAAAAAAAA\
TACCGAAAGACCAAAAATCAGAACTAATTAACTGTTCAGCCCAGTTTGAAGCAAATGCTT\
TTGAAGCACCACTTACATTTGCAAATGCTGATTCAG\
GTTTATTGCATTCTTCTGTGAAAAGAAGCTGTTCACAGAATGATTCTGAAGAACCAACTT\
TGTCCTTAACTAGCTCTTTTGGGACAATTCTGAGGAAATGTTCTAGAAATGAAACATGTT\
CTAATAATACAGTAATCTCTCAGGATCTTGATTATAAAGAAGCAAAATGTAATAAGGAAA\
AACTACAGTTATTTATTACCCCAGAAGCTGATTCTCTGTCATGCCTGCAGGAAGGACAGT\
GTGAAAATGATCCAAAAAGCAAAAAAGTTTCAGATATAAAAGAAGAGGTCTTGGCTGCAG\
CATGTCACCCAGTACAACATTCAAAAGTGGAATACAGTGATACTGACTTTCAATCCCAGA\
AAAGTCTTTTATATGATCATGAAAATGCCAGCACTCTTATTTTAACTCCTACTTCCAAGG\
ATGTTCTGTCAAACCTAGTCATGATTTCTAGAGGCAAAGAATCATACAAAATGTCAGACA\
AGCTCAAAGGTAACAATTATGAATCTGATGTTGAATTAACCAAAAATATTCCCATGGAAA\
AGAATCAAGATGTATGTGCTTTAAATGAAAATTATAAAAACGTTGAGCTGTTGCCACCTG\
AAAAATACATGAGAGTAGCATCACCTTCAAGAAAGGTACAATTCAACCAAAACACAAATC\
TAAGAGTAATCCAAAAAAATCAAGAAGAAACTACTTCAATTTCAAAAATAACTGTCAATC\
CAGACTCTGAAGAACTTTTCTCAGACAATGAGAATAATTTTGTCTTCCAAGTAGCTAATG\
AAAGGAATAATCTTGCTTTAGGAAATACTAAGGAACTTCATGAAACAGACTTGACTTGTG\
TAAACGAACCCATTTTCAAGAACTCTACCATGGTTTTATATGGAGACACAGGTGATAAAC\
AAGCAACCCAAGTGTCAATTAAAAAAGATTTGGTTTATGTTCTTGCAGAGGAGAACAAAA\
ATAGTGTAAAGCAGCATATAAAAATGACTCTAGGTCAAGATTTAAAATCGGACATCTCCT\
TGAATATAGATAAAATACCAGAAAAAAATAATGATTACATGAACAAATGGGCAGGACTCT\
TAGGTCCAATTTCAAATCACAGTTTTGGAGGTAGCTTCAGAACAGCTTCAAATAAGGAAA\
TCAAGCTCTCTGAACATAACATTAAGAAGAGCAAAATGTTCTTCAAAGATATTGAAGAAC\
AATATCCTACTAGTTTAGCTTGTGTTGAAATTGTAAATACCTTGGCATTAGATAATCAAA\
AGAAACTGAGCAAGCCTCAGTCAATTAATACTGTATCTGCACATTTACAGAGTAGTGTAG\
TTGTTTCTGATTGTAAAAATAGTCATATAACCCCTCAGATGTTATTTTCCAAGCAGGATT\
TTAATTCAAACCATAATTTAACACCTAGCCAAAAGGCAGAAATTACAGAACTTTCTACTA\
TATTAGAAGAATCAGGAAGTCAGTTTGAATTTACTCAGTTTAGAAAACCAAGCTACATAT\
TGCAGAAGAGTACATTTGAAGTGCCTGAAAACCAGATGACTATCTTAAAGACCACTTCTG\
AGGAATGCAGAGATGCTGATCTTCATGTCATAATGAATGCCCCATCGATTGGTCAGGTAG\
ACAGCAGCAAGCAATTTGAAGGTACAGTTGAAATTAAACGGAAGTTTGCTGGCCTGTTGA\
AAAATGACTGTAACAAAAGTGCTTCTGGTTATTTAACAGATGAAAATGAAGTGGGGTTTA\
GGGGCTTTTATTCTGCTCATGGCACAAAACTGAATGTTTCTACTGAAGCTCTGCAAAAAG\
CTGTGAAACTGTTTAGTGATATTGAGAATATTAGTGAGGAAACTTCTGCAGAGGTACATC\
CAATAAGTTTATCTTCAAGTAAATGTCATGATTCTGTTGTTTCAATGTTTAAGATAGAAA\
ATCATAATGATAAAACTGTAAGTGAAAAAAATAATAAATGCCAACTGATATTACAAAATA\
ATATTGAAATGACTACTGGCACTTTTGTTGAAGAAATTACTGAAAATTACAAGAGAAATA\
CTGAAAATGAAGATAACAAATATACTGCTGCCAGTAGAAATTCTCATAACTTAGAATTTG\
ATGGCAGTGATTCAAGTAAAAATGATACTGTTTGTATTCATAAAGATGAAACGGACTTGC\
TATTTACTGATCAGCACAACATATGTCTTAAATTATCTGGCCAGTTTATGAAGGAGGGAA\
ACACTCAGATTAAAGAAGATTTGTCAGATTTAACTTTTTTGGAAGTTGCGAAAGCTCAAG\
AAGCATGTCATGGTAATACTTCAAATAAAGAACAGTTAACTGCTACTAAAACGGAGCAAA\
ATATAAAAGATTTTGAGACTTCTGATACATTTTTTCAGACTGCAAGTGGGAAAAATATTA\
GTGTCGCCAAAGAGTCATTTAATAAAATTGTAAATTTCTTTGATCAGAAACCAGAAGAAT\
TGCATAACTTTTCCTTAAATTCTGAATTACATTCTGACATAAGAAAGAACAAAATGGACA\
TTCTAAGTTATGAGGAAACAGACATAGTTAAACACAAAATACTGAAAGAAAGTGTCCCAG\
TTGGTACTGGAAATCAACTAGTGACCTTCCAGGGACAACCCGAACGTGATGAAAAGATCA\
AAGAACCTACTCTATTGGGTTTTCATACAGCTAGCGGGAAAAAAGTTAAAATTGCAAAGG\
AATCTTTGGACAAAGTGAAAAACCTTTTTGATGAAAAAGAGCAAGGTACTAGTGAAATCA\
CCAGTTTTAGCCATCAATGGGCAAAGACCCTAAAGTACAGAGAGGCCTGTAAAGACCTTG\
AATTAGCATGTGAGACCATTGAGATCACAGCTGCCCCAAAGTGTAAAGAAATGCAGAATT\
CTCTCAATAATGATAAAAACCTTGTTTCTATTGAGACTGTGGTGCCACCTAAGCTCTTAA\
GTGATAATTTATGTAGACAAACTGAAAATCTCAAAACATCAAAAAGTATCTTTTTGAAAG\
TTAAAGTACATGAAAATGTAGAAAAAGAAACAGCAAAAAGTCCTGCAACTTGTTACACAA\
ATCAGTCCCCTTATTCAGTCATTGAAAATTCAGCCTTAGCTTTTTACACAAGTTGTAGTA\
GAAAAACTTCTGTGAGTCAGACTTCATTACTTGAAGCAAAAAAATGGCTTAGAGAAGGAA\
TATTTGATGGTCAACCAGAAAGAATAAATACTGCAGATTATGTAGGAAATTATTTGTATG\
AAAATAATTCAAACAGTACTATAGCTGAAAATGACAAAAATCATCTCTCCGAAAAACAAG\
ATACTTATTTAAGTAACAGTAGCATGTCTAACAGCTATTCCTACCATTCTGATGAGGTAT\
ATAATGATTCAGGATATCTCTCAAAAAATAAACTTGATTCTGGTATTGAGCCAGTATTGA\
AGAATGTTGAAGATCAAAAAAACACTAGTTTTTCCAAAGTAATATCCAATGTAAAAGATG\
CAAATGCATACCCACAAACTGTAAATGAAGATATTTGCGTTGAGGAACTTGTGACTAGCT\
CTTCACCCTGCAAAAATAAAAATGCAGCCATTAAATTGTCCATATCTAATAGTAATAATT\
TTGAGGTAGGGCCACCTGCATTTAGGATAGCCAGTGGTAAAATCGTTTGTGTTTCACATG\
AAACAATTAAAAAAGTGAAAGACATATTTACAGACAGTTTCAGTAAAGTAATTAAGGAAA\
ACAACGAGAATAAATCAAAAATTTGCCAAACGAAAATTATGGCAGGTTGTTACGAGGCAT\
TGGATGATTCAGAGGATATTCTTCATAACTCTCTAGATAATGATGAATGTAGCACGCATT\
CACATAAGGTTTTTGCTGACATTCAGAGTGAAGAAATTTTACAACATAACCAAAATATGT\
CTGGATTGGAGAAAGTTTCTAAAATATCACCTTGTGATGTTAGTTTGGAAACTTCAGATA\
TATGTAAATGTAGTATAGGGAAGCTTCATAAGTCAGTCTCATCTGCAAATACTTGTGGGA\
TTTTTAGCACAGCAAGTGGAAAATCTGTCCAGGTATCAGATGCTTCATTACAAAACGCAA\
GACAAGTGTTTTCTGAAATAGAAGATAGTACCAAGCAAGTCTTTTCCAAAGTATTGTTTA\
AAAGTAACGAACATTCAGACCAGCTCACAAGAGAAGAAAATACTGCTATACGTACTCCAG\
AACATTTAATATCCCAAAAAGGCTTTTCATATAATGTGGTAAATTCATCTGCTTTCTCTG\
GATTTAGTACAGCAAGTGGAAAGCAAGTTTCCATTTTAGAAAGTTCCTTACACAAAGTTA\
AGGGAGTGTTAGAGGAATTTGATTTAATCAGAACTGAGCATAGTCTTCACTATTCACCTA\
CGTCTAGACAAAATGTATCAAAAATACTTCCTCGTGTTGATAAGAGAAACCCAGAGCACT\
GTGTAAACTCAGAAATGGAAAAAACCTGCAGTAAAGAATTTAAATTATCAAATAACTTAA\
ATGTTGAAGGTGGTTCTTCAGAAAATAATCACTCTATTAAAGTTTCTCCATATCTCTCTC\
AATTTCAACAAGACAAACAACAGTTGGTATTAGGAACCAAAGTGTCACTTGTTGAGAACA\
TTCATGTTTTGGGAAAAGAACAGGCTTCACCTAAAAACGTAAAAATGGAAATTGGTAAAA\
CTGAAACTTTTTCTGATGTTCCTGTGAAAACAAATATAGAAGTTTGTTCTACTTACTCCA\
AAGATTCAGAAAACTACTTTGAAACAGAAGCAGTAGAAATTGCTAAAGCTTTTATGGAAG\
ATGATGAACTGACAGATTCTAAACTGCCAAGTCATGCCACACATTCTCTTTTTACATGTC\
CCGAAAATGAGGAAATGGTTTTGTCAAATTCAAGAATTGGAAAAAGAAGAGGAGAGCCCC\
TTATCTTAGTGGGAGAACCCTCAATCAAAAGAAACTTATTAAATGAATTTGACAGGATAATAGAAAATCAAG\
AAAAATCCTTAAAGGCTTCAAAAAGCACTCCAGATG\
GCACAATAAAAGATCGAAGATTGTTTATGCATCATGTTTCTTTAGAGCCGATTACCTGTG\
TACCCTTTCG\
CACAACTAAGGAACGTCAAGAGATACAGAATCCAAATTTTACCGCACCTGGTCAAGAATT\
TCTGTCTAAATCTCATTTGTATGAACATCTGACTTTGGAAAAATCTTCAAGCAATTTAGC\
AGTTTCAGGACATCCATTTTATCAAGTTTCTGCTACAAGAAATGAAAAAATGAGACACTT\
GATTACTACAGGCAGACCAACCAAAGTCTTTGTTCCACCTTTTAAAACTAAATCACATTT\
TCACAGAGTTGAACAGTGTGTTAGGAATATTAACTTGGAGGAAAACAGACAAAAGCAAAA\
CATTGATGGACATGGCTCTGATGATAGTAAAAATAAGATTAATGACAATGAGATTCATCA\
GTTTAACAAAAACAACTCCAATCAAGCAGTAGCTGTAACTTTCACAAAGTGTGAAGAAGA\
ACCTTTAG".lower()

s2 = "gtacctctgtctttttttttttgtaaatagtacatatagttttatagatgacgattcctt\
ctgtgtttttttctgctttttaaaatcttcatatcttatatttaatcttaggcatcatct\
gtatacatgattgtttaggtctttaattaccagtgtttagaatcaggtcactcaaacatg\
gtagataagtttgcatagtttgtgtatatccatcactcttgagacagttttattttaagt\
tccggggtacatgtgcaggatgtgcaggtttgttacataagtaaacgtatgccatgttgg\
tttgctgcacctgtcaacccttcacctgagtattaagcccagcatgcattagctattttt\
cctggtgctctccttccccccacacacccccacctcctgacagaccctagtgtgtgttgt\
tcccctccctgtgtccgtgtgttctcattgttcagctcccacttatgagtgagaacatgt\
gatgtttagttttctgttcctgcattagtttgcttaggataatggcttccagctccatct\
gtgtccctgcaaaggacgtgatcttgttcctttttatggctacatggtattccatggtgt\
atagttccacattttatttatccagtctatcattgatgggcatttgggttgattccatgt\
ctgtgctattgtgaatagtgctgcagtgaatgtacaggtggatgtatctttataatacaa\
tgatttatcttcctttgggtatataccccgtaatgggattgctgagtcagatggtatttt\
tggttctaggtctttgaggaattgccacactgtcttccacaacggttgaactaatttaca\
ttccagccaacaacttgagacagtttttgactcataaacattcagagcttggctagctaa\
ttcctgctttaatttaaaaagtgtttattatatgcaaattggacaactcatataaatatg\
tggtgctacttactatgtattttctctaaagcatgttaaaaaaataggctagatatagtg\
gctcatgcctgtaatcttagcactttgggaggctaaggcaggaggatcacttatggtcag\
gagtttaagaacaccctgggcaacatagcgagaccccatctctacaaaaaatttaaaata\
cccaggcatggtggcatgcttctgatgttgtagctactcaggatgctcagacaggaggat\
cacttgagcccaagtgactgaggctgcagtgaaccaaaattgtaccagtgcactccagcc\
tgggccacaaaatgagaccttgtccctgaaaaaaaaaaaagaaaaaaaaaatttaaatag\
aggaaatactagctaagtttaatgtaggccagttctaaaataatgatttattgctgctgt\
tgttacataattttcttaaatattttaaagattgcatactgttactgctctatttctgca\
tctccgtggtgtaactctgtcctctttgttgttgcaacagttcacttagcaactaaactg\
tatgtttacaaagtgattttatctccctatgagaagactttagtgaatagctcagtgaat\
agtagagttggtgagaccacagtacagaactgtttgaagtttgggttaaatttttagagg\
aaaatgtttgatactatgcatatcatagttaaagccaatgaaaaagctaatataggccag\
gcgcagtggctcacgcctataatcccagcactttgggaggccaaggcaggcagatcactt\
aaggtcaagagttcaagaccagcctggccaacatggtaaaaccccatctctatgaaaaaa\
aacaaaaattatccagatgtggtggcatgtgcctgtaatcccagctactcgggacgctaa\
ggcaggagaatcacttgaacctgggagatggaggttgcaatgagctgagatcacgccact\
gcactccagcctgggtgacagaacgagactccatctcaaaaaaaaaaaaaaaaaaagcta\
atacatgtgatcactgatgaaatgcaattaagaactggttagtagaaaattcagagggtc\
aagaaatttaacagagcagttgaactcatttgcctttatcgttgagattagatcatcttt\
caggctgttagtatatggaccctgtttttaaaaattgtggttttgtttttttcaatgtga\
aagaattaagaaaattgttacttttctaattccttttctgtgccttgcttttctgttcac\
accagtattaacagcaatgaaattttttcaattttattttccaataaaaattactttgag\
ttttttttatggtagctagctacttccttgacctagatactaattttgattgagttggta\
actattattaaaaaaacaacttaggtctaatttatcttgagctaaaaaatgtaataactg\
aaaaatagagcatatttaggattctttctgctttaaatttgacattcagttattttcatg\
taatttgtgttttgagcactaccttttaattaatttatttatttttattttttagagact\
gtctcattctgttacctagtctggagtgcactagtgtgatctcagctcaccgtagcctca\
ccctcctgggctcaagcagtccttgcacctcaccctcctgagtaactggcaccacaggca\
tacaccaccacacccagctaatttttatttttcatagagtcatggtctcactatgttgcc\
caggctagtctcgaactcctgggctcaagcagtcttcctgcctcagcctcccaaaagtgc\
tgagattacaggcatgagccactgtgcccaaacactacctttttaacttagtgaaaaata\
tttagtgaatgtgattgatggtactttaattttgtcactttgtgtttttatgtttag\
gtaagtgttcatttttacctttcgtgttgccaatcactatttttaaagtgtttattcagt\
agacttggtatgctaacaattaagagtgttataaactatgtcttttcagccatttttgtg\
tagtcagtttgggggagtatggtttgatatacagatacacagattcagtattcgtataca\
gatttgatatcttggtatacagattcgatatctctgaatctgtataccaagaaatcatgt\
tttaagggtctcaatatattttcaaaaagattattagtataataattgagaaattactgt\
taaaaagttttgagtttctctagaaaatttgaaactcttaacaaaacctgcataatacta\
acttaactgttttcatatacatagcaagttcagactctgacttatatgaactttaaaagt\
tggtttccgggaggccgaggcgggcggatcacgaggtcaggagatcgagaccatcccggc\
taaaacggtgaaaccccgtctctactaaaaaaatacaaaaaattagccgggcgtagtggc\
gggcgcctgtagtcccagctacttgggaggctgaggcaggagaatggcgtgaacctggga\
ggcggagcttgcagtgagccgagatcccgccactgcactccagcctgggcgacagagcga\
gactccgtctcaaaaaaaaaaaaaaagttggtttccgattataccatttactgggtaata\
tatactacttagttacactacttacatagcttcagtttccttatctataaaatgcaaata\
acacctcccatgagggctgggcgtggcgctcatgcctgtaatcccagcactttgggaggc\
cgaggtgggtggatcacctgaggtcaggagtttgagaccagcctgaccaacatggtgaaa\
ccccatctttactaaaaatacaaaaaattagccaagcgtggtggcgcgcacctataatcc\
caactactccagaagctgaggcaggagaatcacctgaacctgggaggtggagggtgcagt\
gagctgacatcacaccactgctctccagcctgggcaacagagcgagactgtctcaaaaaa\
aaaaaaaaaaaagtgtatttaaagcacttagcagtgaacttgacatatagtaggcagaga\
gcattcagtaagtgttggcttgctccctttttttcatttaggaagtgatctaaaaacagt\
attgttagtaaatggtatcttgatcttaatgttatgtggactattttaacttccctttta\
aatgtatatatatctaacaacttagttcaactacagtcatgtgtcatttgacagggatat\
atgttctgagaaatagattgttagatttcatcattgtgggaacatcatagagtatactta\
cacaaacctaggtggtatagcctactatatacctaggctgtatggtatagcttattgctc\
ctaggctgcaaacctatacagcatgttactgtcctgaatactctaggcagttttaacaca\
gtggcaagcatttgtgtatgtgaacatagaaaaggtacagtaaaaatacggtattaaaat\
cttatggggctgggctcagtggctcatgcctgtaatcccagcactttgggaggctgaggc\
aggcggatcacctgaggtcaggagtttgagaccagcctggccaacatggtaaaaccttgt\
ctctactaaaaatataaaaattagctgggcatggtggtggcacacgcctgtaatcccagc\
tactagggaagttgaagcaggagaatcacttgaaccctggaggcagagatttcagtgagc\
caagatcgcaccactgcactcctgcctgggcgacagagcaagactccatctgaaaaaaaa\
aaaaaatcttatgggaccactattaaagtcttataggatgaccattgcatatgtggtcta\
ttgttgaccaaaatgtcattatgtggcaaatgactgcattaggttaaccttatacatacc\
tatattaggtatgtatttggttttgtttttttgtgtgtgtttttttctattagtgtatct\
gactggtaataatcttaaataattgaatctgtttgttagttgcaattaaagcaaatgcca\
aaactccaacatttcagtggataatcttaaataactagttcctttttaaaaaacctataa\
actcataaaaatattttagttattagaactcttcctgtctagaccccatgtattacagag\
agacaccgaagttagtctcctcattcaaaaagtgccttttgcccctaagtcattctggtg\
gatacagatttacttaatcaagtgttgtccaggtcacattcaatataggatttactttat\
ggacaaagtagtacgtttatagtacttaaactatttgctgtcctttagtgtgaaattctg\
aggtatatatgcttaaagatatttgtaattcttttgtggaaaataatggctttatttata\
gcaacccattctgttcttgtgcatactgaagtatattgactttccacctagggaaaaaaa\
aaacaataactcagacttgtaaatgctttcaacggtgttactacttaatttccctcattt\
ctgtaacatataagtgtataacttagtcagcttctggttactggaacagtacaggtcact\
gttaaacaattaaaccacttttataataatctaacacctcctaaagccttgcatggacat\
ttttacttattaaattatacaaatttattccctgtaataaagcatcaaaaagcaaagtac\
ctgttatatattatctcagcatgacatggaaatgcctaccttgaattatggtttaatctt\
accctcttagcctctgtagaatttttaaataagaattgtttctattactagtactttaat\
gtaatttgataattgtaaaaagcctcttaactctaattcaaggacctacataataaatta\
ctccttcagttaatggctgcccccgtgctgaaaaaaaaaaaaaaaaagagagaaaaagtt\
tatttgaagaaattttgttaggccttattgccagtaaacctagagttatatttagtgtca\
gtttttcaaaaagtagcttatctgtggtatctggtagcatctgtttatcctatttaggat\
ttatcctgtttagaccctgttaaatagtggtgttttaaagtggtcaaaacagaacaaaaa\
tgtaattgacattgaagactgactttactctttcaaacattaggtcactatttgttgtaa\
gtatttttgtttaacatttaaagagtcaatactttagctttaaaaaaatggtctatagac\
ttttgagaaataaaactgatattatttgccttaaaaacatatatgaaatatttcttttta\
g\
gtaaaattagctttttatttatatctgttctccctctataggtatggtatataatattct\
gacctcaggtgatccacctgcctctcaaagtgctgggattacagacatgagccactgtgc\
ctaatcaaggacctctttatactcttaaaaattactgaggacctaaaagagcatttgttt\
atgtggaatatatctattgatatttaccatattagaaatgtaaattgattaatgttaaaa\
ttagtaatattatgcgttggtcatttggaagatatgagttcactgagttatgcggatctt\
ccgaaagttgacagttttattatgcagtattaaacaatcactttcattgatgccattacc\
gatcagaaaagtttaagtagtagaaagctgtcaagcttacagagccagatacaagcttcc\
caaaaattctgattttcatctaaaagcttgaatttttccccggcaataagtattgtcact\
tatttttcttgtaggtgacaagcttattttcattcatttttgaaaagatgtctgccgaat\
acccaagtctgaataactatagtttgttggttattctttcaagtaaaaggtatttcatga\
aaaaatagctagtatagctcacaactcaatcatttaagtgtgttttcttgagaaacgcac\
tgaagtatgcaagcataatataccaacagtacaaatatcaacagtgaaaaggacatacat\
aacattttactaataagacagttttgacagcttggattccctaaaatggttgtagatacc\
taacaggattccactgatcatttcttgagaatcattgtcctataatatatacataataat\
ctaaatttacaatatcagtattaactactgacaataaaactactaaggaaaatgtaagaa\
ttgtttgcagtttttgtccttagagtatataggttgagtatccctatctgaaatgcttgg\
gaccaggactatttcagatttcagattttttcagattttgaaatgcttgcatatacaata\
cataatgagatatctggggataggactcaagtctaaacacgaaatttatttaagtttcat\
aaacaccttatacatataacttaaatgtaattttatacaatattttaaataatttttgca\
taagacaatttaaattgtgatccatcacatgaggtcagatgtggaattttctactggcct\
catgttggcactcaaaaagtttcaggtttgtgaccattttggattttcagattagggata\
ctcaacccatatattattaagaatgtttagtcaaaatactgtgttcaaatgtcactcaaa\
ataattcttccggatgtggttaccaatttgataattaggttacattcctttttttccatt\
tgttttcaattttaggatttgtcttttcttatttaattttacatttgaataaataaaaca\
ttacatagttcattcatcagaactacaaaaaggtatacttagagtttttattcacccacc\
tcttgcttactataggtaatcttttttagtgtttttttttcaggattctgtttaataaaa\
ataagcaaatacatgtatatactcattaccctttcttactcaaaagatacagtatataca\
ccattttgcaccttgtttattggttgttgtttacttaagaattatttggagatgactcct\
taatgagtatatagagatcgtcctcattcttttttgtggttacatagtagttgatcatct\
ggctgtgtcagtgtttcctagtttatttaaccaatttccaactagtggacttattgaaga\
tttaattaggttccagttacatactgagaatgaacaatatctaaagcttagcttttaaac\
cttcataagactaaattttaaatttggtatttgcatcagaaattagctaacacctttgag\
ttatgatggttaacatcaactgactaaatttatgctgatttctgttgtatgcttgtactg\
tgagttatttggtgcatagtcattatcaatttgtgaatcaatttattttcatagttaaca\
tttattgagcatctgttacattcactgaaaattgtaaagcctataattgtctcaaatttt\
ttgtgtatttacagtaacatggatattctcttagattttaactaatatgtaatataaaat\
aattgtttcctag\
gtaagacatgtttaaatttttctaaattctaatacagtatgagaaaagtctcgtttttat\
aaatgaacatttctaaaaataatgacactaacgttaagaagttaacacttcccgttttat\
aaaatttataaaatactttggtagtattttatagtgctgttcatatcattattttatttt\
ttaattttatgacagctttgtaaagtagacagattttattctaattttatggatgaagta\
ctaaggttgagaggaattaaggaaattgctccgaatcagttaacaaaaagattgcagata\
ttaaaaatatccttttatctctcctctctaaacctttaaaaaagtactaagatagttttt\
ttaatgtataattcccaaggacaatgatgagaagaaacaacaaaagtttggaagccaaaa\
acataaaggatttagtaagcatgagaaagctaaaacctgacactagagcaaacagagatg\
ctttcccctaaaaaacctgaaaaagattcaaattggcagcaacaggtacttctgaaggtg\
aagtagaaaataggaagattagttgaaattctttttaagaaacatctatatttcctcccc\
cactgcaaataggcggttatccttcttctgccaggaaatcagaaggttgttcttgaaaaa\
gatgaattgagaggattctgaattgaaggtgggctggagggaggggacaccaggcacaat\
tgagggaaagatactaaaatgaaagatcagatacaaatctgtatgtcaagcagtgagacc\
tagctccttcccacacttggttcccaaatgcaggccctctaggcatgagactggaagatt\
tttttttcctagggaatatgcctgacccaatagaaaagaccaaaaaatactgacagttga\
ggatactcagatgaaacagtatagccagtcaccagaccaggaagttaactgttgacatgc\
acagagcttccaggaagctacttagtgcttcacttttaaataagaaaagatagtcaaaga\
taactagtcattggaagaaagctactatgaaacatagtcaccaaagtacaaaatccatag\
cagaaaggaacctagaggaaatcgactatgaaaacttcataaaaacctactaatattctc\
aggtaagaaaagaaaaaatggccgtaaaataagaacaagttgctataaaaagctcttaga\
aattaaaaatatgatagcacaaataaattaactcagtagaaataatggaagaatcatgaa\
agttcccagaatacagaataaaatgaaaaaaggtatgaaaagtcaattctgtggatctat\
catctgaaaatacagagtttgagaaggaaggcacagaagagaaatgaagaaagaaatttt\
aaaataaatacataattttaaaagttctactagtactgaaggacatgagtttccttaatt\
aaaagggcccactgagtgagcacacaagtaaaaatgacccacagtaaggcacatccttgt\
gaatttttagaataatagaggcagacaggaaccttaaattcattagaggaccaagaagtt\
aggtttcaaattgtttcaagccataatagtatgaattctcttattatcaacaatggaatc\
tagaagactgtagatcttatataatacagagaagtgccttcaaaatactgagagaaaatg\
atttccaacctagaatctgaattaagtgtgagggtagacatttttcagatgtgaagtact\
aaaagatctcttgtgcgcttttctcaggaaactaaccaaaacaaatgcatacaccaagaa\
ggaggaaggtataggacttaagaaataagaattcaacatagaagagaggcaaagggagct\
ttcaggatgatattgaagggagatcccagagtagctgtgttgctaagtctagaaaggcag\
ctagactactttggaactgaagaagataagagactttggaagagtttgccttcaagataa\
aaataaagcagtacctgcatgttttaatgtattaggaaacttcttagtaaagatggtgaa\
ttgaggccaggcacagtggcttacacctgtaatccagcacattgggaggctgaggtgggt\
agatcacttgaggccaggagttcgagactagcctggccaacatggtaaaatcccatctct\
actgaaaatacaaaaattagccaggcgtggtggcacacgcctgtaatcccagctactcca\
gaggctgaggcacaagaaccgcttgaacctttgaggtggaggttgtggtaaaattgcacc\
actgcacttcagcctgggtgacagagtgagactctgtctcaaaaaaaaaaaaaaaaaaaa\
aaagatggtgaattgaacatactcatatcctttctttgccttccaaacttttaccaaaac\
atcattgaagaaacttacacacacacaaaaaaaaaacaaggaaaataggaaataacaaag\
taactaaatttctcaaagcatgcagaaggaaactgaatgaaagctggtggtggggacagc\
agagaaccaacgattttacactcaggtctcaaaagactaggaattggtggcttcatttct\
tatctttagaattgggtggtgcagaaggagggagccaaaatggaataagttgaaattatg\
tttaagaagcaatactcgccgggtacggtggctcacatggaggctgaggcgggtgaatca\
cctgaggtcaggaattcgagaccagcctggctaacatggtgaaaccccatctctgttaaa\
aatgcaaaaattagctcggcatggtagcatgcccctgtaatccagctactcaggaggctg\
aggtgggagaactgcttgaacccaggaggtggaggctgcagtgagccaagattgcgccac\
tgcactccagcctggacgacagagcaagaccccacatcaaaaaaaaaaaaaaaaaagcag\
cagcagcaatactcatgaagctgggcaactgtctcctgcccgctctatgaaaagaaccag\
aggcttattctccagagaggatacagtagaaggtgaacacactaggcacagttgaaggca\
gaagcaactacttgaaagcaagaagaagttaatatatgcatattgaatgttgggatctcc\
cctcaccaagcccttttccaccactcagcttccagaacatagacagctaagttttcacta\
gtggaagtttccatttaatcaagctactgtgtagcttgcagtcaacaagttctatctttg\
taccaagtgcttcaaaacagccttttggtccctcactcttaactataaacagacatccaa\
agattatgagacatcagaaaaagcaaaaataaaataaccaaaaaacacattaatgaaaac\
aacttagaagaaacattattcaaggagaagaaaaaatgtttttttaaaaactataatttg\
tgaacagaatgaaaagaggtttatatatatagctaagagtttagatgtgaataaacagta\
agtacatagaaaataagcagattttaaaaattaactcaagagaaagcaaaagttgtaaag\
gaagtacactatttatatactacccattaatggccgggtgtggtggttcacgcctgtaat\
cccagcactttgggaggccgaggcgggtggatcacaaggtcaggagatcgagaccatcct\
ggctaacatggtgaaaccccatctctactaaaaataacaaaacaaaattagccagacgta\
gtggtgggcgcctgtagtcccagctacttgggaggctgaggcaggagaatggcatcaacc\
caggaggcggagctttcagtgagccgagattgcaccactgcactccagcctgggcgagag\
agcgagactccgtctcaaaaaaaaaacaacaaaataaaaaaataaaataaaatatactgc\
ctattaatactacatatactttatactgacttagccgtaatgtaaatgttgaacattgat\
agtgagaggtgaagctggctgggcttctgggtcgtgtggggacttggagaacttttctgt\
ccggctaaaggattgtaaacacaccaatcagcgctctgtgtctagctaaaggtttgtaaa\
cgcaccggtcagcactctgtgtctagctaaaggtttgtaaatgcaccaatcagcactctg\
taaaatagaccaatcagcaggacgtgggcggggccaaataagggaataaaagctggccac\
ctgagccagccccagcagccgctcggctccacttccatgccatggaatctttgttttttc\
actctttgcaatgaatcttgctgctgctcactctttagtgagcactacctttatgagctg\
taacactcaccacgaaggtctgcggcttcactcctgaagtcagcaagaccacgaacccac\
caggaagaagaaacaaccctgtacgtgccatctttgagagctgtaacactcactgggaag\
gtctgcggcttcactcctgaagtcagcaagaccacaaacccaccagaaagaagaaactct\
ggacacatctgaacatcaggaagaacaaactcgggacacactatctttaagaactgtaac\
accatgagggtccacagcttcattcttgaagtcagcaagaccaagaacccaccagaagga\
accaattccggacacagtagaattaaatacgtaatttaggaagatgaaaggcaagagtgt\
gtgtgtagtaaggtagaagctgtgttgacagagctgaattttcattttctgtaggggtac\
ttcaagagaaaaagtcaagaagaaacatgtcacttagacatataaatatgataaaatcat\
ctaaaactgtttaaagtagttgcaaaatcttttctagctgataaatttttaagcctaaaa\
atatcattgaaattattttaatgttacattttattttattttatttatttatttatttat\
tttgatacagagtctcactctgtcgcccaggctggagtacagtggcacgatcttggctca\
ctgcaacctctgcctcctaggttcaagcgactctcctgcttcagcctcccaagtagccgg\
gattacaggcgcgtgccaccatgcccggctaattttttgtatttttactagagaaggggt\
ttcaccgtgttagccaggatggtctggatctcctgacctcgtgatccgcccaccttggcc\
tcccaaggtgctgggattacagacgtgagccactgtaccaggcctaatgttaccttttca\
aaaacacctgattgtggaattgttgaagtcactgagttgtatttctggaatgtgtttttt\
agcaggctgcacatacacatatgtagaaagccaggtgatttttttttcatttcttttttt\
tttatcaaaaacagttgtattaaataagaaaggaaatacgtatttacccgtgtattacct\
taatttatgtgtaaaatgggagaatagtttaatgtatttaacaaacaaacatttgttaaa\
gtacctgctcaaactacctaatatatactatagtgaaagatataaggataaataagtcta\
actcagattgctagcctgggaaccagacatgaaaacaagaattataatgtaatataaatt\
ctagaatagatgtaaaaagtgatctaagaacatagaaaaattatcagctaatcacatgac\
tgctcaatgggaaaagtacttcagacagaatgtaaagaatgcttggttaaagatggcatt\
ccaaatcttggaatttggttgggggacagagggaaacaaaaagaaatggggaggttagga\
ccaaataggaagcttcctgtatgtcatttctgataagttgaaacctaggtaggtgatagg\
ctgtctttggaagtttctaacaagaggaacaaaataagattggtgttttagaagtatacc\
aaagcaaaactgttgcaaggagattagtaaatacaggtcttaacctagcagaggaggtag\
agggtagagaatgattgagatagaaattcagtagatttggccagatagtgataagttgag\
actggcaaattatttccacttagatttaaatagatatcttgagcataacctacaaggcaa\
actccttatactaaaaatattctgaatatttaaaaagaaaggattaaaagatcaatcaat\
agaagtttggggacagaaggtttattcattcttgtgcattagatctcatctagatcacct\
gtttgaagaaatcattccagcaattatcttgtctctctcctgcatggatttttttcctaa\
tagattgttctcatcaccctaagcagttgttgtacatctctcatcttaaaaagaacagcc\
tttcttaagtaatctcaacagtccattttcttctcttaaagcccaactcattagaattgt\
ccctcctcttttcacttatctctttgagtactctcctgaacccagtctagtcagtccttt\
cagtagaactggtccccctgcttacctccctactcctcaatacacagtgaattctcaaca\
aagaagccgggggatccttttaaacataagacagattatgtcatttctttactcagaact\
attccgtggtgtgccatctcagagtagagactaaaagcccttgtcatggtgtacagattc\
ttcatgatctggctgctttgctatttttccagtctgaccttctaatgttccccttgctct\
ccttgctccaggcacacttgtgtctaggccaatcgacatatttgtttgtctgtttccttc\
cactgaaatatacatgcaaaaacaaaattttgttaccgtgttccccagcaaaacaatgtc\
tggcacctggtaggaattcattaaatagttgatggatgggcgaacggataactaaaggaa\
caacttcaagttccaggtatccagggtttggtaaaaggaaatctggggttttcaacaaga\
tatcaagtattaggaagaccacgtatgctgagaaagatgatcacttttggacatgttgag\
tttgaaatgagtgtgaaacatcaaggtacagatgtctgatgctatatgtagtgtaaaatg\
taggaacaaccctaggagaaaaatcgggcatgaggataaaggatattttcattgttaggt\
gataatttaagcaatggaaatgactcacattagcaagggaaagtgtctaaggaagacatc\
cagttttggagactttttttgaggaatcaggaagaggtaaaaccagtaaaagatgaaaga\
ggtacagtgatggtgagaattttaaaagaaggaaaatgtaaactgtcatagctattagga\
aagttgagtagaatgagtttgcgtgcatcccacatgcatctgggaggtcattaacaactt\
tattgagaacagtttctgtagagtagtgggagaaatgagagtttattgagtagagattga\
ggaagtgaaaatagctacattacctattgaagaaggttgactgtggagtgtaacagtgag\
tattagcttgaggcagagataaaggtgagtgagaaaataagagtttcaaaggtaggcaag\
atttttgggctaaataaaaagggcactttaaaaaaggtataaataggtagaagagagaaa\
agggagcgaggtgggataattgaaagaggggatctcctgtggagactgaggtattaggcg\
gagtagagagttcaggtgaagatgtgaaggtgagagaagaggatgggtagacatttccct\
ggtgaaggaggtaaggagtactatgatggaattagaggggacacactgagagggtccaca\
cttgacagactctcttctattatgtgttatgtgaggtagattgtaaagtcaaaggctagc\
cttgaaaaatgtgatattgttttggaatggcaaccatggtgaatacaaaacagttaccag\
aatagtatcaccatgtagcaaatgagggtctgcaacaaaggcatattcctaaatatttat\
atgtgtactagtcaataaacttatatattttctccccattgcag\
gtattgtatgacaatttgtgtgatgaatttttgcctttcagttagatatttccgttgtta\
aataatgtcctgatggttttccccctttggtggtggtaattttaaagccctttttaatgt\
tttagattttctaaatccaaagattaggtttaaattattctaatgtttctttcaaagata\
acttcttgtggacttgttaaaaaaaattagacacacaatctaggactgctgttactggaa\
tatattttctatcatgctactaattttctttttaaaatgtgataaaaatagggccgggcg\
tggtggctcatgcctgtaatcccagaactttgggagactaaggcgggcggatcacctgag\
gtcaggagttcaagaccagcctggccaacatagtgaaaccctgtctctactaaaaataca\
aaataaataaataaataaataaatagctgagcgtggtggcaggcacctgtaatcccagct\
gcttgggaggctgaggcaggagaatcgtttgaacccgggaggcagaggttgcagtgagcc\
gagatcgcgccattgcactccagcctgggcaacaagagtgaaaaactctgtctcaaaaag\
agataaaaatagtaaagatattcatatttatacagctttacaagttgaaacatcctttca\
tttatgaagaattaaaaggggtaccctttttagagaaaaggagagcatgtaaacttcgag\
gaaattgatatgtataattttataaaacagggcttgcgcttttttttttttgagacagag\
tttcgctcttgttgcccaggctggagtgcaatggtgcaacctcggctcaccgcaacctcc\
tcctcccgagttcaagtgattctcctgcctcagcctgctgaatagctgggattacaggca\
tgtgccaccacacctggctacttttgtgttttttttacttttatatattttttttttgtt\
tagtagagacagggtttctccattttggtcaggctggtcttgaactcccgacctcagatg\
atctgcccgcctcagcctcccaaagtgctgggattacaggcgtgagccactgtgcctggc\
caggggttgtgctttttaaatttcaattttatttttgctaagtatttattctttgatag"

path = "files/"

r_f = open(path+"dna.txt", 'r')

s1 = ''
s2 = ''

from_i_to_e = 0;
from_e_to_i = -1;
last = 'e'

for line in r_f:
	l = line.strip()
	if len(l) > 1:
		if l[0].islower():
			s2 += l
			if last == 'e':
				last = 'i'
				from_e_to_i += 1
		else:
			s1 += l
			if last == 'i':
				last = 'e'
				from_i_to_e += 1

from_i_to_e /= len(s2)
from_e_to_i /= len(s1)

r_f.close()

s1 = s1.lower()

nucl = "acgt"

nn1 = {}
nn2 = {}

r = 10**4
Ns = [7]
#probs = [90., 91., 92., 93., 94., 95., 96., 97., 98., 99.]
probs = [85. + i for i in range(15)]
#Ts = [1.5, 1.6, 1.7, 1.8, 1.9, 2.0]
Ts = [2.5]

m = 100
prob_finding = .99
prob_threshold = 1 - np.power(1 - prob_finding, 1/m)

def rec(str, numb):
	if numb == 0:
		return
	numb -= 1
	for n in nucl:
		nn1[str+n] = max(s1.count(str+n), 1)
		nn2[str+n] = max(s2.count(str+n), 1)
		#nn1[str+n] = s1.count(str+n)
		#nn2[str+n] = s2.count(str+n)
		rec(str+n, numb)

def not_substring(s: str, l:list) -> bool:
	for ss in l:
		if len(s) > len(ss[0]):
			if s.find(ss[0]) != -1:
				return False
		elif len(s) < len(ss[0]):
			if ss[0].find(s) != -1:
				l.remove(ss)
				return True
	return True

output_probs = open(path+"output_probs.txt", 'w')

for N in Ns:
	for prob in probs:
		for T in Ts:
			print(N, prob, T)
			#N = 7
			#prob = 95.
			z = stats.norm.ppf(prob/100)
			#T = 2
			
			fileName = path+"features.deep-" + str(N) + ".txt"

			i_f = open(path+"important.features.deep-"+str(N)+".prob-"+str(prob)+".threshold-"+str(T)+".txt", 'w')
			important_features = []

			if os.path.isfile(fileName):
				f = open(fileName,'r')
				ex_sum = 0
				in_sum = 0
				for line in f:
					property = line.strip().split()
					a = property[0]
					p_ex = float(property[1])
					int_ex = z*float(property[2])
					p_in = float(property[3])
					int_in = z*float(property[4])

					if min(p_ex, p_in) > 0:
						if p_ex - int_ex > T * (p_in + int_in) or p_in - int_in > T * (p_ex + int_ex):
							if not_substring(a, important_features):
								print(a+':\tex =',int(p_ex*r),'+-',int(int_ex*r),'\tintr = ',int(p_in*r),'+-',int(int_in*r), sep='')
								important_features.append([a, p_ex, int_ex, p_in, int_in])
				for line in important_features:
					a = line[0]
					p_ex = line[1]
					int_ex = line[2]
					p_in = line[3]
					int_in = line[4]
					ex_sum += p_ex
					in_sum += p_in
					i_f.write(a+' '+str(p_ex)+' '+str(int_ex)+' '+str(p_in)+' '+str(int_in)+'\n')
					#i_f.write(a+':\tex = '+str(p_ex)+'+-'+str(int_ex)+'\tintr = '+str(p_in)+'+-'+str(int_in)+'\n')
				i_f.write("exons prob sum: "+str(ex_sum)+",\tintron prob sum: "+str(in_sum))
				if min(ex_sum, in_sum) > prob_threshold and max(ex_sum, in_sum) < prob_threshold*2:
					output_probs.write("N = "+str(N)+", prob = "+str(prob)+", T = "+str(T)+
						"\texons prob sum: "+str(ex_sum)+",\tintron prob sum: "+str(in_sum)+'\n')
				print("exons prob sum: "+str(ex_sum)+",\tintron prob sum: "+str(in_sum))

				f.close()
			else:
				rec('', N)

				#ex_int = {}
				#in_int = {}
				#important_features = {}

				f = open(fileName, 'w')

				nn1 = collections.OrderedDict(sorted(nn1.items()))
				nn2 = collections.OrderedDict(sorted(nn2.items()))

				ex_sum = 0
				in_sum = 0

				for a in nn1:
					#nn1[a] = int(nn1[a]/(len(s1)-len(a)+1)*(10**N))
					#nn2[a] = int(nn2[a]/(len(s2)-len(a)+1)*(10**N))
					p_ex = nn1[a]/(len(s1)-len(a)+1)
					p_in = nn2[a]/(len(s2)-len(a)+1)
					#nn1[a] = p_ex
					#nn2[a] = p_in
					#if max(nn1[a],nn2[a]) > 100 and abs(nn1[a]-nn2[a])/max(nn1[a],nn2[a]) > 0.5:
					#	print(a+': ex =',nn1[a],'intr = ',nn2[a])
					int_ex = z * np.sqrt(p_ex * (1 - p_ex)) / np.sqrt(len(s1)-len(a));
					int_in = z * np.sqrt(p_in * (1 - p_in)) / np.sqrt(len(s2)-len(a));
					#ex_int[a] = int_ex
					#in_int[a] = int_in
					f.write(a+'\t'+str(p_ex)+'\t'+str(int_ex/z)+'\t'+str(p_in)+'\t'+str(int_in/z)+'\n')

					if min(p_ex, p_in) > 0:
						if p_ex - int_ex > T * (p_in + int_in) or p_in - int_in > T * (p_ex + int_ex):
							if not_substring(a, important_features):
								#ex_sum += p_ex
								#in_sum += p_in
								print(a+':\tex =',int(p_ex*r),'+-',int(int_ex*r),'\tintr = ',int(p_in*r),'+-',int(int_in*r), sep='')
								#f.write(a+': ex = '+str(int(p_ex*r))+'+-'+str(int(int_ex*r))+' intr = '+str(int(p_in*r))+'+-'+str(int(int_in*r))+'\n')
								#i_f.write(a+':\tex = '+str(p_ex)+'+-'+str(int_ex)+'\tintr = '+str(p_in)+'+-'+str(int_in)+'\n')
								important_features.append([a, p_ex, int_ex, p_in, int_in])
								#important_features[a] = {"exon_prob": p_ex, "exon_spread": ex_int, "intron_prob": p_in, "intron_spread": int_in}
					#elif nn2[a] - in_int[a] > T * (nn1[a] + ex_int[a]):
					#	print(a+': ex = ',int(p_ex*r),'+-',int(int_ex*r),' intr = ',int(p_in*r),'+-',int(int_in*r), sep='')
				for line in important_features:
					a = line[0]
					p_ex = line[1]
					int_ex = line[2]
					p_in = line[3]
					int_in = line[4]
					ex_sum += p_ex
					in_sum += p_in
					#i_f.write(a+':\tex = '+str(p_ex)+'+-'+str(int_ex)+'\tintr = '+str(p_in)+'+-'+str(int_in)+'\n')
					i_f.write(a+' '+str(p_ex)+' '+str(int_ex)+' '+str(p_in)+' '+str(int_in)+'\n')
				# i_f.write("exons prob sum: "+str(ex_sum)+",\tintron prob sum: "+str(in_sum))
				if min(ex_sum, in_sum) > prob_threshold:
					output_probs.write("N = "+str(N)+", prob = "+str(prob)+", T = "+str(T)+
					   "\texons prob sum: "+str(ex_sum)+",\tintron prob sum: "+str(in_sum)+'\n')
				print("exons prob sum: "+str(ex_sum)+",\tintron prob sum: "+str(in_sum))

				f.close()

			i_f.close();

output_probs.close()
			# print(int(a*100),int(c*100),int(g*100),int(t*100))

			# 37 16 18 27 ex
			# 32 17 19 30 intr

			#print(nn1)
			#print(nn2)