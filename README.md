# BIMM182
### Q2
To generate random sequences and align them with given parameters, run AlignSequences.py:

`python3 AlignSequences.py <number of pairs> <length of seq> -m <match> -s <mismatch> -d <indel>`

>Ran `python3 AlignSequences.py 500 1000 -m 1 -s -30 -d 0 -o Para1.txt` for first set of parameters. 

>Ran `python3 AlignSequences.py 500 1000 -m 1 -s -30 -d -20 -o Para2.txt` for the second set.

The output files, `Para1.txt` and `Para2.txt` are the collections of lengths of alignments for each random pair. The data is used to generate histogram.

To plot histogram, run `python3 PlotHist.py <length_file.txt> -b <bin size>`, remember to change the filename accordingly. The script reads the output file that contains the alignment lengths.

The first set of parameters:
`python3 PlotHist.py Para1.txt -b 20`
![image](./Histograms/Parameter1.png)

The second set of parameters:
`python3 PlotHist.py Para2.txt -b 5`
![image](./Histograms/Parameter2.png)

The lengths of the optimal local alignments are different. 
For the first set of parameters, the lack of a penalty for introducing gaps allows the alignment to extend by adding gaps, thus increasing the length of the alignment. The resulting alignments might include several gaps that separate islands of matches. Therefore, even in sequences with low similarity,  mismatches can be avoided by adding gaps to improve the score.

With parameter2, the large indel penalty makes the algorithm conservative about introducing gaps. This tends to produce shorter alignments with fewer or no gaps. The algorithm will often terminate an alignment segment instead of introducing a gap or continuing through a mismatch to keep the alignment score high.

Tried Parameter1: 
```
500 pairs length 20, alignment length median/average = 27

500 pairs length 200, alignment length median/average = 271

500 pairs length 600, alignment length median/average = 812

500 pairs length 1000, alignment length median/average = 1350
```
From Excel, `l_p1(n)` is a linear function `y - 1.35x + 0.72`, with `R^2=1`, where x is random sequence length, y is alignment length
![image](./Histograms/Para1_Trend.png)

Tried Parameter2: 
```
500 pairs length 20, alignment length median/average = 4.5

500 pairs length 200, alignment length median/average = 7.2

500 pairs length 600, alignment length median/average = 8.2

500 pairs length 1000, alignment length median/average = 10
```
From Excel, `l_p2(n)` is a linear function `y = 0.005x + 5.21`, with `R^2=0.8939`, where x is random sequence length, y is alignment length
![image](./Histograms/Para2_Trend.png)

#### Extra Credit:
Parameter1: the slope of 1.35 means that the alignment length (y) grows faster than the sequence length (x). This attributes to 0 gap penalty that allows even distant matches to contribute to the same local alignment block

Parameter2: the slope of 0.05 indicates a much slower rate of increase in alignment length with increasing sequence length. The high penalties for both mismatches and indels mean that only the strongest matching segments contribute to the optimal local alignment. As sequence length increases, the addition to alignment length is minimal because most potential alignment extensions would incur penalties that outweigh their benefits.

### Q3
Ran `python3 length_function.py -3 -1 -0.75 -0.5 -0.33 -0.25 0`

Found an abrupt change occurs at -1 to 0
![image](Histograms/Alignment_lengths.png)

Tried `python3 length_function.py -3 -1 -0.75 -0.5 -0.33 -0.25 0`

The abrupt change occurs at around -1.38 penalty.
![image](Histograms/Alignment_lengths_zoom_in.png)

#### Extra Credit:
For high penalties, mismatches and indels are discouraged to keep high alignment score. The algorithm favors very short and high-confidence alignments that have good matchings. The abrupt change occurs when the cost of introducing a mismatch or a gap is outweighed by the benefit of continuing the alignment. It's at this threshold where the algorithm shifts from a conservative to a more permissive behavior.

### Q4
Memory usage for **linear-space** local alignment.

`python3 Linear_locAL.py p4seqs.fasta -m 1 -s -10 -d -1 -a`

```
jiayic@Jiayis-MBP BIMM182_A2 % ps aux | grep locAL.py

jiayic           27653  99.3  0.3 409760336  21072 s020  R+    6:12PM   0:45.14 /Library/Frameworks/Python.framework/Versions/3.12/Resources/Python.app/Contents/MacOS/Python Linear_locAL.py p4seqs.fasta -m 1 -s -10 -d -1 -a

jiayic@Jiayis-MBP BIMM182_A2 % ps -o pid,rss,vsz,comm -p 27653

  PID    RSS      VSZ COMM
27653  11712 409891408 /Library/Frameworks/Python.framework/Versions/3.12/Resour
```
Running result:
```
jiayic@Jiayis-MBP BIMM182_A2 %  python3 Linear_locAL.py p4seqs.fasta -m 1 -s -10 -d -1 -a
Alignment:
TGTCCTACGC-TAG--CC-GC-TA-ATGGT-ATCGAAC--CCTG-TGTTTCCCGCGAGCAATCGTTA--TGAGCGA-TGG---GACCTCATCCCGGATTAGGATACTTCACGCTTACGAACTCTCAGGGAC--AGTTCCG
TGT-CTAC-CATA-TACCTG-GT-CAT-GTGA-C-AACAGCC--ATG-TT--C-CGAG---T-G---GCTGAGCG-CT--TTCGACCTCATCCCGGATTAGGATACTTCACGCTTACGAACTCTCAGGGACAGAGTTCCG
Score of the best local alignment: 56
Length of the best local alignment: 140
```


### Q5
I used Python. I spent about 10h in this homework. I did not discuss with anyone. I looked at Piazza posts and Discussion slides. Hint for Q4 is helpful. I looked up some code online ([Matplotlib reference](https://matplotlib.org/stable/api/_as_gen/matplotlib.pyplot.hist.html)) and used my CSE181 local alignment code (line 26 through 70 in `locAL.py`).

