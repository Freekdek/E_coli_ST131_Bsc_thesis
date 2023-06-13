# Whole genome typing methods for differentiating between epidemiologically related and unrelated E. coli ST131 isolates.
![Whole genome typing tools](/git_icon.png "Whole genome typing tools icons")

### Introduction 

Current typing methods have difficulty to distinguish epidemiologically related and unrelated E. coli 
ST131 samples. In this study, whole genome typing methods Ridom SeqSphere, chewBBACA, 
PopPUNK and Snippy were compared for their ability to distinguish related from unrelated E. coli
ST131 samples and form clusters based on cutoffs determined in this study.

### Methods 

Samples collected from the same patient within 3 months were defined as related and unrelated 
when collected from different institutions. Samples were sequenced with Illumina and typed using 
the stable E. coli wgMLST scheme in Ridom Seqsphere, a specific E. coli ST131 wgMLST scheme in
Ridom SeqSphere and chewBBACA, PopPUNK, and Snippy. The medians of genetic distances were 
compared using the Mann-Whitney U test and the correlation between the typing methods was
calculated using the Spearman Rank correlation. The cutoff threshold and percentile 
misclassifications were determined at three cutoffs: the lowest cutoff where all epidemiologically 
related were genetically classified as related (cutoff 1), the highest cutoff where all epidemiologically 
unrelated were classified as genetically unrelated (cutoff 2) and the maximum accuracy (cutoff 3).
The percentile misclassifications as either related or unrelated were calculated for cutoffs 1 and 2 
respectively, and both misclassifications for cutoff 3. The cutoffs were then applied to a real-world 
dataset and clusters were evaluated using the Simpson's diversity index and Silhouette score.

### Results 

All typing methods showed a significant genetic distance difference between related and unrelated 
samples. The percentage of misclassifications at cutoff 1 was between 35.9% - 47.1% for the 
wgMLST schemes methods, 18.3% for PopPUNK and 3.3% for Snippy. Cutoffs 2 & 3 showed a lower
percentage of misclassifications, 20.0% and 16.0% for chewBBACA and 8.0% for the remaining typing 
methods. Cutoff 1 had a poor performance compared to cutoff 2, which formed well-defined 
clusters, as shown by a high Simpson's diversity index and Silhouette score greater than 0.85

### Conclusion 
Contrary to expectation, the typing methods are capable of distinguishing strictly defined related 
from unrelated samples. Cutoffs 2 and 3 seem to be best suited for clustering and the typing 
methods have a better clustering performance than using a stable E. coli scheme in Ridom 
SeqSphere. It is recommended to use Snippy, which had the best performance in misclassifications 
and the highest Silhouette score at cutoffs 2 & 3. Although the clusters require to be validated based 
on epidemiological data. For future studies, the adjusted Rand index can be used to evaluate the 
clusters based on epidemiological data.

See <a href="https://github.com/Freekdek/E_coli_ST131_Bsc_thesis/blob/main/Afstudeerverslag Freek de Kreek Microvida.pdf">here</a> for the thesis.


## Poster presentation (clustering excluded)
<object data="https://github.com/Freekdek/E_coli_ST131_Bsc_thesis/blob/main/Poster_Presentation_Freek_de_Kreek.pdf" width="700px" height="700px">
    <embed src="https://github.com/Freekdek/E_coli_ST131_Bsc_thesis/blob/main/Poster_Presentation_Freek_de_Kreek.pdf.pdf">
        <p>This browser does not support PDFs. Please download the PDF to view it: <a href="https://github.com/Freekdek/E_coli_ST131_Bsc_thesis/blob/main/Poster_Presentation_Freek_de_Kreek.pdf">Download PDF</a>.</p>
    </embed>
</object>)

## Data availability
For the Python scripts used for analysis in this study, please refer to <a href="https://github.com/Freekdek/E_coli_ST131_Bsc_thesis/blob/main/scripts/">this link</a>.

Results can also be viewed <a href="https://github.com/Freekdek/E_coli_ST131_Bsc_thesis/blob/main/results/">here</a>.

## Flowchart
![Whole genome typing tools](/flowchart_eind.drawio.svg.png "Flowchart E coli ST131")
