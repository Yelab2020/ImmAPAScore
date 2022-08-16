# ImmAPAScore
ImmAPAScore is an integrated algorithm estimating the relationship of APA events and tumor immunity to characterize the regulatory landscape of APA events in tumor. We provided all main codes of constructing the ImmAPAScroe algorithm, filtering and ranking the ImmAPAs. The main process of calculating ImmAPAScores for APA events and find ImmAPAs was described below.

**Construction of ImmAPAScore to identify immune-related APA in human
cancer**

In order to identify APA events that have potential effects on
immune-related pathways, and referring to the previous correlation
analysis method of immune-related features, we developed an
algorithm named ImmAPAScore based on the APA event profile and gene
expression of tumor samples to measure whether APA is related to tumor
immunity refer to previous association analysis method of
immune related features. All genes were ranked based on their
correlation with a specific APA event. The ranked gene list for the
specific APA event will be performed Gene Set Enrichment Analysis
(GSEA) for each immune-related pathway. Each APA event in each
cancer type was estimated and obtained ImmAPAScore which assesses the
correlation of APA and each immune-related pathway.

For each APA event, we ranked all genes based on the correlation of
their expression with the PDUI of APA event in each cancer type. The
expression of gene $x$, APA events $y$ and tumor purity scores $p$
across all samples in single cancer type were defined as
$G(x)\  = \ (g_{1},\ g_{2},\ \ldots,\ g_{x},\ \ldots,\ g_{n})$,
$A(y)\  = \ (a_{1},a_{2},\ \ldots,\ a_{y},\ \ldots,\ a_{n})$ and
$P = \ (p_{1},\ p_{2},\ \ldots,\ p_{z},\ \ldots,\ p_{n})$. We introduced
the tumor purity score as co-variable, and calculated the correlation
(Corr) between gene $x$ and APA event $y$ as follow:

$$Corr(xy)\  = \frac{R_{\text{GA}} - R_{\text{GP\ }}*\text{\ R}_{\text{AP}}}{\sqrt{1 - R_{\text{GP}}^{2}}\ \ *\sqrt{1 - R_{\text{AP}}^{2}}}$$

where $R_{\text{GA}}$, $R_{\text{GP\ }}$, and $\text{\ R}_{\text{AP}}$
are the correlation coefficients between the expression of gene $x$ and
the PDUI of APA events $y$, the expression of gene $x$ and tumor purity,
and the PDUI of APA events $y$ and tumor purity, respectively. In
addition, we obtained the P- value for the Corr, defined as $P(xy)$. For
each APA event-gene pair, we calculated the correlation score (CS) as
follows:

$$CS(xy) = - log10(P(xy))\ *\ sign(Corr(xy))$$

All genes were ranked based on correlation score (CS) for each APA event
and then performed GSEA on 17 immune-related pathways. We got the
enrichment score ($ES(y,i)$) for APA event y and immune-related pathway
$i$, and the P value of the enrichment score ($ES(y,i)$). Finally, we
calculated the ImmAPAScore as follows:

$$ImmAPAScore(y,i) = \left\{ \begin{matrix}
1 - 2P,\ \ if\ ES(y,i) > 0 \\
2p - 1,\ \ if\ ES(y,i) < 0 \\
\end{matrix} \right.\ $$

**Identification of tumor-immunity-related APA events**

Based on the ImmAPAScore of each APA event in all cancer types, we first
excluded APA events that related with immune-related pathway in less
than 10 cancer types and the total number of APA events-ImmPath pairs
less than 15 in all cancer types . And we calculated the Immune
Enrichment (IE) score of each APA event across 31 cancer types based on
the ImmAPAScore. IE score was calculated as follows:

$$IE(y,i) = \sum_{c = 1}^{31}{I(}ImmAPAScore(y,i))*sign(ImmAPAScore(y,i)),$$

$$where\ I(x)\  = \left\{ \begin{matrix}
0,\ \& if\ \ |ImmAPAScore(y,i)|\ \  < \ 0.995 \\
1,\ \& if\ \ |ImmAPAScore(y,i)|\ \  \geq \ 0.995 \\
\end{matrix} \right.\ $$

To get those more significant tumor-immunity-related APA events, we
retained APA event with absolutely IE score greater than 5 in more than
5 immune pathways, and finally 541 APA events was defined as ImmAPA.

**Identification of top-ranked APAs**

**﻿**To focus on more significant tumor immunity related ImmAPAs, we
added $NES\ (y,i)$ (Normalized Enrichment Score) of the GSEA result for
each APA event $y$ and immune-related pathways $i$ (there was no
significant APA event related with Interferons Receptors related
pathway) across all cancer types as *R1*. And we also ranked the PDUI
difference of APA events between tumor and normal tissue ($R2$).
Finally, we ranked all $R1$ of 541 ImmAPAs in each immune-related
pathway and combined with R2 to get the final ranking $R3\ $.
The rank scores were calculated as follows:

$$R1 = \sum_{c = 1}^{31}{I(}ImmAPAScore(y,i))*\ NES(y,i),\ where\ I(x)\  = \left\{ \begin{matrix}
0,\ \ if\ \ |ImmAPAScore(y,i)|\ \  < \ 0.995 \\
1,\ \ if\ \ |ImmAPAScore(y,i)|\ \  \geq \ 0.995 \\
\end{matrix} \right.\ $$

$$R3 = \frac{\sum_{i = 1}^{16}{R1i\  + \ R2\ *\ 16\ }}{32}$$


**Construction of ICB_APASig score by machine learning**

We performed a machine learning-based algorithm to construct ICB_APASig score. Briefly, 
(1) we identified ICB-related ImmAPAs through comparing the differentially expressed ImmAPAs between Tumor Immune Dysfunction and Exclusion (TIDE) prediction-dependent non-responders and responders in TCGA metastatic SKCM, based on Wilcoxon signed-rank test.
(2) We performed univariate Cox regression analysis to identify prognosis relevant DE ImmAPAs by assessing the association of OS and the PDUI level of ICB-related ImmAPAs;
(3) We performed LASSO Cox regression model analysis and selected the optimal combination from ICB-related ImmAPAs in (2). 
(4) We performed the multivariate Cox regression analysis on ICB-related ImmAPAs in (3) to obtain ImmAPAs with significant p-value (p < 0.05) as the final signature, named "ICB_APASig". 
(5) The ICB_APASig score of each sample was constructed based on the PDUI level of $GPNMB$ (PDUI<sub>$GPNMB$</sub> ) and $COL1A1$ (PDUI<sub>$COL1A1$</sub>) and multivariate Cox regression coefficient as the following equations:    

$$ICB\underline{ \ }APASig \ Score = 100*(-0.029 ∗ \ PDUI_{\text{GPNMB\ }} + 0.056 ∗ PDUI_{\text{COL1A1\ }}).$$
