# Metabarcoding assessment of anthropisation gradients in Aven d’Orgnac cave

Zélia BONTEMPS(1), Jeanne DORE(1), Eric MAHIEU(2), Marine LEPRETRE(2), Odette BELUCHE(2), Karine LABADIE(2), Wafa ACHOUAK(3,4), Mylène HUGONI(1,5,6), Yves PERRETTE(7), Jean-Jacques DELANNOY(6,7), Thierry HEULIN(3,4) and Yvan MOËNNE-LOCCOZ(1*) 

(1) Univ Lyon, Université Claude Bernard Lyon 1, CNRS, INRAE, VetAgro Sup, UMR 5557 Ecologie Microbienne, F-69622 Villeurbanne, France
(2) Genoscope, Institut de Biologie François Jacob, Commissariat à l'énergie atomique (CEA), Université de Paris-Saclay, F-91057 Evry, France
(3) Aix Marseille Univ, CEA, CNRS, BIAM, LEMiRE, F-13115 Saint Paul-Lez-Durance, France
(4) Aix Marseille Univ, CNRS, FR 3098 ECCOREV, F-13545 Aix-en-Provence, France.
(5) Univ Lyon, Université Claude Bernard Lyon 1, CNRS, INSA de Lyon, UMR Microbiologie Adaptation et Pathogénie, F-69622 Villeurbanne, France
(6) Institut Universitaire de France (IUF), Paris, France
(7) Environnements Dynamiques et Territoires de la Montagne, Université Savoie Mont-Blanc, EDYTEM, F-73370 Le Bourget-du-Lac, France

[DOI]()

## Abstract
Background  
Many caves worldwide are open to the public and tourism-related anthropization can alter the cave microbiome, but the spatial distribution of these effects is poorly understood. Here, we tested the hypothesis that cave microbial diversity would change according to anthropization features. To this end, we selected a double anthropization gradient within a French cave (Aven d’Orgnac), by considering a longitudinal (tourist, trek and speleologist zones) and a transversal gradient (central and lateral positions with regard to the path). The microbial community on cave surfaces was assessed by Illumina MiSeq sequencing. 
Results  
Results showed that microbial biodiversity was lower in the tourist zone than in less anthropized zones, both for bacteria and microeukaryotes, but the difference between central and lateral positions was not significant. Microbial community structure was statistically different in the tourist zone than in the other zones, and the difference between central and lateral locations was also significant but of less magnitude. In parallel, particularities in community composition were evidenced, as Alphaproteobacteria and Bacteroidia classes were favored and Sordariomycetes class counterselected in the tourist zone. In addition, 15.2% of bacterial OTUs and 4.8% of microeukaryotic OTUs were observed only in the tourist zone, whereas 1.6% of bacterial OTUs and 5.7% of microeukaryotic OTUs were found exclusively in the speleologist zone. 
Conclusion  
Overall, anthropization affected microbial diversity, and its effects were of higher significance along the longitudinal than the transversal gradient, and for bacteria than for microeukaryotes.   

  
Keywords: Aven d’Orgnac; Cave anthropization; Unbalanced microbiota; Metabarcoding


Code used for the Orgnac project.
## Scripts
* Non-Metric Dimensional Scaling (NMDS)
* Diversity index
* Statistics

## R version
R 4.0.3

## Library 
* library(phyloseq)
* library(minpack.lm)
* library(Hmisc)
* library(vegan)
* library(stats4)
* library(ade4)
* library(picante)
* library(iCAMP)

