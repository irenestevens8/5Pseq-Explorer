# 5PSeq Explorer (RRID: SCR_027718): Interactive Analysis of Co-translational mRNA Decay and Ribosome Dynamics

**Irene Stevens [ORCID: 0000-0003-3823-1499] <sup>1</sup>, Vicent Pelechano [ORCID: 0000-0002-9415-788X] <sup>1</sup>**


<sup>1</sup> Science for Life Laboratory (SciLifeLab), Department of Microbiology, Tumor and Cell Biology (MTC), Karolinska Institutet, Stockholm, Sweden.     


# **Abstract**
**Background **: Co-translational mRNA decay occurs when 5’to 3’ exonucleases follow the last translating ribosome, generating in vivo ribosome protected fragments. Degradome sequencing (5PSeq) therefore offers unique insights into ribosome dynamics. Despite its potential, resources for systematic analysis of 5′P mRNA decay intermediates and associated features, such as ribosome stalls and collisions, are scarce. 
**Findings **: We introduce 5PSeq Explorer, a web-based platform built from 775 uniformly processed 5PSeq datasets across 23 species suitable for exploring ribosome dynamics in vivo at codon, amino acid, and transcript levels.  
**Conclusions **: By integrating normalised counts, structured metadata, and scalable visualisation tools, 5PSeq Explorer provides a framework for studying the crosstalk between RNA decay and ribosome dynamics. To ensure reproducibility and accessibility, we offer both a public web interface and a Docker-based plug-and-play local version. 
**URL **: https://fivepseq-explorer.serve.scilifelab.se/app/fivepseq-explorer  

# **Data Freeze**
Raw count files for 775 samples, together with metadata and RNA composition is available
[5PSeq Explorer Data Freeze (Version 1)]( https://doi.org/10.48723/zchv-5x22). 

> **Please cite:**  
> Stevens, I. (2025). *5PSeq Explorer Data Freeze (Version 1)* [Data set].  
> Karolinska Institutet.  
> https://doi.org/10.48723/zchv-5x22

# **Source Code**
The web version 5PSeq Explorer source code is available at [Source code]( https://github.com/irenestevens8/5Pseq-Explorer/blob/main/web-app.R). 

# **Docker Containers**
For instructions how download and run precompiled versions of the software, please refer to the 5PSeq Explorer User Guide below. 
[Web 5PSeq Explorer]( https://hub.docker.com/repository/docker/stevensirene/fivepseq-explorer/tags/v1.3/sha256-24bb23231fd947c4574b2db344e0bcfb5726dbf6f504164dae80040a5cec5023).
[Local version 5PSeq Explorer]( https://hub.docker.com/repository/docker/stevensirene/5pseq-explorer-local/tags/v1.0/sha256-db7491fac8cbda52f9bdd6459cb2f178b4227c608370f5d7f464e0d77e61c87a).

# **Scripts for reproduction of manuscript figures**
All manuscript figures may be reproduced from the R scripts within this repository, with the exception of those not prepared in R.


## **Main Figures**

**Figure 3**: [Rscript Markdown]().                    
**Supplementary Figure 4**: [Rscript Markdown]().    

## **Additional data pertaining to the study**

[Codon optimality data]().          
[Codon stability]()

# **Supplementary Tables** 

Table S1: [Public datasets used in this study. ]()        

Table S2: [Metadata ]() 

Table S3: [RNA Composition ]()  


# **Contact** 
For questions regarding data, code and web interface, please contact:
[Irene Stevens](irene.stevens@ki.se)   

For questions regarding the 5PSeq method, please contact:                      
[Vicente Pelechano](vicente.pelechano.garcia@ki.se)   
                         

# **Links to relevant websites** 
[5PSeq Explorer User Guide](https://github.com/irenestevens8/fivepseq/tree/Candida)       
[SciCrunch Registry](https://rrid.site/data/record/nlx_144509-1/SCR_027718/resolver?q=SCR_027718&l=SCR_027718&i=rrid:scr_027718
[Fivepseq package](https://github.com/irenestevens8/fivepseq/tree/Candida)           
[Fivepseq package User Guide](https://fivepseq.readthedocs.io/en/latest/)                             
[Pelechano Lab website](https://pelechanolab.com/)                     

