# TransSite
predicting the change in protein transport-sequence activity upon the presence of genetic variations. (beta)

The presence of single residue mutations in the protein localisation peptides have not discussed or analysed at length. We have created a python code, 

that take into account the change in the biochemical properties (including charge), molecular weight and isoelectric pointm the factors most important

in the protien transportation, and provide a collective change value in these parameters in a simple yet effective way. 

PS: the code is in its beta version.

# Installation

  ## pip

    pip install transsite
    
  ## conda
      
      conda install transsite

  ## Docker
  
  ## source code
  
      git clone https://github.com/AhmedArslan/TransSite.git

# Usage

  ## single protein entry:
  
       transsite <protein> <ref> <position> <mut>
       
       e.g. 
       
       transsite Lactb M 1 S
  ## file:
      tab sepated file with the same format as the single protien entry i.g. <protein> <ref> <position> <mut>
 
 # Output
 
    <protein> <mutation-position> <charge-difference> <molecular-weight-difference> <substitutaion-score> <conservation-score> <prediction>
    
 # contact:
 
    aarslan@stanford.edu
 
 # Cite
 
 if you use this tool, please cite: 
 
  
