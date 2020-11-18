# TransSite
predicting the change in protein transport-sequence activity upon the presence of genetic variations. (beta)

The presence of single residue mutations in the protein localisation peptides have not discussed or analysed at length. We have created a python code, 

that take into account the change in the biochemical properties (including charge), molecular weight and isoelectric pointm the factors most important

in the protien transportation, and provide a collective change value in these parameters in a simple yet effective way. 

PS: the code is in its beta version.

# installation

  ## pip

    pip install transsite
    
  ## conda
      
    conda install transsite

  ## docker [TODO]
  
  ## source code
  
    git clone https://github.com/AhmedArslan/TransSite.git

# usage

  ## single protein:
  
    transsite <protein> <ref> <position> <mut>
       
     e.g. 
       
    transsite Lactb M 1 S
  ## file:
    tab sepated file with the same format as the single protien entry i.g. <protein> <ref> <position> <mut>
 
 # output
 
    <protein> <mutation-position> <charge-difference> <molecular-weight-difference> <substitutaion-score> <conservation-score> <prediction>
    
 # contact:
 
 if you run into any problem do raise an "issue", also see wiki for proper use and parameters or contact at aarslan@stanford.edu
 
 # cite
 
 if you use this tool, please cite: 
 
  
