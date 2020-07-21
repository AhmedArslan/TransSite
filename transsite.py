#!/usr/bin/env pypy3 python3

#1 predict net charge

import glob
import sys
# seq = sys.argv[1]

protein = sys.argv[1]
aa1 = sys.argv[2]
position = sys.argv[3]
aa2 = sys.argv[4]

def netcharge(seq):
   charge = -0.002
   AACharge={"C":-.045,"D":-.999,"E":-.998,"H":.091\
            ,"K":1,"R":1,"Y":-.001}
   for aa in seq:
      if aa in AACharge:
         charge=charge+AACharge[aa]
      else:
         pass
   return charge

#2 predict_isoelectric_point

scales = {
"EMBOSS":     {'Cterm': 3.6, 'pKAsp': 3.9,  'pKGlu': 4.1, 'pKCys': 8.5, 'pKTyr': 10.1, 'pk_his': 6.5, 'Nterm': 8.6, 'pKLys': 10.8, 'pKArg': 12.5},
"DTASelect":  {'Cterm': 3.1, 'pKAsp': 4.4,  'pKGlu': 4.4, 'pKCys': 8.5, 'pKTyr': 10.0, 'pk_his': 6.5, 'Nterm': 8.0, 'pKLys': 10.0, 'pKArg': 12.0},
"Solomon":    {'Cterm': 2.4, 'pKAsp': 3.9,  'pKGlu': 4.3, 'pKCys': 8.3, 'pKTyr': 10.1, 'pk_his': 6.0, 'Nterm': 9.6, 'pKLys': 10.5, 'pKArg': 12.5}, 
"Sillero":    {'Cterm': 3.2, 'pKAsp': 4.0,  'pKGlu': 4.5, 'pKCys': 9.0, 'pKTyr': 10.0, 'pk_his': 6.4, 'Nterm': 8.2, 'pKLys': 10.4, 'pKArg': 12.0},
"Rodwell":    {'Cterm': 3.1, 'pKAsp': 3.68, 'pKGlu': 4.25,'pKCys': 8.33,'pKTyr': 10.07,'pk_his': 6.0, 'Nterm': 8.0, 'pKLys': 11.5, 'pKArg': 11.5},
"Patrickios": {'Cterm': 4.2, 'pKAsp': 4.2,  'pKGlu': 4.2, 'pKCys': 0.0, 'pKTyr':  0.0, 'pk_his': 0.0, 'Nterm': 11.2,'pKLys': 11.2, 'pKArg': 11.2},
"Wikipedia":  {'Cterm': 3.65,'pKAsp': 3.9,  'pKGlu': 4.07,'pKCys': 8.18,'pKTyr': 10.46,'pk_his': 6.04,'Nterm': 8.2, 'pKLys': 10.54,'pKArg': 12.48},
"Grimsley":   {'Cterm': 3.3, 'pKAsp': 3.5,  'pKGlu': 4.2, 'pKCys': 6.8, 'pKTyr': 10.3, 'pk_his': 6.6, 'Nterm': 7.7, 'pKLys': 10.5, 'pKArg': 12.04},
'Lehninger':  {'Cterm': 2.34,'pKAsp': 3.86, 'pKGlu': 4.25,'pKCys': 8.33,'pKTyr': 10.0, 'pk_his': 6.0, 'Nterm': 9.69,'pKLys': 10.5, 'pKArg': 12.4},
'Bjellqvist': {'Cterm': 3.55,'pKAsp': 4.05, 'pKGlu': 4.45,'pKCys': 9.0, 'pKTyr': 10.0, 'pk_his': 5.98,'Nterm': 7.5, 'pKLys': 10.0, 'pKArg': 12.0},   
'IPC_peptide':{'Cterm': 2.383, 'pKAsp': 3.887, 'pKGlu': 4.317, 'pKCys': 8.297, 'pKTyr': 10.071, 'pk_his': 6.018, 'Nterm': 9.564, 'pKLys': 10.517, 'pKArg': 12.503},    # IPC peptide
'IPC_protein':{'Cterm': 2.869, 'pKAsp': 3.872, 'pKGlu': 4.412, 'pKCys': 7.555, 'pKTyr': 10.85,  'pk_his': 5.637, 'Nterm': 9.094, 'pKLys': 9.052,  'pKArg': 11.84},     # IPC protein 
'Toseland':   {'Cterm': 3.19,'pKAsp': 3.6,  'pKGlu': 4.29,'pKCys': 6.87,'pKTyr': 9.61, 'pk_his': 6.33,'Nterm': 8.71, 'pKLys': 10.45, 'pKArg':  12},
'Thurlkill':  {'Cterm': 3.67,'pKAsp': 3.67, 'pKGlu': 4.25,'pKCys': 8.55,'pKTyr': 9.84, 'pk_his': 6.54,'Nterm': 8.0, 'pKLys': 10.4, 'pKArg': 12.0},
'Nozaki':     {'Cterm': 3.8, 'pKAsp': 4.0,  'pKGlu': 4.4, 'pKCys': 9.5, 'pKTyr': 9.6,  'pk_his': 6.3, 'Nterm': 7.5, 'pKLys': 10.4, 'pKArg': 12},   
'Dawson':     {'Cterm': 3.2, 'pKAsp': 3.9,  'pKGlu': 4.3, 'pKCys': 8.3, 'pKTyr': 10.1, 'pk_his': 6.0, 'Nterm': 8.2, 'pKLys': 10.5, 'pKArg':  12},   
          }

aaDict = {'Asp':'D', 'Glu':'E', 'Cys':'C', 'Tyr':'Y', 'His':'H', 
          'Lys':'K', 'Arg':'R', 'Met':'M', 'Phe':'F', 'Leu':'L', 
          'Val':'V', 'Ala':'A', 'Gly':'G', 'Gln':'Q', 'Asn':'N',
          'Ile':'I', 'Trp':'W', 'Ser':'S', 'Thr':'T', 'Sec':'U',
          'Pro':'P', 'Xaa':'X', 'Sec':'U', 'Pyl':'O', 'Asx':'B',
          'Xle':'J', }

acidic = ['D', 'E', 'C', 'Y']
basic = ['K', 'R', 'H']

pKcterminal = {'D': 4.55, 'E': 4.75} 
pKnterminal = {'A': 7.59, 'M': 7.0, 'S': 6.93, 'P': 8.36, 'T': 6.82, 'V': 7.44, 'E': 7.7} 

#N-teminus, middle, C-terminus
promost ={
'K':[10.00,  9.80, 10.30],
'R':[11.50, 12.50, 11.50],
'H':[ 4.89,  6.08,  6.89],
'D':[ 3.57,  4.07,  4.57],
'E':[ 4.15,  4.45,  4.75],
'C':[ 8.00,  8.28,  9.00],
'Y':[ 9.34,  9.84, 10.34],
'U':[ 5.20,  5.43,  5.60], # ref (http://onlinelibrary.wiley.com/doi/10.1002/bip.21581/pdf)
}
           
promost_mid = {
"G":[7.50, 3.70],
"A":[7.58, 3.75],
"S":[6.86, 3.61],
"P":[8.36, 3.40],
"V":[7.44, 3.69],
"T":[7.02, 3.57],
"C":[8.12, 3.10],
"I":[7.48, 3.72],
"L":[7.46, 3.73],
"J":[7.46, 3.73],
"N":[7.22, 3.64],
"D":[7.70, 3.50],
"Q":[6.73, 3.57],
"K":[6.67, 3.40],
"E":[7.19, 3.50],
"M":[6.98, 3.68],
"H":[7.18, 3.17],
"F":[6.96, 3.98],
"R":[6.76, 3.41],
"Y":[6.83, 3.60],
"W":[7.11, 3.78],
"X":[7.26, 3.57], 
"Z":[6.96, 3.535], 
'B':[7.46, 3.57],  
'U':[5.20, 5.60], 
'O':[7.00, 3.50],     
}

def predict_isoelectric_point_ProMoST(seq):
    '''Calculate isoelectric point using ProMoST model'''
    NQ = 0.0
    pH = 6.51             
    pHprev = 0.0         
    pHnext = 14.0        
    E = 0.01            
    temp = 0.01
    while 1:
            if seq[0] in promost.keys():   QN1=-1.0/(1.0+pow(10,(promost[seq[0]][2]-pH)))
            else: QN1=-1.0/(1.0+pow(10,(promost_mid[seq[0]][1]-pH)))
            #print 
            if seq[-1] in promost.keys():  QP2= 1.0/(1.0+pow(10,(pH-promost[seq[-1]][0])))
            else: QP2=1.0/(1.0+pow(10,(pH-promost_mid[seq[-1]][0])))
            
            QN2=-seq.count('D')/(1.0+pow(10,(promost['D'][1]-pH)))           
            QN3=-seq.count('E')/(1.0+pow(10,(promost['E'][1]-pH)))           
            QN4=-seq.count('C')/(1.0+pow(10,(promost['C'][1]-pH)))           
            QN5=-seq.count('Y')/(1.0+pow(10,(promost['Y'][1]-pH)))        
            QP1= seq.count('H')/(1.0+pow(10,(pH-promost['H'][1])))                            
            QP3= seq.count('K')/(1.0+pow(10,(pH-promost['K'][1])))           
            QP4= seq.count('R')/(1.0+pow(10,(pH-promost['R'][1])))                
        
            NQ=QN1+QN2+QN3+QN4+QN5+QP1+QP2+QP3+QP4  

            if NQ<0.0:                                
                    temp = pH
                    pH = pH-((pH-pHprev)/2.0)
                    pHnext = temp
                  
            else:
                    temp = pH
                    pH = pH + ((pHnext-pH)/2.0)
                    pHprev = temp
                   

            if (pH-pHprev<E) and (pHnext-pH<E): 
                    return pH         

#3 molecular_weight
def calculate_molecular_weight(seq):
    """molecular weight"""
    massDict = {'D':115.0886, 'E': 129.1155, 'C': 103.1388, 'Y':163.1760, 'H':137.1411, 
           'K':128.1741, 'R': 156.1875,  'M': 131.1926, 'F':147.1766, 'L':113.1594,
           'V':99.1326,  'A':  71.0788,  'G':  57.0519, 'Q':128.1307, 'N':114.1038, 
           'I':113.1594, 'W': 186.2132,  'S':  87.0782, 'P': 97.1167, 'T':101.1051, 
           'U':141.05,   'h2o':18.01524, 'X':        0, 'Z':128.6231, 'O':255.31, 
           'B':114.5962, 'J': 113.1594,}
    molecular_weight = massDict['h2o']
    for aa in seq:
            molecular_weight+=massDict[aa]       
    return molecular_weight

#4 biochemical properties

def blo(aa1, aa2):
  with open("/path/to/Blosum62.txt") as bp:
    for p in bp:
      p = p.rstrip().split('\t')
      if p[0] == aa1 and p[-1] == aa2:
        return p[0], "to", p[-1], p[1]

def mts(protein, aa1, position, aa2):
	
  files = glob.glob("/path/to/conservation_mouse3/*")
  for f in files:
    f1 = f.split('/')
    f2 = f1[-1].split(".")
    if protein == f2[0]:
      with open(f) as ij:
        i = ij.readline().split()
        print('Net charge:',round(netcharge(i[0])))
        print("Isoelectric point:", round(predict_isoelectric_point_ProMoST(i[0]), 2))
        print('Molecular weight:', round(calculate_molecular_weight(i[0])),'Da')
        print(i)

        for j in ij.readlines()[0:]:
          j = j.rstrip().split('\t')
          if position == j[0] and aa1 == j[1]:
            s = i[0][int(j[0])-2]+aa2+i[0][int(j[0]):-1]
            print(s)
            print('Variant form Net charge:',round(netcharge(s)))
            print("Variant form Isoelectric point:", round(predict_isoelectric_point_ProMoST(i[0]), 2))
            print('Variant form Molecular weight:', round(calculate_molecular_weight(s)),'Da')
            print("Biochemical properties change from: ", blo(aa1, aa2))
            if float(j[-1]) > 1:
              print('variable site: ' +j[-1]+ " score")
            if float(j[-1]) < 1 and float(j[-1]) > 0.9:
              print('moderataely conserveted site: ' +j[-1]+ " score")
            if float(j[-1]) < -0.1:
              print('higly conserveted site: ' +j[-1]+ " score")


if __name__=='__main__':

  import argparse
  sys.setrecursionlimit(2000)

  parser = argparse.ArgumentParser()

  parser.add_argument('-p','--protein', required=True, help='Protein name is required')
  parser.add_argument('-ref', '--reference',  required=True, help='reference amino acid is missing')
  parser.add_argument('-pos', '--position',  required=True, help='reference amino acid position is missing')
  parser.add_argument('-alt', '--alteration',  required=True, help='altered amino acid is missing')
  args = vars(parser.parse_args())

  if args:
    try:
      mts(protein, aa1, position, aa2)

    except IOError:
      pass
  else:
    print ("to run MP-Mup seek help")
    #Molecular Transportation via mutational Upshot (MP-MUp)?

    #python3 mpmup.py -p [name] [ref] [position] [alt]
