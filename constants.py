"""
Constants and mappings for the result files + further processing
"""

# "N+0.984"  Asn deamidation: 114.042927 +  0.984016
# "Q+0.984"  # Gln deamidation: 128.058578 +  0.984016

# N-Terminal
# "+43.006" # Carbamylation
# "-17.027" # NH3 loss
# "+43.006-17.027" # Carbamylation and NH3 loss

modification_dict = {'+0.984':'',           
                     '+43.006':'',
                     '-17.027':'',
                     '+43.006-17.027':''}

novor_modification_dict = {'(1)':'+15.995',
                            '(2)':'+57.021',
                            '(N-term|0)':'+42.011'}

pepnovoplus_result_mod_dict = {'+16':'+15.995',
                            '^+42':'+42.011',
                            '+57':'+57.021'}

unimod_dict = {'[UNIMOD:1]':'+42.011',
                    '[UNIMOD:4]':'+57.021',
                    '[UNIMOD:35]':'+15.995',
                    '[UNIMOD:765]':'',           # Met-loss case is considered by not starting with M
                    '[UNIMOD:766]':'42.011'}

aa_dict = {
  "G": 57.021464,
  "A": 71.037114,
  "S": 87.032028,
  "P": 97.052764,
  "V": 99.068414,
  "T": 101.047670,
  "C+57.021": 160.030649, # 103.009185 + 57.021464
  "L": 113.084064,
  "I": 113.084064,
  "N": 114.042927,
  "D": 115.026943,
  "Q": 128.058578,
  "K": 128.094963,
  "E": 129.042593,
  "M": 131.040485,
  "H": 137.058912,
  "F": 147.068414,
  "R": 156.101111,
  "Y": 163.063329,
  "W": 186.079313,
  # Amino acid modifications.
  "M+15.995": 147.035400,    # Met oxidation:   131.040485 + 15.994915
#  "N+0.984": 115.026943     # Asn deamidation: 114.042927 +  0.984016
#  "Q+0.984": 129.042594     # Gln deamidation: 128.058578 +  0.984016
  # N-terminal modifications.
  "+42.011": 42.010565,      # Acetylation
#  "+43.006": 43.005814      # Carbamylation
#  "-17.027": -17.026549     # NH3 loss
#  "+43.006-17.027": 25.980265      # Carbamylation and NH3 loss
}

tool_name_plot_name_dict = {'msgfplus_percolator':'MSGF+ with Percolator', 
                            'pepnovoplus':'PepNovo+',  
                            'novor':'Novor', 
                            'casanovo':'CasaNovo' , 
                            'pi_helixnovo':'Pi-HelixNovo', 
                            'contranovo':'ContraNovo', 
                            'instanovo':'InstaNovo'}