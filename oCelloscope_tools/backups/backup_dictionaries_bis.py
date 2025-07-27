import numpy as np

antimycotics_dict_A = {'A': 'amphotericin B', 'B': 'voriconazole', 'C': 'posaconazole', 'D': 'itraconazole', 'E': 'isavuconazole', 'F': 'olorofim'}
antimycotics_dict_C = {'A': 'amphotericin B', 'B': '5-flucytosine', 'C': 'voriconazole', 'D': 'posaconazole', 'E': 'anidulafungin', 'F': 'itraconazole', 'G': 'fluconazole', 'H': 'micafungin'}

# EUCAST dictionaries: threshold for growth inhibition compared to growth control
#eucast_A = {'A': 0.05, 'B': 0.5, 'C': 0.5, 'D': 0.5, 'E': 0.5, 'F': 0.5}
#eucast_C = {'A': 0.05, 'B': 0.5, 'C': 0.5, 'D': 0.5, 'E': 0.5, 'F': 0.5, 'G': 0.5, 'H': 0.5}

eucast_A = {'A': 0.05, 'B': 0.05, 'C': 0.05, 'D': 0.05, 'E': 0.05, 'F': 0.05}
eucast_C = {'A': 0.90, 'B': 0.5, 'C': 0.5, 'D': 0.5, 'E': 0.5, 'F': 0.5, 'G': 0.5, 'H': 0.5}

breakpoint_dict_A = {'A_ASFU': [5, 4], 'B_ASFU': [5, 3], 'C_ASFU': [7, 5], 'D_ASFU': [5, 3], 'E_ASFU': [5, 3], 'F_ASFU': [1, 2]}

"""
MIC_well_number >= breakpoint_list[0]:             ref_cat = 'S' 
MIC_well_number <= breakpoint_list[1]:             ref_cat = 'R' 
""" 
# first item of list is the well number (on EUCAST layout) that equals <= S 
# second item of list is the well number (EUCAST layout) that equals >= R 
# well numbers that fall in between these (for example for voriconazole and A. fumigatus), correpsond to ATU or intermediate category 
# note: vorico MIC for A. fumigatus >1 mg/mL is considered resistant, but MIC 2 mg/mL is still considered ATU: both categories are overlapping 
# as it is encoded here, 2 mg/mL is in a separate category from R, unlike EUCAST - however ATU_VME (predicted S instead of ATU) can be penalized as VME while defining thresholds, therefore considering it as R

breakpoint_dict_C = {
'A_CDAL': [1, 2], 'B_CDAL': [np.NaN, np.NaN], 'C_CDAL': [0.06, 0.5], 'D_CDAL': [0.06, 0.12], 'E_CDAL': [0.03, 0.06], 'F_CDAL': [0.06, 0.12], 'G_CDAL': [2, 8], 'H_CDAL': [0.016, 0.064],
'A_CDDU': [1, 2], 'B_CDDU': [np.NaN, np.NaN], 'C_CDDU': [0.06, 0.5], 'D_CDDU': [0.06, 0.12], 'E_CDDU': [np.NaN, np.NaN], 'F_CDDU': [0.06, 0.12], 'G_CDDU': [2, 8], 'H_CDDU': [np.NaN, np.NaN],
'A_CDGL': [1, 2], 'B_CDGL': [np.NaN, np.NaN], 'C_CDGL': [np.NaN, np.NaN], 'D_CDGL': [np.NaN, np.NaN], 'E_CDGL': [0.06, 0.12], 'F_CDGL': [np.NaN, np.NaN], 'G_CDGL': [0.001, 32], 'H_CDGL': [0.03, 0.06],
'A_CDKR': [1, 2], 'B_CDKR': [np.NaN, np.NaN], 'C_CDKR': [np.NaN, np.NaN], 'D_CDKR': [np.NaN, np.NaN], 'E_CDKR': [0.06, 0.12], 'F_CDKR': [np.NaN, np.NaN], 'G_CDKR': [np.NaN, np.NaN], 'H_CDKR': [np.NaN, np.NaN],
'A_CDPA': [1, 2], 'B_CDPA': [np.NaN, np.NaN], 'C_CDPA': [0.125, 0.5], 'D_CDPA': [0.06, 0.12], 'E_CDPA': [4, 8], 'F_CDPA': [0.125, 0.25], 'G_CDPA': [2, 8], 'H_CDPA': [2, 4],
'A_CDTR': [1, 2], 'B_CDTR': [np.NaN, np.NaN], 'C_CDTR': [0.125, 0.5], 'D_CDTR': [0.06, 0.12], 'E_CDTR': [0.06, 0.12], 'F_CDTR': [0.125, 0.25], 'G_CDTR': [2, 8], 'H_CDTR': [np.NaN, np.NaN],
'A_CDAU': [1, 2], 'B_CDAU': [np.NaN, np.NaN], 'C_CDAU': [np.NaN, np.NaN], 'D_CDAU': [np.NaN, np.NaN], 'E_CDAU': [2, 4], 'F_CDAU': [np.NaN, np.NaN], 'G_CDAU': [16, 32], 'H_CDAU': [np.NaN, np.NaN],
'A_other': [np.NaN, np.NaN], 'B_other': [np.NaN, np.NaN], 'C_other': [np.NaN, np.NaN], 'D_other': [np.NaN, np.NaN], 'E_other': [np.NaN, np.NaN], 'F_other': [np.NaN, np.NaN], 'G_other': [2, 8], 'H_other': [np.NaN, np.NaN]}

breakpoint_dict_C = {
'A_CDAL': [4, 3], 'B_CDAL': [np.NaN, np.NaN], 'C_CDAL': [8, 5], 'D_CDAL': [8, 7], 'E_CDAL': [9, 8], 'F_CDAL': [8, 7], 'G_CDAL': [7, 5], 'H_CDAL': [10, 8],
'A_CDDU': [4, 3], 'B_CDDU': [np.NaN, np.NaN], 'C_CDDU': [8, 5], 'D_CDDU': [8, 7], 'E_CDDU': [np.NaN, np.NaN], 'F_CDDU': [8, 7], 'G_CDDU': [7, 5], 'H_CDDU': [np.NaN, np.NaN],
'A_CDGL': [4, 3], 'B_CDGL': [np.NaN, np.NaN], 'C_CDGL': [np.NaN, np.NaN], 'D_CDGL': [np.NaN, np.NaN], 'E_CDGL': [8, 7], 'F_CDGL': [np.NaN, np.NaN], 'G_CDGL': [11, 3], 'H_CDGL': [9, 8],
'A_CDKR': [4, 3], 'B_CDKR': [np.NaN, np.NaN], 'C_CDKR': [np.NaN, np.NaN], 'D_CDKR': [np.NaN, np.NaN], 'E_CDKR': [8, 7], 'F_CDKR': [np.NaN, np.NaN], 'G_CDKR': [np.NaN, np.NaN], 'H_CDKR': [np.NaN, np.NaN],
'A_CDPA': [4, 3], 'B_CDPA': [np.NaN, np.NaN], 'C_CDPA': [7, 5], 'D_CDPA': [8, 7], 'E_CDPA': [2, 1], 'F_CDPA': [7, 6], 'G_CDPA': [7, 5], 'H_CDPA': [3, 2],
'A_CDTR': [4, 3], 'B_CDTR': [np.NaN, np.NaN], 'C_CDTR': [7, 5], 'D_CDTR': [8, 7], 'E_CDTR': [8, 7], 'F_CDTR': [7, 6], 'G_CDTR': [7, 5], 'H_CDTR': [np.NaN, np.NaN],
'A_CDAU': [4, 3], 'B_CDAU': [np.NaN, np.NaN], 'C_CDAU': [np.NaN, np.NaN], 'D_CDAU': [np.NaN, np.NaN], 'E_CDAU': [3, 2], 'F_CDAU': [np.NaN, np.NaN], 'G_CDAU': [4, 3], 'H_CDAU': [np.NaN, np.NaN],
'A_other': [np.NaN, np.NaN], 'B_other': [np.NaN, np.NaN], 'C_other': [np.NaN, np.NaN], 'D_other': [np.NaN, np.NaN], 'E_other': [np.NaN, np.NaN], 'F_other': [np.NaN, np.NaN], 'G_other': [7, 5], 'H_other': [np.NaN, np.NaN]}

MIC_dict_A = {'A-94329739 60 scan areas.xlsx':
                {'A': 'A5', 'B': 'B4', 'C': 'C5', 'D': 'D0', 'E': 'E2', 'F': 'Z0'}, # twijfel vorico: B1, 2, 3? - geen foto
            'A-93624221 72 scan areas.xlsx':
                {'A': 'A8', 'B': 'B7', 'C': 'C7', 'D': 'D8', 'E': 'E6', 'F': 'F8'}, # olorofim F8 MIC obv foto; eerder afgelezen F9
#             'A-94277936 60 scan areas.xlsx':
#                 {'A': 'A3', 'B': 'B0', 'C': 'C5', 'D': 'D6', 'E': 'E0', 'F': 'Z0'}, # slechts 11u data beschikbaar; geen foto
            'A-94512364 60 scan areas.xlsx':
                {'A': 'A7', 'B': 'B4', 'C': 'C5', 'D': 'D0', 'E': 'E2', 'F': 'F7'}, # slechts 33u data beschikbaar; isa eerst 4 mg/L, maar lijkt eerder 8 mg/L obv foto
            'A-94647258 60 scan areas.xlsx':
                {'A': 'A6', 'B': 'B0', 'C': 'C3', 'D': 'D0', 'E': 'E0', 'F': 'Z0'},
            'A-ATCC-204305 72 scan areas.xlsx':
                {'A': 'A7', 'B': 'B7', 'C': 'C9', 'D': 'D9', 'E': 'E6', 'F': 'F9'},
            'A-CYP-15-232 72 scan areas.xlsx':
                {'A': 'A6', 'B': 'B2', 'C': 'C4', 'D': 'D0', 'E': 'E1', 'F': 'F7'}, # obv foto isavu van E2 naar E1 en MIC oloro van F8 naar F7
            'A-CYP-15-183 72 scan areas.xlsx':
                {'A': 'A7', 'B': 'B2', 'C': 'C4', 'D': 'D0', 'E': 'E1', 'F': 'F8'}, # McF 0.38
            'A-CYP-15-234 72 scan areas.xlsx':
                {'A': 'A7', 'B': 'B3', 'C': 'C5', 'D': 'D0', 'E': 'E2', 'F': 'F7'}, # obv foto gewijzigd van F8 naar F7 (doch wat twijfel of zelfs niet F5 moet zijn)
            'A-CYP-15-229 72 scan areas.xlsx':
                {'A': 'A6', 'B': 'B3', 'C': 'C4', 'D': 'D0', 'E': 'E2', 'F': 'F6'}, # geen filter, cultuur 24u, afgelezen 69u - evt B2 ipv B3; evt E1 ipv E2, van F8 naar F6 obv foto
            'A-CYP-15-202 72 scan areas.xlsx':
                {'A': 'A7', 'B': 'B3', 'C': 'C5', 'D': 'D0', 'E': 'E2', 'F': 'F6'}, # geen filter; obv foto A8 naar A7; B4 naar B3; E3 naar E2; F8 naar F6 (potentieel F5 of F4)
            'A-CYP-15-231 96 scan areas.xlsx':
                {'A': 'A8', 'B': 'B4', 'C': 'C5', 'D': 'D0', 'E': 'E2', 'F': 'F8'}, # obv foto E2 ipv E3; F8 ipv F9 (moet misschien zelfs F7 zijn)
            'A-CYP-15-235 96 scan areas.xlsx':
                {'A': 'A7', 'B': 'B3', 'C': 'C4', 'D': 'D0', 'E': 'E2', 'F': 'F8'}, # cultuur 8d oud; obv foto E3 naar E2, F9 naar F8
            'A-CYP-15-207 96 scan areas.xlsx':
                {'A': 'A8', 'B': 'B4', 'C': 'C5', 'D': 'D0', 'E': 'E3', 'F': 'F8'}, # McF 0.57; obv foto F9 naar F8
            'A-CYP-15-208 96 scan areas.xlsx':
                {'A': 'A7', 'B': 'B3', 'C': 'C4', 'D': 'D0', 'E': 'E3', 'F': 'F9'}, # filter doch veel hyfen; obv foto B4 naar B3, C5 naar C4
            'A-CYP-15-198 72 scan areas.xlsx':
                {'A': 'A7', 'B': 'B1', 'C': 'C4', 'D': 'D0', 'E': 'E1', 'F': 'F6'}, # B2 naar B1, C5 naar C4, F7 naar F6
            'A-94991776 60 scan areas.xlsx':
                {'A': 'A7', 'B': 'B5', 'C': 'C6', 'D': 'D6', 'E': 'E5', 'F': 'Z0'}, # Afgelezen na 71u; plaat langdurig op 42 °C
#             'A-94742083 96 scan areas.xlsx':
#                 {'A': 'A6', 'B': 'B2', 'C': 'C4', 'D': 'D0', 'E': 'E2', 'F': 'Z0'}, # Slechts 2u gelopen! Ruwe data nog zoeken in map; nu MICs overgenomen van LWS en vanop foto; klinische TR34 positieve stam!
            'A-94943434 60 scan areas.xlsx':
                {'A': 'A6', 'B': 'B6', 'C': 'C7', 'D': 'D6', 'E': 'E6', 'F': 'Z0'}, # McF 0.67; afgelezen 53u
           ### DEEL 2
            'A-CYP-15-193 60 scan areas.xlsx':
                {'A': 'A6', 'B': 'B0', 'C': 'C6', 'D': 'D0', 'E': 'E0', 'F': 'Z0'}, # McF 0.46; afgelezen 49u; Eagle type res voor itra           
            'A-96584145 60 scan areas.xlsx':
                {'A': 'A6', 'B': 'B2', 'C': 'C6', 'D': 'D0', 'E': 'E2', 'F': 'Z0'}, # McF 0.60; afgelezen 49u
            'A-96264824 60 scan areas.xlsx':
                {'A': 'A6', 'B': 'B0', 'C': 'C4', 'D': 'D6', 'E': 'E0', 'F': 'Z0'}, # McF?; afgelezen 47u30           
            'A-96043965 72 scan areas.xlsx':
                {'A': 'A7', 'B': 'B6', 'C': 'C7', 'D': 'D7', 'E': 'E5', 'F': 'F7'}, # McF 0.50; afgelezen 48u            
            ############### VAN A-95746873 is er enkel .txt file, geen Excel ################
            'A-95746873 60 scan areas.xlsx':
                {'A': 'A6', 'B': 'B6', 'C': 'C8', 'D': 'D7', 'E': 'E5', 'F': 'Z0'}, # McF 0.70; afgelezen 58u45!!!                     
            'A-CYP-15-209 72 scan areas.xlsx':
                {'A': 'A6', 'B': 'B0', 'C': 'C4', 'D': 'D4', 'E': 'E0', 'F': 'F8'}, # McF 0.25; afgelezen 50u?           
            'A-95707146 72 scan areas.xlsx':
                {'A': 'A7', 'B': 'B5', 'C': 'C7', 'D': 'D8', 'E': 'E5', 'F': 'F7'}, # McF?; afgelezen 47u30 - geen foto als ref beschikbaar ondanks dubieuze MICs voor vor en posaco         
            'A-95551158 72 scan areas.xlsx':
                {'A': 'A7', 'B': 'B6', 'C': 'C8', 'D': 'D8', 'E': 'E5', 'F': 'F6'}, # McF 0.42; afgelezen 44u30       
            'A-95761404 72 scan areas.xlsx':
                {'A': 'A7', 'B': 'B4', 'C': 'C6', 'D': 'D6', 'E': 'E4', 'F': 'F7'}, # McF 0.52; afgelezen 58u               
            'A-CYP-15-224 72 scan areas.xlsx':
                {'A': 'A6', 'B': 'B0', 'C': 'C4', 'D': 'D3', 'E': 'E0', 'F': 'F7'}, # McF 0.35; afgelezen 50u - PLAAT NIET TERUGGEVONDEN - TWIJFEL OLOROFIM? F1 vs F7?             
            'A-CYP-15-245 72 scan areas.xlsx':
                {'A': 'A6', 'B': 'B3', 'C': 'C4', 'D': 'D0', 'E': 'E2', 'F': 'F6'}, # McF 0.42; afgelezen 33u          
            'A-CYP-15-250 72 scan areas.xlsx':
                {'A': 'A6', 'B': 'B0', 'C': 'C4', 'D': 'D0', 'E': 'E0', 'F': 'Z0'}, # McF 0.43; afgelezen 47u 
            'A-CYP-15-161 72 scan areas.xlsx':
                {'A': 'A6', 'B': 'B0', 'C': 'C5', 'D': 'D6', 'E': 'E0', 'F': 'F7'}, # McF 0.56; afgelezen 56u 
            'A-CYP-15-243 72 scan areas.xlsx':
                {'A': 'A6', 'B': 'B0', 'C': 'C5', 'D': 'D6', 'E': 'E0', 'F': 'F7'}, # McF 0.65; afgelezen na 4d waarvan 2d koelkast
            'A-94939738 60 scan areas.xlsx':
                {'A': 'A6', 'B': 'B2', 'C': 'C5', 'D': 'D0', 'E': 'E2', 'F': 'Z0'}, # Toegevoegd op 28 dec, geen andere info of foto teruggevonden
           }

MIC_dict_C = {'C-Dhondt 96 scan areas.xlsx':
                {'A': 'A5', 'B': 'B9', 'C': 'C0', 'D': 'D9', 'E': 'E5', 'F': 'F6', 'G': 'G3', 'H': 'H7'}, # 
            'C-Alkandri2 96 scan areas.xlsx':
                {'A': 'A4', 'B': 'B9', 'C': 'C3', 'D': 'D6', 'E': 'E3', 'F': 'F5', 'G': 'G0', 'H': 'H6'}, # 
            'C-Enam 96 scan areas.xlsx':
                {'A': 'A4', 'B': 'B9', 'C': 'C5', 'D': 'D7', 'E': 'E1', 'F': 'F6', 'G': 'G2', 'H': 'H0'}, # MICs staal 21/11 gebruikt, gezien afgelezen 24u
            'C-Vansteelandt2 96 scan areas.xlsx':
                {'A': 'A5', 'B': 'B8', 'C': 'C3', 'D': 'D6', 'E': 'E5', 'F': 'F6', 'G': 'G1', 'H': 'H7'}, # MICs staal 19/11 gebruikt, idem
            'C-Klein 72 scan areas.xlsx':
                {'A': 'A5', 'B': 'B9', 'C': 'C2', 'D': 'D5', 'E': 'E5', 'F': 'F4', 'G': 'G1', 'H': 'H8'}, # Nog correleren met foto
#             'C-94768530-CDGL 96 scan areas.xlsx':
#                 {'A': 'A5', 'B': 'B13', 'C': 'C6', 'D': 'D4', 'E': 'E8', 'F': 'F5', 'G': 'G5', 'H': 'H10'}, # Nog correleren met foto
            'C-94768530-CDGL 96 scan areas.xlsx':
                {'A': 'A5', 'B': 'B11', 'C': 'C6', 'D': 'D4', 'E': 'E8', 'F': 'F5', 'G': 'G5', 'H': 'H10'}, # Nog correleren met foto
            'C-94878906-CDKR 96 scan areas.xlsx':
                {'A': 'A4', 'B': 'B5', 'C': 'C5', 'D': 'D7', 'E': 'E9', 'F': 'F7', 'G': 'G3', 'H': 'H7'}, 
            'C-ATCC-CDPA 96 scan areas.xlsx':
                {'A': 'A5', 'B': 'B9', 'C': 'C9', 'D': 'D9', 'E': 'E4', 'F': 'F8', 'G': 'G7', 'H': 'H4'}, 
#             'C-ATCC-CDAL 96 scan areas':
#                 {'A': 'A5', 'B': 'B6', 'C': 'C13', 'D': 'D10', 'E': 'E13', 'F': 'F10', 'G': 'G10', 'H': 'H13'}, # MICs afgelezen op plaat 54u geïncubeerd - geen andere refmics beschikbaar?
            'C-ATCC-CDAL 96 scan areas.xlsx': 
                {'A': 'A5', 'B': 'B6', 'C': 'C11', 'D': 'D10', 'E': 'E11', 'F': 'F10', 'G': 'G10', 'H': 'H11'}, # MICs niet overeenkomend met EXCEL? Eerdere nota: MICs afgelezen op plaat 54u geïncubeerd - geen andere refmics beschikbaar?
            'C-ATCC-CDTR3 96 scan areas.xlsx':
                {'A': 'A5', 'B': 'B10', 'C': 'C7', 'D': 'D7', 'E': 'E7', 'F': 'F7', 'G': 'G9', 'H': 'H10'}, # rMICs afgelezen op plaat 34 u geïncubeerd
            'C-MCC3-10 96 scan areas.xlsx':
                {'A': 'A5', 'B': 'B9', 'C': 'C9', 'D': 'D9', 'E': 'E6', 'F': 'F7', 'G': 'G7', 'H': 'H8'}, # 0.66 mcF; afgelezen 26u20
            'C-MRU224 96 scan areas.xlsx':
                {'A': 'A5', 'B': 'B0', 'C': 'C9', 'D': 'D0', 'E': 'E3', 'F': 'F0', 'G': 'G0', 'H': 'H2'}, # MOEILIJKE AFLEZING AZOLES: schommelen rond 50% inhibitie - 0.5 mcF; afgelezen 23u30
            'C-B11220-AR381 96 scan areas.xlsx':
                {'A': 'A6', 'B': 'B9', 'C': 'C9', 'D': 'D10', 'E': 'E9', 'F': 'F10', 'G': 'G6', 'H': 'H9'}, # afgelezen 26u30; 0.55 mcF
            'C-B8441-AR387 96 scan areas.xlsx':
                {'A': 'A6', 'B': 'B10', 'C': 'C9', 'D': 'D9', 'E': 'E3', 'F': 'F8', 'G': 'G7', 'H': 'H5'}, # 0.40 mcF, afgelezen 26u30
            'C-ML315-AR1037 96 scan areas.xlsx':
                {'A': 'A5', 'B': 'B9', 'C': 'C4', 'D': 'D6', 'E': 'E4', 'F': 'F6', 'G': 'G3', 'H': 'H0'}, # 0.62 mcF (of 0.52?); afgelezen 23u30
            'C-96299088 96 scan areas.xlsx':
                {'A': 'A5', 'B': 'B11', 'C': 'C7', 'D': 'D8', 'E': 'E9', 'F': 'F8', 'G': 'G7', 'H': 'H9'}, # 0.56 mcF; afgelezen 22u30
            'C-95738583 96 scan areas.xlsx':
                {'A': 'A6', 'B': 'B10', 'C': 'C6', 'D': 'D5', 'E': 'E9', 'F': 'F6', 'G': 'G5', 'H': 'H13'}, # 0.51 mcF; afgelezen 26u        
            'C-95690277 96 scan areas.xlsx':
                {'A': 'A6', 'B': 'B13', 'C': 'C10', 'D': 'D10', 'E': 'E13', 'F': 'F13', 'G': 'G8', 'H': 'H13'}, # NOG FOTO bekijken; amfo-B MIC op curve lijkt A8, maar is net niet onder cutoff 95% inhibitie, vandaag A6; afgelezen 25u
            'C-CDKR-IQC 96 scan areas.xlsx':
                {'A': 'A5', 'B': 'B4', 'C': 'C7', 'D': 'D7', 'E': 'E9', 'F': 'F7', 'G': 'G4', 'H': 'H8'}, # NOG FOTO bekijken
            'C-ATCC-CDGL 96 scan areas.xlsx':
                {'A': 'A5', 'B': 'B13', 'C': 'C6', 'D': 'D4', 'E': 'E9', 'F': 'F4', 'G': 'G4', 'H': 'H10'}} # NOG FOTO bekijken}

#all_parameters = ['BCA', 'BCA_rate', 'BCA_mean_rate', 'BCA_cumulative', 'BCA_cumulative_rate', 'BCA_cumulative_mean_rate', 'BCA_cumulative_norm', 'BCA_ratio', 'BCA_rate_ratio', 'BCA_mean_rate_ratio', 'BCA_cumulative_ratio', 'BCA_cumulative_rate_ratio', 'BCA_cumulative_mean_rate_ratio', 'BCA_cumulative_norm_ratio', 'BCANormalized', 'BCANormalized_rate', 'BCANormalized_mean_rate', 'BCANormalized_cumulative', 'BCANormalized_cumulative_rate', 'BCANormalized_cumulative_mean_rate', 'BCANormalized_cumulative_norm', 'BCANormalized_ratio', 'BCANormalized_rate_ratio', 'BCANormalized_mean_rate_ratio', 'BCANormalized_cumulative_ratio', 'BCANormalized_cumulative_rate_ratio', 'BCANormalized_cumulative_mean_rate_ratio', 'BCANormalized_cumulative_norm_ratio', 'SESAfungi', 'SESAfungi_rate', 'SESAfungi_mean_rate', 'SESAfungi_cumulative', 'SESAfungi_cumulative_rate', 'SESAfungi_cumulative_mean_rate', 'SESAfungi_cumulative_norm', 'SESAfungi_ratio', 'SESAfungi_rate_ratio', 'SESAfungi_mean_rate_ratio', 'SESAfungi_cumulative_ratio', 'SESAfungi_cumulative_rate_ratio', 'SESAfungi_cumulative_mean_rate_ratio', 'SESAfungi_cumulative_norm_ratio', 'SESAfungiNormalized', 'SESAfungiNormalized_rate', 'SESAfungiNormalized_mean_rate', 'SESAfungiNormalized_cumulative', 'SESAfungiNormalized_cumulative_rate', 'SESAfungiNormalized_cumulative_mean_rate', 'SESAfungiNormalized_cumulative_norm', 'SESAfungiNormalized_ratio', 'SESAfungiNormalized_rate_ratio', 'SESAfungiNormalized_mean_rate_ratio', 'SESAfungiNormalized_cumulative_ratio', 'SESAfungiNormalized_cumulative_rate_ratio', 'SESAfungiNormalized_cumulative_mean_rate_ratio', 'SESAfungiNormalized_cumulative_norm_ratio', 'TA', 'TA_rate', 'TA_mean_rate', 'TA_cumulative', 'TA_cumulative_rate', 'TA_cumulative_mean_rate', 'TA_cumulative_norm', 'TA_ratio', 'TA_rate_ratio', 'TA_mean_rate_ratio', 'TA_cumulative_ratio', 'TA_cumulative_rate_ratio', 'TA_cumulative_mean_rate_ratio', 'TA_cumulative_norm_ratio', 'TANormalized', 'TANormalized_rate', 'TANormalized_mean_rate', 'TANormalized_cumulative', 'TANormalized_cumulative_rate', 'TANormalized_cumulative_mean_rate', 'TANormalized_cumulative_norm', 'TANormalized_ratio', 'TANormalized_rate_ratio', 'TANormalized_mean_rate_ratio', 'TANormalized_cumulative_ratio', 'TANormalized_cumulative_rate_ratio', 'TANormalized_cumulative_mean_rate_ratio', 'TANormalized_cumulative_norm_ratio']

all_parameters = ['BCA', 'BCA_rate', 'BCA_mean_rate', 'BCA_cumulative', 'BCA_cumulative_rate', 'BCA_cumulative_mean_rate', 'BCA_cumulative_norm', 'BCA_ratio', 'BCA_rate_old', 'BCA_rate_old_time_corrected', 'BCA_rate_ratio', 'BCA_mean_rate_ratio', 'BCA_cumulative_ratio', 'BCA_cumulative_rate_ratio', 'BCA_cumulative_mean_rate_ratio', 'BCA_cumulative_norm_ratio', 'BCA_rate_old_ratio', 'BCA_rate_old_time_corrected_ratio', 'BCANormalized', 'BCANormalized_rate', 'BCANormalized_mean_rate', 'BCANormalized_cumulative', 'BCANormalized_cumulative_rate', 'BCANormalized_cumulative_mean_rate', 'BCANormalized_cumulative_norm', 'BCANormalized_ratio', 'BCANormalized_rate_old', 'BCANormalized_rate_old_time_corrected', 'BCANormalized_rate_ratio', 'BCANormalized_mean_rate_ratio', 'BCANormalized_cumulative_ratio', 'BCANormalized_cumulative_rate_ratio', 'BCANormalized_cumulative_mean_rate_ratio', 'BCANormalized_cumulative_norm_ratio', 'BCANormalized_rate_old_ratio', 'BCANormalized_rate_old_time_corrected_ratio', 'SESAfungi', 'SESAfungi_rate', 'SESAfungi_mean_rate', 'SESAfungi_cumulative', 'SESAfungi_cumulative_rate', 'SESAfungi_cumulative_mean_rate', 'SESAfungi_cumulative_norm', 'SESAfungi_ratio', 'SESAfungi_rate_old', 'SESAfungi_rate_old_time_corrected', 'SESAfungi_rate_ratio', 'SESAfungi_mean_rate_ratio', 'SESAfungi_cumulative_ratio', 'SESAfungi_cumulative_rate_ratio', 'SESAfungi_cumulative_mean_rate_ratio', 'SESAfungi_cumulative_norm_ratio', 'SESAfungi_rate_old_ratio', 'SESAfungi_rate_old_time_corrected_ratio', 'SESAfungiNormalized', 'SESAfungiNormalized_rate', 'SESAfungiNormalized_mean_rate', 'SESAfungiNormalized_cumulative', 'SESAfungiNormalized_cumulative_rate', 'SESAfungiNormalized_cumulative_mean_rate', 'SESAfungiNormalized_cumulative_norm', 'SESAfungiNormalized_ratio', 'SESAfungiNormalized_rate_old', 'SESAfungiNormalized_rate_old_time_corrected', 'SESAfungiNormalized_rate_ratio', 'SESAfungiNormalized_mean_rate_ratio', 'SESAfungiNormalized_cumulative_ratio', 'SESAfungiNormalized_cumulative_rate_ratio', 'SESAfungiNormalized_cumulative_mean_rate_ratio', 'SESAfungiNormalized_cumulative_norm_ratio', 'SESAfungiNormalized_rate_old_ratio', 'SESAfungiNormalized_rate_old_time_corrected_ratio', 'TA', 'TA_rate', 'TA_mean_rate', 'TA_cumulative', 'TA_cumulative_rate', 'TA_cumulative_mean_rate', 'TA_cumulative_norm', 'TA_ratio', 'TA_rate_old', 'TA_rate_old_time_corrected', 'TA_rate_ratio', 'TA_mean_rate_ratio', 'TA_cumulative_ratio', 'TA_cumulative_rate_ratio', 'TA_cumulative_mean_rate_ratio', 'TA_cumulative_norm_ratio', 'TA_rate_old_ratio', 'TA_rate_old_time_corrected_ratio', 'TANormalized', 'TANormalized_rate', 'TANormalized_mean_rate', 'TANormalized_cumulative', 'TANormalized_cumulative_rate', 'TANormalized_cumulative_mean_rate', 'TANormalized_cumulative_norm', 'TANormalized_ratio', 'TANormalized_rate_old', 'TANormalized_rate_old_time_corrected', 'TANormalized_rate_ratio', 'TANormalized_mean_rate_ratio', 'TANormalized_cumulative_ratio', 'TANormalized_cumulative_rate_ratio', 'TANormalized_cumulative_mean_rate_ratio', 'TANormalized_cumulative_norm_ratio', 'TANormalized_rate_old_ratio', 'TANormalized_rate_old_time_corrected_ratio']

species_dict_A = {'A-94329739 60 scan areas.xlsx': 'ASFU', 
            'A-93624221 72 scan areas.xlsx': 'ASFU', 
            'A-94277936 60 scan areas.xlsx': 'ASFU', 
            'A-94512364 60 scan areas.xlsx': 'ASFU',
            'A-94647258 60 scan areas.xlsx': 'ASFU',
            'A-ATCC-204305 72 scan areas.xlsx': 'ASFU', 
            'A-CYP-15-232 72 scan areas.xlsx': 'ASFU', 
            'A-CYP-15-183 72 scan areas.xlsx': 'ASFU',
            'A-CYP-15-234 72 scan areas.xlsx': 'ASFU', 
            'A-CYP-15-229 72 scan areas.xlsx': 'ASFU',
            'A-CYP-15-202 72 scan areas.xlsx': 'ASFU',
            'A-CYP-15-231 96 scan areas.xlsx': 'ASFU',
            'A-CYP-15-235 96 scan areas.xlsx': 'ASFU',
            'A-CYP-15-207 96 scan areas.xlsx': 'ASFU',
            'A-CYP-15-208 96 scan areas.xlsx': 'ASFU',
            'A-CYP-15-198 72 scan areas.xlsx': 'ASFU',
            'A-94991776 60 scan areas.xlsx': 'ASFU',
            'A-94742083 96 scan areas.xlsx': 'ASFU',
            'A-94943434 60 scan areas.xlsx': 'ASFU',
            'A-CYP-15-193 60 scan areas.xlsx': 'ASFU',
            'A-96584145 60 scan areas.xlsx': 'ASFU',
            'A-96264824 60 scan areas.xlsx': 'ASFU',
            'A-96043965 72 scan areas.xlsx': 'ASFU',
            'A-95746873 60 scan areas.xlsx': 'ASFU',
            'A-CYP-15-209 72 scan areas.xlsx': 'ASFU',
            'A-95707146 72 scan areas.xlsx': 'ASFU',
            'A-95551158 72 scan areas.xlsx': 'ASFU',
            'A-95761404 72 scan areas.xlsx': 'ASFU',
            'A-CYP-15-224 72 scan areas.xlsx': 'ASFU',
            'A-CYP-15-245 72 scan areas.xlsx': 'ASFU', 
            'A-CYP-15-250 72 scan areas.xlsx': 'ASFU',
            'A-CYP-15-161 72 scan areas.xlsx': 'ASFU',
            'A-CYP-15-243 72 scan areas.xlsx': 'ASFU',
            'A-94939738 60 scan areas.xlsx': 'ASFU'
           }

# ['C-95690277 96 scan areas.xlsx', 'C-MRU224 96 scan areas.xlsx']

species_dict_C = {'C-Dhondt 96 scan areas.xlsx': 'CDAU',
            'C-Alkandri2 96 scan areas.xlsx': 'CDAU',
            'C-Enam 96 scan areas.xlsx': 'CDAU',
            'C-Vansteelandt2 96 scan areas.xlsx': 'CDAU',
            'C-Klein 72 scan areas.xlsx': 'CDAU',
            'C-94768530-CDGL 96 scan areas.xlsx': 'CDGL',
            'C-94768530-CDGL 96 scan areas.xlsx': 'CDGL',
            'C-94878906-CDKR 96 scan areas.xlsx': 'CDKR',
            'C-ATCC-CDPA 96 scan areas.xlsx': 'CDPA',
            'C-ATCC-CDAL 96 scan areas': 'CDAL',
            'C-ATCC-CDAL 96 scan areas.xlsx': 'CDAL',
            'C-ATCC-CDTR3 96 scan areas.xlsx': 'CDTR',
            'C-MCC3-10 96 scan areas.xlsx': 'CDAU',
            'C-MRU224 96 scan areas.xlsx': 'CDAU',
            'C-B11220-AR381 96 scan areas.xlsx': 'CDAU',
            'C-B8441-AR387 96 scan areas.xlsx': 'CDAU',
            'C-ML315-AR1037 96 scan areas.xlsx': 'CDAU',
            'C-96299088 96 scan areas.xlsx': 'CDTR', ### id checken
            'C-95738583 96 scan areas.xlsx': 'CDGL', ### id checken
            'C-95690277 96 scan areas.xlsx': 'CDDU', ### id checken
            'C-CDKR-IQC 96 scan areas.xlsx': 'CDKR',
            'C-ATCC-CDGL 96 scan areas.xlsx': 'CDGL'}
