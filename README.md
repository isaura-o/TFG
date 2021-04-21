# TFG
Bachelor degree Project

**Title**: Comorbidity in epidemic spreading diseases.

To run the programs via terminal:

    * Erdos-Reyni network: ./net_gen_v1  (Average degree) (Number of nodes. ex: 1234)  (seed for random values. ex: 2)
    
    * Scale-free network: ./conf_net_gen_v4 (Average degree) (Gamma= degree exponent. ex: 2.3) (Number of nodes. ex. 1234) (seed. ex: 2)
    
    * Network cleaner: ./net_cleaner_v3 (network file we want to clean: erdos-reyni or scale-free networks)

    * Network statistics: ./mesur_v5 (network file we want the statistics of)
                             !No network file cleaned needed since the program has the net_cleaner program as subroutine.

    * SIS-SIR dynamics file:  ./sisAB_v3  (network file previously cleaned) (recovery rate A, (dA)) (recovery rate B, (dB)) (infection rate A, (lA)) (infection rate B, (lB))
      
