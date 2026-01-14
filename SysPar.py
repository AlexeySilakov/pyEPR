import numpy as np

class sysPar():
    def __init__(self, **kwargs):
       # a pool of Easyspin-inspired Sys variables.
       # Basic spin definitions
       self.S = None
       self.Nucs = None
       self.useFor = None
       self.nNucs = None
       self.gn = None
       self.I  = None
       self.NatAbund = None # in %
       self.gamma = None # in MHz per mT
       # Electron Zeeman interaction
       self.g = None
       self.gFrame = None
       # Hyperfine interaction
       self.A = None # in MHz
       self.AFrame = None # in degree
       self.Apa = None # in degree
       # Nuclear quadrupole
       self.Q = None # in MHz
       self.Q_K = None
       self.Q_eta = None
       self.QFrame = None # in degree
       self.Qpa = None # in degree
       # Zero-field splitting
       self.D = None
       self.D_ = None
       self.DFrame = None
       # High-order Stevens operators
       self.Bk = {}         # e.g. {2: [...], 4: [...]}
       self.BkFrame = {}    # e.g. {2: (α, β, γ)}
       # Electron–electron couplings
       self.J = None
       self.dip = None
       self.dvec = None
       self.ee = None
       self.eeFrame = None
       self.ee2 = None
       # Orbital / angular momentum interactions
       self.L = None
       self.gL = None
       self.soc = None
       self.CF = {}         # Crystal field coefficients: CF0, CF2...CF12
       # Nuclear shielding
       self.sigma = None
       self.sigmaFrame = None
       # Broadening / strain parameters
       self.lw = None
       self.lwEndor = None
       self.gStrain = None
       self.AStrain = None
       self.DStrain = None
       self.HStrain = None
       self.RelativeScale = None

       # Override defaults with any passed keyword arguments
       for key, value in kwargs.items():
           if hasattr(self, key):
               setattr(self, key, value)
               if key == 'Nucs': self.setNucs(value)
           else:
               raise AttributeError(f"😭 Unknown parameter '{key}' for SysPar")

        
    def get(self, param):
        return getattr(self, param)
    def getToolTips(self, param=None):
        dc = {'Nucs':'need isotope number and atom letters such as 57Fe. Sys.gn and Sys.I are then pulled from the database',
              'nNucs': 'number of equivalent nuclei of this type to consider',
              'useFor': 'flags to simulation program to either only use this nucleus for simulation, orientation selection or both',
              'gamma': 'nuclear g-value in MHz per mT',
              'Apa': 'Euler angles in degree to rotate A-tensor to g-tensor frame (old Easyspin notation)',
              'AFrame': 'Euler angles in degree to rotate  g-tensor to A-tensor frame (new Easyspin notation)'}
        if param is None:
            return dc
        else:
            if param in dc.keys():
                return dc[param]
            else:
                return ''

    def set(self, **kwargs):
        for key, value in kwargs.items():
            if hasattr(self, key):
                setattr(self, key, value)
                if key == 'Nucs': self.setNucs(value)
            else:
                raise AttributeError(f"😭 Unknown parameter '{key}' for SysPar")
    def setDict(self, Dict):
        if type(Dict)==dict:
            for key in Dict.keys():
                if hasattr(self, key):
                    setattr(self, key, Dict[key])
                    if key == 'Nucs': self.setNucs(Dict[key]);
                else:
                    raise AttributeError(f"😭 Unknown parameter '{key}' for SysPar")
    def setNucs(self, Nucs=None, nNucs=None):
        self.gn,self.I,self.gamma=self.isotopes(Nucs)
        self.Nucs = Nucs
        if type(self.A)==type(None):
            self.A = []
            self.Apa = []
            for ii in range(len(Nucs)):
                self.A.append(np.zeros(3))
                self.Apa.append(np.zeros(3))
        elif type(self.A)==np.ndarray:
            raise AttributeError("🤔 HF coupling must be a list of np.array")
        
        if type(self.Q_K)!=type(None):
            if len(self.Q_K)==len(Nucs):
                if type(self.Q_eta)!=type(None):
                    Q_eta = [0.0]*len(Nucs)
                self.Q = []
                for ii in range(len(Nucs)):
                    self.Q.append(np.array([(1.0+Q_eta[ii])*self.Q_K[ii], 
                                            (1.0-Q_eta[ii])*self.Q_K[ii], 
                                            -2.0*self.Q_K[ii]]))
            else:
                raise AttributeError("🤔 lenght of Q_K is different from the lenght of Nucs")
        if type(self.Qpa)==type(None):
            self.Qpa = []
            for ii in range(len(Nucs)): self.Qpa.append(np.zeros(3))
        elif type(self.Q)==np.ndarray:
            raise AttributeError("🤔 Qpa  must be a list of np.array")        

        if type(self.Q)==type(None):
            self.Q = []
            for ii in range(len(Nucs)): self.Q.append(np.zeros(3))
        elif type(self.Q)==np.ndarray:
            raise AttributeError("🤔 Q  must be a list of np.array")

        if type(self.useFor)==type(None):
            self.useFor = []
            for ii in range(len(Nucs)): self.useFor.append('both')

        if not(self.nNucs) and (not nNucs):
            self.nNucs = np.ones_like(self.gn)
        elif nNucs:
            self.nNucs = nNucs
        
    def getDefaultDict(self):
        return {
                'spin(1)':{
                'S':1/2, 'g':np.array([2, 2, 2],dtype=float),
                'lw':0.1,
                },
                'nuc(1)':{
                'Nucs':['1H', list(self.isotopes().keys())],
                'A': np.array([1, 1, -2],dtype=float),
                'Apa':np.array([0, 0, 0],dtype=float),
                'Q_K':0.0,
                'Q_eta':0.0,
                'Qpa':np.array([0, 0, 0],dtype=float),
                'useFor':['all', ['sim.only', 'ori.sel.only', 'all']]
                },
                }
    def getDefaultDictEPR(self):
        return {
                'spin(1)':{
                'S':1/2, 'g':np.array([2, 2, 2],dtype=float),
                'lw':0.1,
                'RelativeScale': 1,
                },
                'nuc(1)':{
                'Nucs':['1H', list(self.isotopes().keys())],
                'A': np.array([1, 1, -2],dtype=float),
                'Apa':np.array([0, 0, 0],dtype=float),
                'Q_K':0.0,
                'Q_eta':0.0,
                'Qpa':np.array([0, 0, 0],dtype=float),
                'useFor':['all', ['sim.only', 'ori.sel.only', 'all']]
                },
                }    
    def setFromCtrl(self, inDict):
        self.set_All_None()
        nSpins = 0
        nNucs = 0
        spDic = []
        nucDic = []
        othDic = []
        for key in inDict.keys():
            if 'spin' in key:
                spDic.append(key)
            elif 'nuc' in key:
                nucDic.append(key)
            else:
                othDic.append(key)
        
        nSpins = len(spDic)
        nNucs = len(nucDic)
        nOthers = len(othDic)
        
        for key in inDict.keys():
            if key=='functions':
                continue
            
            if type(inDict[key])!=dict: 
                
                raise AttributeError("😭 primary input dic key is expected to be a dictionary")
            try:            
                inval = int(key.split('(')[1].split(')')[0])-1
            except Exception:
                raise AttributeError(f"😭 unexpected property tag {key}. We always expect name(number) format")
            if 'spin' in key:
                nPars = nSpins
            elif 'nuc' in key:
                nPars = nNucs
            else: 
                nPars = nOthers   
                
            for sk in inDict[key].keys():

                if not hasattr(self, sk):
                    raise AttributeError(f"😭 property {sk} is not a part of the Sys class")
                attr = getattr(self, sk)

                if type(attr)==type(None):
                    attr = [0]*nPars

                if type(inDict[key][sk])==list: ### from control it may come with choices as a second part of the list
                    attr[inval]=inDict[key][sk][0]
                else:
                    attr[inval]=inDict[key][sk]
                setattr(self, sk, attr)
        if type(self.Nucs)!=type(None):
            self.setNucs(self.Nucs)
        
                
    def getAll(self):
        return {'S': 'S_float', 'lw': 'S_float', 
         'D': 'S_array3', 'D_': 'S_array3', 'DFrame': 'S_array3',
         'Bk': 'S_array', 'BkFrame': 'S_array3', 
         'J': 'SS_array3', 'dip': 'SS_float', 
         'dvec': 'SS_array3', 
         'ee': 'SS_array3', 'eeFrame': 'SS_array3', 
         'ee2': 'SS_array3',
         'Nucs': 'Nucs_string',
         'nNucs': 'Nucs_int', 'gn': 'Nucs_float', 
         'I': 'Nucs_float', 'NatAbund': 'Nucs_float', 
         'gamma': 'Nucs_float', 'g': 'Nucs_array3',
         'gFrame': 'Nucs_array3', 'A': 'Nucs_array3', 
         'AFrame': 'Nucs_array3', 'Apa': 'Nucs_array3',
         'Q': 'Nucs_array3', 'QFrame': 'Nucs_array3', 
         'Qpa': 'Nucs_array3', 
         'Q_K':'Nucs_float', 'Q_eta': 'Nucs_float',
         'lwEndor': 'Nucs_float', 
         'gStrain': 'S_array3', 'AStrain': 'Nucs_array3', 
         'DStrain': 'S_array3', 'HStrain': 'S_array3'}
         
         # 'L': None, 'gL': None, 'soc': None, 'CF': {}, 'sigma': None, 'sigmaFrame': None, 'nn': None, 'nnFrame': None, 'Ham': {}, 
         # , 'initState': None}
    def set_All_None(self):  
        allvars = self.getAll()
        for key in allvars.keys():
            if hasattr(self, key):
                setattr(self, key, None)


    def isotopes(self,Nucs=None):
#   For now generated by ChatGPT5. Most common ones were verified against BrukerAlmanach tables. Data sources:
# - CRC Handbook of Chemistry and Physics, 104th Ed. (2023–2024), Nuclear Moments Tables.
# - NIST CODATA Recommended Values of the Fundamental Physical Constants (2018/2022).
# - IAEA Nuclear Data Services: Live Chart of Nuclides.
# - IUPAC Commission on Isotopic Abundances and Atomic Weights.
# in the future to completely rewrite this table using original sourses to not allow for AI halucinations.
        isotopes_list={
"1H"	:{"abundance": 100	, "I": 1/2	,"gn":5.58569	, "Q": 0.000000 , "nu350mT":14.90},
"2H"	:{"abundance": 0.015, "I": 1	,"gn":0.85744	, "Q": 0.00288	, "nu350mT":2.29 },
"3He"	:{"abundance": 0.0001, "I": 31/2	,"gn":5.95792	, "Q": 0.000000 , "nu350mT":15.90},
"6Li"	:{"abundance": 7.42	, "I": 1	,"gn":0.82205	, "Q": -0.00064	, "nu350mT":2.19 },
"7Li"	:{"abundance": 92.6	, "I": 3/2	,"gn":2.17096	, "Q": -0.04000	, "nu350mT":5.79 },
"9Be"	:{"abundance": 100	, "I": 3/2	,"gn":-0.78500	, "Q": 0.05300	, "nu350mT":2.09 },
"11B"	:{"abundance": 80.4	, "I": 3/2	,"gn":1.79242	, "Q": 0.04000	, "nu350mT":4.78 },
"10B"	:{"abundance": 19.6	, "I": 3	,"gn":0.60022	, "Q": 0.08608	, "nu350mT":1.60 },
"13C"	:{"abundance": 1.11	, "I": 1/2	,"gn":1.40482	, "Q": 0.000000 , "nu350mT":3.75 },
"14N"	:{"abundance": 99.6	, "I": 1	,"gn":0.40378	, "Q": 0.01930	, "nu350mT":1.08 },
"15N"	:{"abundance": 0.365, "I": 1/2	,"gn":-0.56638	, "Q": 0.000000 , "nu350mT":1.51 },
"17O"	:{"abundance": 0.037, "I": 5/2	,"gn":-0.75752	, "Q": -0.02600	, "nu350mT":2.02 },
"19F"	:{"abundance": 100	, "I": 1/2	,"gn":5.25773	, "Q": 0.000000 , "nu350mT":14.03},
"21Ne"	:{"abundance": 0.257, "I": 3/2	,"gn":-0.44120	, "Q": 0.10290	, "nu350mT":1.18 },
"23Na"	:{"abundance": 100	, "I": 3/2	,"gn":1.47839	, "Q": 0.10800	, "nu350mT":3.94 },
"25Mg"	:{"abundance": 10.1	, "I": 5/2	,"gn":-0.34218	, "Q": 0.22000	, "nu350mT":0.91 },
"27Al"	:{"abundance": 100	, "I": 5/2	,"gn":1.45660	, "Q": 0.15000	, "nu350mT":3.89 },
"29Si"	:{"abundance": 4.7	, "I": 1/2	,"gn":-1.11060	, "Q": 0.000000 , "nu350mT":2.96 },
"31P"	:{"abundance": 100	, "I": 1/2	,"gn":2.26320	, "Q": 0.000000 , "nu350mT":6.04 },
"33S"	:{"abundance": 0.76	, "I": 3/2	,"gn":0.42911	, "Q": -0.06400	, "nu350mT":1.14 },
"35Cl"	:{"abundance": 75.5	, "I": 3/2	,"gn":0.54792	, "Q": -0.08249	, "nu350mT":1.46 },
"37Cl"	:{"abundance": 24.5	, "I": 3/2	,"gn":0.45608	, "Q": -0.06493	, "nu350mT":1.22 },
"39K"	:{"abundance": 93.1	, "I": 3/2	,"gn":0.26099	, "Q": 0.05400	, "nu350mT":0.70 },
"43Ca"	:{"abundance": 0.145, "I": 7/2	,"gn":-0.37641	, "Q": 0.000000 , "nu350mT":1.00 },
"45Sc"	:{"abundance": 100	, "I": 7/2	,"gn":1.35906	, "Q": -0.22000	, "nu350mT":3.63 },
"47Ti"	:{"abundance": 7.25	, "I": 5/2	,"gn":-0.31539	, "Q": 0.29000	, "nu350mT":0.84 },
"49Ti"	:{"abundance": 5.51	, "I": 7/2	,"gn":-0.31548	, "Q": 0.24000	, "nu350mT":0.84 },
"50V"	:{"abundance": 0.24	, "I": 6	,"gn":0.55659	, "Q": 0.20900	, "nu350mT":1.48 },
"51V"	:{"abundance": 99.8	, "I": 7/2	,"gn":1.46836	, "Q": -0.05150	, "nu350mT":3.92 },
"53Cr"	:{"abundance": 9.55	, "I": 3/2	,"gn":-0.31470	, "Q": -1.29545	, "nu350mT":0.84 },
"55Mn"	:{"abundance": 100	, "I": 5/2	,"gn":1.38190	, "Q": 0.33000	, "nu350mT":3.69 },
"57Fe"	:{"abundance": 2.19	, "I": 1/2	,"gn":0.18060	, "Q": 0.000000 , "nu350mT":0.48 },
"59Co"	:{"abundance": 100	, "I": 7/2	,"gn":1.31800	, "Q": 0.42000	, "nu350mT":3.52 },
"61Ni"	:{"abundance": 1.19	, "I": 3/2	,"gn":-0.50001	, "Q": 0.16200	, "nu350mT":1.33 },
"65Cu"	:{"abundance": 30	, "I": 3/2	,"gn":1.58800	, "Q": -0.19500	, "nu350mT":4.24 },
"63Cu"	:{"abundance": 69.1	, "I": 3/2	,"gn":1.48400	, "Q": -0.22200	, "nu350mT":3.96 },
"67Zn"	:{"abundance": 4.11	, "I": 5/2	,"gn":0.35031	, "Q": 0.15000	, "nu350mT":0.93 },
"71Ga"	:{"abundance": 39.6	, "I": 3/2	,"gn":1.70818	, "Q": 0.000000 , "nu350mT":4.56 },
"69Ga"	:{"abundance": 60.4	, "I": 3/2	,"gn":1.34439	, "Q": 0.000000 , "nu350mT":3.59 },
"73Ge"	:{"abundance": 7.76	, "I": 9/2	,"gn":-0.19544	, "Q": -0.19000	, "nu350mT":0.52 },
"75As"	:{"abundance": 100	, "I": 3/2	,"gn":0.95965	, "Q": 0.29000	, "nu350mT":2.56 },
"77Se"	:{"abundance": 5.12	, "I": 1/2	,"gn":1.06930	, "Q": 0.000000 , "nu350mT":2.85 },
"79Br"	:{"abundance": 50.5	, "I": 3/2	,"gn":1.40427	, "Q": 0.29300	, "nu350mT":3.75 },
"81Br"	:{"abundance": 49.5	, "I": 3/2	,"gn":1.51371	, "Q": 0.27000	, "nu350mT":4.04 },
"83Kr"	:{"abundance": 11.6	, "I": 9/2	,"gn":-0.21570	, "Q": 0.26000	, "nu350mT":0.58 },
"87Rb"	:{"abundance": 27.9	, "I": 3/2	,"gn":1.82760	, "Q": 0.13000	, "nu350mT":4.88 },
"85Rb"	:{"abundance": 72.2	, "I": 5/2	,"gn":0.53928	, "Q": 0.27300	, "nu350mT":1.44 },
"87Sr"	:{"abundance": 7.02	, "I": 9/2	,"gn":-0.24291	, "Q": 0.15000	, "nu350mT":0.65 },
"89Y"	:{"abundance": 100	, "I": 1/2	,"gn":-0.27364	, "Q": 0.000000 , "nu350mT":0.73 },
"91Zr"	:{"abundance": 11.2	, "I": 5/2	,"gn":-0.52145	, "Q": 0.000000 , "nu350mT":1.39 },
"93Nb"	:{"abundance": 100	, "I": 9/2	,"gn":1.37120	, "Q": -0.28000	, "nu350mT":3.66 },
"97Mo"	:{"abundance": 9.46	, "I": 5/2	,"gn":-0.37340	, "Q": 0.20000	, "nu350mT":1.00 },
"95Mo"	:{"abundance": 15.7	, "I": 5/2	,"gn":-0.36560	, "Q": -0.01900	, "nu350mT":0.98 },
"99Ru"	:{"abundance": 12.7	, "I": 3/2	,"gn":-0.24900	, "Q": 0.07600	, "nu350mT":0.66 },
"101Ru"	:{"abundance": 17.1	, "I": 5/2	,"gn":-0.27900	, "Q": 0.44000	, "nu350mT":0.74 },
"103Rh"	:{"abundance": 100	, "I": 1/2	,"gn":-0.17680	, "Q": 0.000000 , "nu350mT":0.47 },
"105Pd"	:{"abundance": 22.2	, "I": 5/2	,"gn":-0.25600	, "Q": 0.66000	, "nu350mT":0.68 },
"107Ag"	:{"abundance": 51.8	, "I": 1/2	,"gn":-0.22725	, "Q": 0.000000 , "nu350mT":0.61 },
"109Ag"	:{"abundance": 48.2	, "I": 1/2	,"gn":-0.26174	, "Q": 0.000000 , "nu350mT":0.70 },
"111Cd"	:{"abundance": 12.8	, "I": 1/2	,"gn":-1.19043	, "Q": 0.000000 , "nu350mT":3.18 },
"113Cd"	:{"abundance": 12.3	, "I": 1/2	,"gn":-1.24540	, "Q": 0.000000 , "nu350mT":3.32 },
"113In"	:{"abundance": 4.28	, "I": 9/2	,"gn":1.22864	, "Q": 0.84600	, "nu350mT":3.28 },
"115In"	:{"abundance": 95.7	, "I": 9/2	,"gn":1.23129	, "Q": 0.86100	, "nu350mT":3.29 },
"117Sn"	:{"abundance": 7.61	, "I": 1/2	,"gn":-2.00208	, "Q": 0.000000 , "nu350mT":5.34 },
"119Sn"	:{"abundance": 8.58	, "I": 1/2	,"gn":-2.09456	, "Q": 0.000000 , "nu350mT":5.59 },
"115Sn"	:{"abundance": 0.35	, "I": 1/2	,"gn":-1.83770	, "Q": 0.000000 , "nu350mT":4.90 },
"121Sb"	:{"abundance": 57.3	, "I": 5/2	,"gn":1.34550	, "Q": -0.33000	, "nu350mT":3.59 },
"123Sb"	:{"abundance": 42.8	, "I": 7/2	,"gn":0.72876	, "Q": -0.68000	, "nu350mT":1.94 },
"123Te"	:{"abundance": 0.87	, "I": 1/2	,"gn":-1.47360	, "Q": 0.000000 , "nu350mT":3.93 },
"125Te"	:{"abundance": 6.99	, "I": 1/2	,"gn":-1.77660	, "Q": 0.000000 , "nu350mT":4.74 },
"127I"	:{"abundance": 100	, "I": 5/2	,"gn":1.12530	, "Q": -0.78900	, "nu350mT":3.00 },
"131Xe"	:{"abundance": 21.2	, "I": 3/2	,"gn":0.46124	, "Q": -0.12000	, "nu350mT":1.23 },
"129Xe"	:{"abundance": 26.4	, "I": 1/2	,"gn":-1.55595	, "Q": 0.000000 , "nu350mT":4.15 },
"133Cs"	:{"abundance": 100	, "I": 7/2	,"gn":0.73785	, "Q": -0.00300	, "nu350mT":1.97 },
"135Ba"	:{"abundance": 6.59	, "I": 3/2	,"gn":0.55884	, "Q": 0.20000	, "nu350mT":1.49 },
"137Ba"	:{"abundance": 11.3	, "I": 3/2	,"gn":0.62515	, "Q": 0.34000	, "nu350mT":1.67 },
"139La"	:{"abundance": 99.9	, "I": 7/2	,"gn":0.79520	, "Q": 0.20000	, "nu350mT":2.12 },
"138La"	:{"abundance": 0.089, "I": 5	,"gn":0.74278	, "Q": 0.51000	, "nu350mT":1.98 },
"141Pr"	:{"abundance": 100	, "I": 5/2	,"gn":1.60000	, "Q": -0.04100	, "nu350mT":4.27 },
"143Nd"	:{"abundance": 12.2	, "I": 7/2	,"gn":-0.30760	, "Q": -0.56000	, "nu350mT":0.82 },
"145Nd"	:{"abundance": 8.3	, "I": 7/2	,"gn":-0.19000	, "Q": -0.29000	, "nu350mT":0.51 },
"147Sm"	:{"abundance": 15	, "I": 7/2	,"gn":-0.23220	, "Q": -0.18000	, "nu350mT":0.62 },
"149Sm"	:{"abundance": 13.8	, "I": 7/2	,"gn":0.19150	, "Q": 0.05600	, "nu350mT":0.51 },
"153Eu"	:{"abundance": 52.2	, "I": 5/2	,"gn":0.61340	, "Q": 3.92000	, "nu350mT":1.64 },
"151Eu"	:{"abundance": 47.8	, "I": 5/2	,"gn":1.38900	, "Q": 1.53000	, "nu350mT":3.71 },
"157Gd"	:{"abundance": 15.7	, "I": 3/2	,"gn":-0.22530	, "Q": 1.34000	, "nu350mT":0.60 },
"155Gd"	:{"abundance": 14.7	, "I": 3/2	,"gn":-0.17230	, "Q": 1.30000	, "nu350mT":0.46 },
"159Tb"	:{"abundance": 100	, "I": 3/2	,"gn":1.34200	, "Q": 1.34000	, "nu350mT":3.58 },
"161Dy"	:{"abundance": 18.9	, "I": 5/2	,"gn":-0.18900	, "Q": 2.47000	, "nu350mT":0.50 },
"163Dy"	:{"abundance": 25	, "I": 5/2	,"gn":0.26600	, "Q": 2.51000	, "nu350mT":0.71 },
"165Ho"	:{"abundance": 100	, "I": 7/2	,"gn":1.19200	, "Q": 2.73000	, "nu350mT":3.18 },
"167Er"	:{"abundance": 22.9	, "I": 7/2	,"gn":-0.16180	, "Q": 2.82700	, "nu350mT":0.43 },
"169Tm"	:{"abundance": 100	, "I": 1/2	,"gn":-0.46600	, "Q": 0.00000 	, "nu350mT":1.24 },
"173Yb"	:{"abundance": 16.1	, "I": 5/2	,"gn":-0.27195	, "Q": 2.80000	, "nu350mT":0.73 },
"171Yb"	:{"abundance": 14.3	, "I": 1/2	,"gn":0.98850	, "Q": 0.000000	, "nu350mT":2.64 },
"175Lu"	:{"abundance": 37.4	, "I": 7/2	,"gn":0.63943	, "Q": 5.68000	, "nu350mT":1.71 },
"179Hf"	:{"abundance": 13.8	, "I": 9/2	,"gn":-0.14240	, "Q": 5.10000	, "nu350mT":0.38 },
"177Hf"	:{"abundance": 18.5	, "I": 7/2	,"gn":0.22670	, "Q": 4.50000	, "nu350mT":0.60 },
"181Ta"	:{"abundance": 100	, "I": 7/2	,"gn":0.67729	, "Q": 3.44000	, "nu350mT":1.81 },
"183W"	:{"abundance": 14.4	, "I": 1/2	,"gn":0.23557	, "Q": 0.000000	, "nu350mT":0.63 },
"187Re"	:{"abundance": 62.9	, "I": 5/2	,"gn":1.28780	, "Q": 2.22000	, "nu350mT":3.44 },
"185Re"	:{"abundance": 37.1	, "I": 5/2	,"gn":1.27480	, "Q": 2.33000	, "nu350mT":3.40 },
"187Os"	:{"abundance": 1.64	, "I": 1/2	,"gn":0.13110	, "Q": 0.000000 , "nu350mT":0.35 },
"189Os"	:{"abundance": 16.1	, "I": 3/2	,"gn":0.48800	, "Q": 0.80000	, "nu350mT":1.30 },
"191Ir"	:{"abundance": 37.3	, "I": 3/2	,"gn":0.09700	, "Q": 0.78000	, "nu350mT":0.26 },
"193Ir"	:{"abundance": 62.7	, "I": 3/2	,"gn":0.10700	, "Q": 0.70000	, "nu350mT":0.29 },
"195Pt"	:{"abundance": 34	, "I": 1/2	,"gn":1.21900	, "Q": 0.000000 , "nu350mT":3.25 },
"197Au"	:{"abundance": 100	, "I": 3/2	,"gn":0.09797	, "Q": 0.59400	, "nu350mT":0.26 },
"201Hg"	:{"abundance": 13.2	, "I": 3/2	,"gn":-0.37348	, "Q": 0.42000	, "nu350mT":1.00 },
"199Hg"	:{"abundance": 16.8	, "I": 1/2	,"gn":1.01177	, "Q": 0.000000 , "nu350mT":2.70 },
"205Tl"	:{"abundance": 70.5	, "I": 1/2	,"gn":3.27540	, "Q": 0.000000 , "nu350mT":8.74 },
"203Tl"	:{"abundance": 29.5	, "I": 1/2	,"gn":3.24451	, "Q": 0.000000 , "nu350mT":8.66 },
"207Pb"	:{"abundance": 22.6	, "I": 1/2	,"gn":1.17480	, "Q": 0.000000 , "nu350mT":3.13 },
"209Bi"	:{"abundance": 100	, "I": 9/2	,"gn":0.93800	, "Q": -0.46000	, "nu350mT":2.50 },
"229Th"	:{"abundance": 0	, "I": 5/2	,"gn":0.16000	, "Q": 4.40000	, "nu350mT":0.43 },
"235U"	:{"abundance": 0.72	, "I": 7/2	,"gn":-0.10000	, "Q": 4.30000	, "nu350mT":0.27 },
 }



        if Nucs==None:
            return isotopes_list
        else:
            gn=[]
            I = []
            gamma = []
            if isinstance(Nucs, list):
                for nn in Nucs:
                    if nn in isotopes_list:
                        gn.append(isotopes_list[nn]["gn"])
                        I.append(isotopes_list[nn]["I"])
                        gamma.append(isotopes_list[nn]["nu350mT"]/350.0)
                    else:
                        raise AttributeError(f"Unknown nuclei '{nn}' for sysPar")
            if isinstance(Nucs, str):
                if Nucs in isotopes_list:
                    gn.append(isotopes_list[Nucs]["gn"])
                    I.append(isotopes_list[Nucs]["I"])
                    gamma.append(isotopes_list[Nucs]["nu350mT"]/350.0)
                else:
                    raise AttributeError(f"Unknown nuclei '{Nucs}' for sysPar")
            return gn,I,gamma
        
class expPar():
    def __init__(self, **kwargs):
        # Question - isn't it safer to have all values to be None, so that used does not accidentally forget something
        self.mwFreq = 0 # MHz
        self.nPoints = 256 # MHz
        self.MaxFreq = 10 # MHz
        self.gField = None
        self.Field = 350
        self.BMin = 300 # mT
        self.BMax = 400 # mT
        self.tau = 120 # ns
        self.Tinv = 40 # ns
        self.ExciteWidth = 0 # MHz
        self.phi = None # for single crystal
        self.theta = None # for single crystal
        self.tDead=None
        self.Harmonic = int(1)
        for key, value in kwargs.items():
            if hasattr(self, key):
                setattr(self, key, value)
                # if key == 'Nucs': self.setNucs(value);
            else:
                raise AttributeError(f"Unknown parameter '{key}' for expPar")
    def set(self, **kwargs):
        for key, value in kwargs.items():
            if hasattr(self, key):
                setattr(self, key, value)
            else:
                raise AttributeError(f"Unknown parameter '{key}' for expPar")
    def getToolTips(self, param=None):
        ss = 'If "-1", take from data when available'
        dc = {'mwFreq':'Microwave Frequnecy in [GHz]. '+ss,
              'nPoints': 'number of points to calculate. '+ss,
              'Field':'Magnetic Field in [mT]. '+ss,
              'gField':'Magnetic Field in g-value. Actual Field calculated from mwFreq. ',
              'BMin':'Start magnetic Field in [mT]. '+ss,
              'BMax':'End magnetic Field in [mT]. '+ss,
              'tau':'tau, typically a delay between first two MW pulses [ns]. ',
              'ExciteWidth': 'Excitation width in [MHz]',
              }
        if param is None:
            return dc
        else:
            if param in dc.keys():
                return dc[param]
            else:
                return ''

    def getDefaultDict(self):
        return {'mwFreq': float(-1.0), # MHz
                'nPoints': int(256), # MHz
                'MaxFreq': float(-1), # MHz
                'Field': float(-1.0),
                'tau': float(-1.0), # ns
                'ExciteWidth': float(20.0), # MHz
                'tDead':float(0.0),
                }
    def getAll(self):
        return {'mwFreq': -1.0, # MHz
                'nPoints': 256, # MHz
                'MaxFreq': 10.0, # MHz
                'gField': -1.0,
                'Field': 350.0,
                'BMin': 300.0, # mT
                'BMax': 400.0, # mT
                'tau': 120.0, # ns
                'Tinv': 40.0, # ns
                'ExciteWidth': 0.0, # MHz
                'phi': None, # for single crystal
                'theta': None, # for single crystal
                'tDead':None, 
                'Harmonic':1}
    def setDict(self, Dict):
        if type(Dict)==dict:
            for key in Dict.keys():
                if hasattr(self, key):
                    setattr(self, key, Dict[key])
                else:
                    raise AttributeError(f"😭 Unknown parameter '{key}' for expPar")
    
    def setFromCtrl(self, inDict):
        self.setDict(inDict)
        ## just because the other one has it. I don't think there is anything to do here. 
            
# ------------------------------------------------------------------
if __name__ == "__main__":          
    Sys = sysPar(S=[1/2])
    isotopes = Sys.isotopes()
    Sys.setNucs(['1H', '2H'])
    Exp = expPar(mwFreq = 34.6)
    print(Sys.__dict__)