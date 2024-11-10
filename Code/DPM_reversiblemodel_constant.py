from DPM_reversiblemodel_lib import *
κ = 40
σ0, σhalf, σfull = 0, 0.5, 1
α_min, α_step, α_max = 0.0, 0.1, 1
α_default = [0.3, 0.4, 0.6, 0.7]
α_min_nomut, α_step_nomut, α_max_nomut = 0.0, 0.1, 1.0
α_min_norever, α_step_norever, α_max_norever = 0.0, 0.1, 0.0

θ_step, θ_max = 0.1, 1.
θ_step_nomut = 0.1
θ_max_norever = 0

g_full = [0.001, 0.002642, 0.00698, 0.018439, 0.048714, 0.128696, 0.34]
g_norever = [0.001, 0.002642, 0.00698, 0.018439, 0.048714, 0.128696, 0.34]
g_nomut = [0.001, 0.002642, 0.00698, 0.018439, 0.048714, 0.128696, 0.34]

S = [0.000560, 0.005379, 0.051674, 0.496387, 4.768310, 45.804544, 440]
S_nomut = [0.000560, 0.005379, 0.051674, 0.496387, 4.768310, 45.804544, 440]
S_norever = [0.000560, 0.005379, 0.051674, 0.496387, 4.768310, 45.804544, 440]

SR = [9.146e-4, 8.747e-3, 8.365e-2, 0.8]
SR_mutonly = [0, 1e-5, 9.564e-5, 9.146e-4, 8.747e-3, 8.365e-2, 0.8]
SR_reveronly = [0, 1e-5, 9.564e-5, 9.146e-4, 8.747e-3, 8.365e-2, 0.8]

T = [1e-7, 2.154e-6, 4.642e-5, 1e-3]
T_nomut = [0]

R_ratio = [1e-7, 1e-5, 1e-3, 0.1, 0.9]
R_ratio_nomut = [0]

PAR_valpool = {'S_01_2/S_01_1': [0.0004, 0.001474, 0.005429, 0.02, 0.073681, 0.271442, 1]}

par_save_block_size = int(1e4)
pathsave_parraw = './reversible_parraw'
pathsave_par = './reversible_par'
pathsave_nomut_par = './reversible_nomut_par'
pathsave_norever_par = './reversible_norever_par'
rate_key = ['Storr1_σ0',
            'rr1toS_σ0',
            'Storr2_σ0',
            'rr2toS_σ0',
            'Storr1_σfull',
            'rr1toS_σfull',
            'Storr2_σfull',
            'rr2toS_σfull']
par_key = ['α1', 'θ1', 'κ1', 'μ1',
           'α2', 'θ2', 'κ2', 'μ2',
           'g',
           'R1ratio',
           'R2ratio',
           'T1',
           'T2',
           'S_01_1',
           'S_01_2',
           'S_00_1',
           'S_00_2',
           'S_11_1',
           'S_11_2']
pd_colname = ['alpha1', 'theta1', 'kappa1', 'mu1',
              'alpha2', 'theta2', 'kappa2', 'mu2',
              'g',
              'R1ratio',
              'R2ratio',
              'T1',
              'T2',
              'S_01_1',
              'S_01_2',
              'S_00_1',
              'S_00_2',
              'S_11_1',
              'S_11_2']

par_key_all = par_key + rate_key
pd_colname_full = ['g',
                   'R1ratio',
                   'R2ratio',
                   'T1',
                   'T2',
                   'S_01_1',
                   'S_01_2',
                   'S_00_1',
                   'S_00_2',
                   'S_11_1',
                   'S_11_2',
                   'alpha1', 'theta1', 'kappa1', 'mu1',
                   'alpha2', 'theta2', 'kappa2', 'mu2']

par_sensitivity = ['alpha1',         # 0
                   'theta1',         # 1
                   'mu1',            # 2
                   'alpha2',         # 3
                   'theta2',         # 4
                   'mu2',            # 5
                   'T1',             # 6
                   'T2',             # 7
                   'g',              # 8
                   'S_01_1',         # 9
                   'S_01_2_S_01_1',  # 10
                   'SR_00_1',        # 11
                   'SR_00_2',        # 12
                   'SR_11_1',        # 13
                   'SR_11_2',        # 14
                   'R1ratio',        # 15
                   'R2ratio']        # 16

par_sensitivity_use2 = ['T1',
                        'T2',
                        'S_01_2_S_01_1',
                        'SR_00_1',
                        'SR_11_1',
                        'SR_11_2',
                        'R1ratio',
                        'R2ratio']

par_sensitivity2 = ['αθμ']
par_sensitivity_index_same = [6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16]

par_sensitivity3 = ['T1', 'T2', 'SR_00_1', 'SR_11_1', 'SR_11_2']
par_sensitivity_index3 = [[6], [7], [11, 12], [13], [14]]
par_sensitivity_index = [[0, 1, 2, 3, 4, 5], [6], [7], [10], [11, 12], [13], [14], [15], [16]]


block_format = '{0:0' + str(int(5)) + '}'
#            0       1       2       3       4       5       6       7       8
celltype = ['0101', '0001', '0100', '0000', '0111', '0011', '1101', '1100', '1111']
ind_0101, ind_0001, ind_0100, ind_0000, ind_0111, ind_0011, ind_1101, ind_1100, ind_1111 = 0, 1, 2, 3, 4, 5, 6, 7, 8
ind_mutonly = [ind_0101, ind_0111, ind_1101, ind_1111]
ind_reveronly = [ind_0101, ind_0001, ind_0100, ind_0000]
ind_rever = [ind_0001, ind_0100, ind_0000, ind_0011, ind_1100]
ind_nonmut = [0, 1, 2, 3]
ind_mut = [4, 5, 6, 7, 8]
ind_R1mut = [6, 7]
ind_R2mut = [4, 5]
ind_multi_mut = 8
Num_type = 4  # nonmut, R1 mut, R2 mut, R12 multi mut
Num_drug = 2
Limit_mortality = 1e13
Limit_radiologicdetection = 1e9
Stepsize = 42
duration_5year = 1820
daysperweek = 7
Numstep = math.ceil(duration_5year/Stepsize)
cellflow = ['n0101_pro',            # 0
            'n0001_pro',            # 1
            'n0100_pro',            # 2
            'n0000_pro',            # 3
            'n0111_pro',            # 4
            'n0011_pro',            # 5
            'n1101_pro',            # 6
            'n1100_pro',            # 7
            'n1111_pro',            # 8
            'n0101_to_n0001',       # 9
            'n0101_to_n0100',       # 10
            'n0001_to_n0101',       # 11
            'n0100_to_n0101',       # 12
            'n0001_to_n0000',       # 13
            'n0100_to_n0000',       # 14
            'n0000_to_n0001',       # 15
            'n0000_to_n0100',       # 16
            'n1101_to_n1100',       # 17
            'n1100_to_n1101',       # 18
            'n0111_to_n0011',       # 19
            'n0011_to_n0111',       # 20
            'n0101_to_n1101',       # 21
            'n0101_to_n0111',       # 22
            'n0001_to_n1101',       # 23
            'n0001_to_n0011',       # 24
            'n0100_to_n0111',       # 25
            'n0100_to_n1100',       # 26
            'n0000_to_n1100',       # 27
            'n0000_to_n0011',       # 28
            'n1101_to_n1111',       # 29
            'n0111_to_n1111',       # 30
            'n1100_to_n1111',       # 31
            'n0011_to_n1111']       # 32
cellflow_pro = [0, 1, 2, 3, 4, 5, 6, 7, 8]
cellflow_tran = [9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20]
cellflow_muta = [21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32]

duration_steady = 365  # one year
Simtimestep = 1
n_total = 5e9
ODE_method = 'Radau'
SCIPY_INTEGRATE_SOLVE_IVP = {'RTOL': 1e-3, 'ATOL': 1e-3, 'max_step': 1}

n0101_init, n0001_init, n0100_init, n0000_init = 10., 0., 0., 0.
n0111_init, n0011_init, n1101_init, n1100_init, n1111_init = 0., 0., 0., 0., 0.
x0_init = list([n0101_init, n0001_init, n0100_init, n0000_init, n0111_init, n0011_init, n1101_init, n1100_init, n1111_init])

σ_notreat = (0, 0)
σ1_full = (1, 0)
σ2_full = (0, 1)
σ_half = (0.5, 0.5)
dose_combination = [σ1_full, σ2_full, σ_half]

strategy_default = ['mono drug1',
                    'mono drug2',
                    'mono half',
                    'best cycle',
                    'strategy0',
                    'strategy1',
                    'strategy1 cycle',
                    'strategy2',
                    'strategy2 cycle']

nomut_cycleweek = np.arange(1, 27, 1)
strategy_nomut = strategy_default[:3] + [str(i) + 'w ' + i_σ for i in nomut_cycleweek for i_σ in ['σ1 first', 'σ2 first']]
ind_strategy_nocycle = [i for i, i_strategy in enumerate(strategy_default) if 'cycle' not in i_strategy]

strategy2threshold = 1e11
max_singledrugduration = 52  # 52 weeks, one year
weekpercycle = [3, 6]
strategy_cycle_weekpercycle = list(np.array(weekpercycle).astype('int').astype('str'))
strategy_cycle_weekpercycle = [i_weekpercycle + 'w' for i_weekpercycle in strategy_cycle_weekpercycle]
strategy_cycle_name = [i_weekpercycle+' '+i_first for i_first in ['σ1 1th', 'σ2 1th'] for i_weekpercycle in strategy_cycle_weekpercycle]

Num_week = Stepsize/7
weekpercycle_DPM = []
for i_Weekpercycle in range(1, int(Num_week/2) + 1):
    if Num_week % (2 * i_Weekpercycle) == 0:
        weekpercycle_DPM.append(i_Weekpercycle)
weekpercycle_DPM = [3]   # [3, 6]

teatmentname = ['σ1 full', 'σ2 full', 'half', '6w σ1 1th', '6w σ2 1th', '3w σ1 1th', '3w σ2 1th']
teatmentnamerandom = ['σ1 full', 'σ2 full', 'half', '3w σ1 1th', '3w σ2 1th']
Heading_param_mutonly = ['paramID'] + pd_colname
Heading_param_full = ['paramID'] + pd_colname_full
Heading_dosage = ['paramID', 'Strategy name']
Heading_dosage.extend([f'drug1, drug2 at t={t_i}' for t_i in np.arange(0, duration_5year, Stepsize)])
Heading_stopt = ['paramID'] + strategy_default
Heading_stopt_nomut = ['paramID'] + [i.replace('σ', 'sigma') for i in strategy_nomut]


Ntotal_ratio_plot = 1.5
MINN_plot = .1
Ylim_add_plot = 2
color_drug_plot = ['#98F5FF', '#FFEC8B']
Ylimmin_plot = 1e-6
Xlimmin_plot = -1
Lable_drug = ('Drug 1', 'Drug 2')
# Legend_colnum = 6
# Legend_order = np.reshape(np.array(range(0, 6)), (-1, Legend_colnum)).flatten('F')
# Legend_boxanchor = (0.5, 1.15)

color33 = [(0.036188642032125906, 0.9970485376113385, 0.17062970257010457),
           (0.9572745069684092, 0.0026458856838154077, 0.8947132352989293),
           (0.0, 0.5, 1.0), (1.0, 0.5, 0.0),
           (0.4650786908927126, 0.28911543467877776, 0.4651251824393656),
           (0.7932912084606683, 0.9760764054747108, 0.2742753143314446),
           (0.16034308790342755, 0.9905065495673819, 0.7926128496933004),
           (0.7360996736431331, 0.5368648196244276, 0.9653251592002877),
           (0.7484884279650402, 0.01785175530855554, 0.0001472118174313808),
           (0.0722252182964418, 0.49753381359827553, 0.013455046067000542),
           (0.0, 0.0, 1.0), (0.0, 0.0, 0.5),
           (0.2674871447609737, 0.6741114547661433, 0.42667732903103983),
           (0.5134503055831275, 0.17007948498248915, 0.9712254609354182),
           (0.954189560309339, 0.33786797930460233, 0.49785537376047284),
           (0.6319756718712118, 0.8650077989058742, 0.7630786882823456),
           (0.5556769834755146, 0.615800119805989, 0.0007521762773736729),
           (0.9792961988274446, 0.6933932055153016, 0.5186869529783017),
           (0.007916064401007228, 0.38989726477481457, 0.524396358200471),
           (0.808475888375974, 0.03255926554107724, 0.48196970251731464),
           (0.4034440223161161, 0.23424457962924616, 0.02451197095279478),
           (0.39765113895930704, 0.9794004364201219, 0.44550410045212685),
           (0.4257302371687287, 0.9426136898688422, 0.0034006254793605972),
           (0.37513821104438705, 0.590519308790788, 0.7933297797236271),
           (0.6190970699028584, 0.5738435462491948, 0.49860007827406894),
           (0.44204279461280627, 0.005640086337699746, 0.3001790110701217),
           (0.3597622071361718, 0.016381438816429394, 0.6668172170887553),
           (0.08968390856245623, 0.251345807820759, 0.8074066156339954),
           (0.9683795014161362, 0.305627077459085, 0.9504645667010929),
           (0.029406140806420122, 0.24105186037094328, 0.21447063961061563),
           (0.987383455987808, 0.8229263375051603, 0.013278263931654322),
           (0.694493751649723, 0.35232067194883365, 0.18875637627816422),
           (0.0, 1.0, 0.5)]  # distinctipy.get_colors(33)

color9 = ['k', 'b', 'g', 'r', 'c', 'm', 'y', 'blueviolet', '#808080']
color11 = ['black', 'gray', 'darkorange', 'red', 'darkgreen', 'deepskyblue', 'blueviolet', 'gold', 'c', 'm', 'y']
color13 = ['black', 'gray', 'darkorange', 'red', 'darkgreen', 'deepskyblue', 'blueviolet', 'gold', 'c', 'm', 'y', '#DAA520', '#808080']
