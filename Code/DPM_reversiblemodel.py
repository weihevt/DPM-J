from DPM_reversiblemodel_fun import *
from DPM_reversiblemodel_plot import *

'''Run simulation for the reversible resistance only model'''
DPM_reversiblemodel_fun_run_csvfolder(pathload='./reversible_parcsv_rever/',
                                      pathsave='./reversible_csvresult_rever/',
                                      Strategy_name=strategy_default[:7])

'''Run simulation for the irreversible resistance only model'''
DPM_reversiblemodel_fun_run_csvfolder(pathload='./reversible_parcsv_mut',
                                      pathsave='./reversible_csvresult_mut',
                                      Strategy_name=strategy_default,
                                      mutonly=True)

'''Run simulation for the joint model'''
DPM_reversiblemodel_fun_run_csvfolder(pathload='./reversible_parcsv/',
                                      pathsave='./reversible_csvresult/',
                                      Strategy_name=strategy_default)


'''Analysis of simulation results from the reversible resistance only model'''
pathload = './reversible_csvresult_rever/'
DPM_reversiblemodel_fun_ayalysis(pathload=pathload,
                                 strategy_name=strategy_default[:7],
                                 mutonly=False, reveronly=True)

'''Analysis of simulation results from the irreversible resistance only model'''
pathload = './reversible_csvresult_mut/'
DPM_reversiblemodel_fun_ayalysis(pathload=pathload,
                                 strategy_name=strategy_default,
                                 mutonly=True, reveronly=False)

'''Analysis of simulation reuslts from the joint model'''
pathload = './reversible_csvresult/'
DPM_reversiblemodel_fun_ayalysis(pathload=pathload,
                                 strategy_name=strategy_default,
                                 mutonly=False, reveronly=False)
