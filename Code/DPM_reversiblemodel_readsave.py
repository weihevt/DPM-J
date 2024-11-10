from DPM_reversiblemodel_constant import *
from DPM_reversiblemodel_miscellaneous import DPM_reversiblemodel_miscellaneous_findtreat, DPM_reversiblemodel_stime


def DPM_reversiblemodel_readsave_save_csv(PAR, result, filename_param, filename_stopt, filename_dosage, filename_pop, reveronly=False):
    PAR_index = str(int(PAR['paramID']))
    PAR_val = list(PAR.values())
    PAR_val = PAR_val[:len(pd_colname) + 1]
    PAR_str = ','.join([str(i) for i in PAR_val])
    if filename_param:
        with open(filename_param, 'a', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(PAR_str.split(','))

    cured_survival = duration_5year + Stepsize
    non_report = 'nan'

    strategy_result, survival = dict(zip(strategy_default, {'stopt': None, 'dose': None})), []
    if not reveronly:
        for i, i_strategy in enumerate(strategy_default):
            strategy_result[i_strategy] = DPM_reversiblemodel_readsave_assign(
                result[i], strategy_default[i], cured_survival, non_report, PAR_index, filename_stopt, filename_dosage, filename_pop)
    else:
        for i_N in result:
            if all(i_N[:, -1] < 1):
                survival_strategy = cured_survival
            else:
                survival_strategy = DPM_reversiblemodel_stime(i_N)
                if survival_strategy == duration_5year and sum(i_N[:, -1]) < Limit_mortality:
                    survival_strategy += 1
            survival.append(survival_strategy)

    if filename_stopt:
        survival = [strategy_result[i_strategy]['stopt'] for i_strategy in strategy_default]
        survival = str(PAR_index) + ',' + ','.join(str(i) for i in survival)
        with open(filename_stopt, 'a', newline='') as f:
            writer = csv.writer(f)
            writer.writerow(survival.split(','))

    if filename_dosage:
        with open(filename_dosage, 'a', newline='') as f:
            writer = csv.writer(f)
            for i_dose in [strategy_result[i_strategy]['dose'] for i_strategy in strategy_default]:
                writer.writerow(i_dose.split(','))

    # if filename_pop:
    #     with open(filename_pop, 'a', newline='') as f:
    #         writer = csv.writer(f)
    #         for i_pop in [strategy_result[i_strategy]['pop'] for i_strategy in strategy_default]:
    #             writer.writerow(i_pop)

    return


def DPM_reversiblemodel_readsave_assign(strategy, strategy_text, cured_survival, non_report, PAR_index,
                                        filename_stopt, filename_dosage, filename_pop):
    if strategy is not None:
        N, d = (strategy['N'], strategy['d']) if strategy_text in strategy_default[3:] else (strategy, None)
        if strategy_text in strategy_default[4:]:
            if filename_dosage:
                d = DPM_reversiblemodel_miscellaneous_findtreat(d, DPM_reversiblemodel_stime(N)) if isinstance(N, np.ndarray) else \
                    DPM_reversiblemodel_miscellaneous_findtreat(d, d.shape[1])
                d = [i.replace('σ', 'sigma') for i in d]
                d = d + ['-1'] * (Numstep - len(d))
                d = d[:Numstep]
                assert len(d) == Numstep
                d = ','.join([i for i in d])
                dose_strategy = ','.join([PAR_index, strategy_text, d])
            else:
                dose_strategy = ''
        elif strategy_text in strategy_default[:3]:
            dose_strategy = ','.join([PAR_index, strategy_text]) if filename_dosage else ''
        elif strategy_text == strategy_default[3]:
            d = [i.replace('σ', 'sigma') for i in d]
            dose_strategy = ','.join([PAR_index, strategy_text, *d]) if filename_dosage else ''
        else:
            raise ValueError('type of d is not correct.')

        if filename_stopt:
            if isinstance(N, np.ndarray):
                if all(N[:, -1] < 1):
                    survival_Strategy = cured_survival
                else:
                    survival_Strategy = DPM_reversiblemodel_stime(N)
                    if survival_Strategy == duration_5year and sum(N[:, -1]) < Limit_mortality:
                        survival_Strategy += 1
            elif isinstance(N, int):
                survival_Strategy = N
            else:
                raise ValueError('type of N is not correct.')
        else:
            survival_Strategy = ''

        if filename_pop:
            # N = np.around(N[:, ::Stepsize])
            N = np.around(N[:, [0, -1]])
            N = N.astype(int).T.tolist()
            N = [tuple(i_N) for i_N in N]
            N = [str(i_N) for i_N in N] + ['-1'] * (Numstep - len(N))
            pop_strategy = [PAR_index, strategy_text]
            pop_strategy.extend(N[:Numstep])
        else:
            pop_strategy = ''

    else:
        survival_Strategy = non_report
        dose_strategy = PAR_index + ',' + strategy_text + ',' + non_report
        pop_strategy = [PAR_index, strategy_text, str(non_report)]

    return {'stopt': survival_Strategy, 'dose': dose_strategy, 'pop': pop_strategy}


def DPM_reversiblemodel_readsave_dr(filename):
    df = pd.read_csv(filename)
    return df


def DPM_reversiblemodel_readsave_stoptime_csv(filename, Strategy_name):
    stoptime_key = ['paramID']
    stoptime_key.extend(Strategy_name)
    df_stoptime = pd.read_csv(filename, usecols=stoptime_key, dtype='float')
    assert not df_stoptime.isnull().values.any()
    data = df_stoptime.to_dict('list')
    return data


def DPM_reversiblemodel_readsave_dosage_csv(filename, strategy_name):
    df_dosage = pd.read_csv(filename, low_memory=False)
    strategy_name = [i for i in strategy_name if i not in strategy_default[:3]]
    df_dosage = df_dosage.loc[df_dosage['Strategy name'].isin(strategy_name)]
    # assert not df_dosage.isnull().values.any()

    paramID = []
    dosage_out, dosage_out_first2step = dict(), dict()
    for i_strategy in strategy_name:
        i_dosage = df_dosage.loc[df_dosage['Strategy name'].str.lower() == i_strategy.lower()]
        paramID.append(tuple(i_dosage['paramID'].values))
        i_dosage = i_dosage.drop(columns=['paramID', 'Strategy name'])
        i_dosage = i_dosage.fillna('-1')
        i_dosage.replace(-1, '-1', inplace=True)
        i_dosage_first2step = i_dosage.iloc[:, :2]
        dosage_out[i_strategy] = list(i_dosage.apply(';'.join, axis=1).astype(str).values)
        dosage_out_first2step[i_strategy] = list(i_dosage_first2step.apply(';'.join, axis=1).astype(str).values)
    # all paramID are the same
    assert all(ele == paramID[0] for ele in paramID)
    dosage_out['paramID'] = paramID[0]

    if 'strategy0' in strategy_name:
        bool_diff_dosage = [x != y for x, y in zip(dosage_out_first2step['strategy0'], dosage_out_first2step[strategy_name[-1]])]
        bool_same_dosage = [x == y for x, y in zip(dosage_out_first2step['strategy0'], dosage_out_first2step[strategy_name[-1]])]
        assert sum(bool_diff_dosage) + sum(bool_same_dosage) == len(dosage_out['paramID'])
    else:
        bool_diff_dosage = bool_same_dosage = None

    dosage_idx = dict(zip(['paramID', 'bool_diff_dosage', 'bool_same_dosage'], [dosage_out['paramID'], bool_diff_dosage, bool_same_dosage]))
    return dosage_out, dosage_idx


def DPM_reversiblemodel_readsave_pop_csv(filename, strategy_name):
    df_pop = pd.read_csv(filename, low_memory=False)
    df_pop = df_pop.loc[df_pop['Strategy name'].isin(strategy_name)]
    # assert not df_dosage.isnull().values.any()

    paramID, pop_out = [], dict()
    for i_strategy in strategy_name:
        i_pop = df_pop.loc[df_pop['Strategy name'].str.lower() == i_strategy.lower()]
        paramID.append(tuple(i_pop['paramID'].astype(int).values))
        i_pop = i_pop.drop(columns=['paramID', 'Strategy name'])
        i_pop = i_pop.fillna('-1')
        i_pop.replace(-1, '-1', inplace=True)
        pop_out[i_strategy] = list(i_pop.apply(';'.join, axis=1).astype(str).values)
    # all paramID are the same
    assert all(ele == paramID[0] for ele in paramID)
    pop_out['paramID'] = paramID[0]
    return pop_out
