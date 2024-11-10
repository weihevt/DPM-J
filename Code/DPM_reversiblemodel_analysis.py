from DPM_reversiblemodel_constant import *
from DPM_reversiblemodel_plot import \
    DPM_reversiblemodel_plot_hist, \
    DPM_reversiblemodel_plot_linearreg, \
    DPM_reversiblemodel_plot_stackplot
from DPM_reversiblemodel_miscellaneous import \
    DPM_reversiblemodel_miscellaneous_indextrueinlist


def DPM_reversiblemodel_analysis_KM(data, duration):
    E = [1 if i_val <= duration else 0 for i_val in data]
    epf = ExponentialFitter().fit(data, E)
    hazard = {'mean': epf.hazard_.mean(), 'ci': epf.confidence_interval_hazard_.mean()}

    kmf = KaplanMeierFitter()
    kmf.fit(data, E)
    km_interval = kmf.confidence_interval_survival_function_
    km = kmf.survival_function_

    t = km.index.values
    val = km.iloc[:].values.flatten()
    interval_lower = km_interval.iloc[:, 0].values.flatten()
    interval_upper = km_interval.iloc[:, 1].values.flatten()
    median_survival = kmf.median_survival_time_

    idx = np.where(t <= duration)[0]
    t, val, interval_lower, interval_upper = t[idx], val[idx], interval_lower[idx], interval_upper[idx]
    t, val, interval_lower, interval_upper = \
        np.append(t, duration), np.append(val, val[-1]), np.append(interval_lower, interval_lower[-1]), \
        np.append(interval_upper, interval_upper[-1])

    return {'t': t, 'val': val, 'median_survival': median_survival, 'interval_lower': interval_lower,
            'interval_upper': interval_upper, 'hazard': hazard}


def DPM_reversiblemodel_analysis_hz(stoptime, strategy_name, Simduration):
    p = DPM_analysis_pairwise_logrank_test(stoptime, strategy_name, Simduration)
    result = dict(zip(strategy_name, [None]*len(strategy_name)))
    hz = dict()
    for i_strategy in strategy_name:
        hz[i_strategy] = deepcopy(result)
    with tqdm(total=len(list(itertools.product(strategy_name, strategy_name))), ncols=150,
              desc='Calculating hz and sigbetter') as pbar:
        for i in list(itertools.product(strategy_name, strategy_name)):
            if i[0] != i[1]:
                test, ref, name = stoptime[i[0]], stoptime[i[1]], ', '.join(i)
                hz[i[0]][i[1]] = DPM_reversiblemodel_analysis_HZ(ref, test, Simduration)
            pbar.update(1)
    return dict(zip(['hz', 'p'], [hz, p]))


def DPM_reversiblemodel_analysis_sigbetter(stoptime_ref, stoptime_test, paramID):
    ind = np.logical_and(np.array(stoptime_test) > np.array(stoptime_ref) + 30*2,
                         np.array(stoptime_test) > 1.25 * np.array(stoptime_ref))
    return list(itertools.compress(paramID, ind))


def DPM_reversiblemodel_analysis_dose(dose, strategyname, info_stoptime, stoptime, pathsave):
    def DPM_reversiblemodel_analysis_dose_1(data_, ind_, bins_, xticks_, title_, label_, xlim_):
        fig, ax = plt.subplots()
        DPM_reversiblemodel_plot_hist(fig, ax, data_, bins=bins_, xticks=xticks_, xlabel=xlabel, title=title_, facecolor=color_total,
                                      label='Total')
        DPM_reversiblemodel_plot_hist(fig, ax, data_[ind_], bins=bins_, xticks=xticks_, xlabel=xlabel, title=title_, facecolor=color_ind,
                                      label=label_)
        plt.xlim(xlim_)
        plt.legend(loc='upper right')
        plt.rcParams.update({'font.size': 13})
        plt.tight_layout()
        return

    result = {i_strategy: [] for i_strategy in strategyname if i_strategy not in ('mono drug1', 'mono drug2', 'mono half')}
    num_change, norm_num_change, num_cycleinDPM, norm_num_cycleinDPM, σ1dose, norm_σ1dose, σ2dose, norm_σ2dose, \
        bestcycle_period, bestcycle_firstdrug = (deepcopy(result), deepcopy(result), deepcopy(result), deepcopy(result), deepcopy(result),
                                                 deepcopy(result), deepcopy(result), deepcopy(result), [], [])
    for i_strategy in strategyname:
        if i_strategy in ('mono drug1', 'mono drug2', 'mono half'):
            continue
        else:
            i_dose_strategy = dose[i_strategy]
        for i, i_dose in enumerate(i_dose_strategy):
            i_dose = i_dose.split(';')
            i_num_change, i_num_cycleinDPM, i_σ1dose, i_σ2dose = 0, 0, 0, 0
            if '-1' in i_dose:
                i_dose = i_dose[:i_dose.index('-1')]
            if i_strategy == 'best cycle':
                i_bestcycle_period = [int(j_dose.split(' ')[0][:-1]) for j_dose in i_dose]
                i_bestcycle_firstσ = [j_dose.split(' ')[1] for j_dose in i_dose]
                i_bestcycle_firstσ = [j_dose.replace('sigma', 'σ') for j_dose in i_bestcycle_firstσ]
                bestcycle_period.append(list(set(i_bestcycle_period)))
                bestcycle_firstdrug.append(list(set(i_bestcycle_firstσ)))
            else:
                i_current = i_dose[0]
                if i_strategy == 'strategy0':
                    assert set(i_dose).issubset({'sigma1 full', 'sigma2 full'})
                if i_strategy in ['strategy1', 'strategy2', 'strategy3', 'strategy4']:
                    assert set(i_dose).issubset({'sigma1 full', 'sigma2 full', 'half'})
                for j, j_step in enumerate(i_dose):
                    if j_step != i_current:
                        i_num_change += 1
                        i_current = j_step
                    if j_step not in ['sigma1 full', 'sigma2 full', 'half']:
                        i_num_cycleinDPM += 1
                    if j_step == 'sigma1 full':
                        i_σ1dose += 1
                    elif j_step == 'sigma2 full':
                        i_σ2dose += 1
                    else:
                        i_σ1dose += 0.5
                        i_σ2dose += 0.5
                num_change[i_strategy].append(i_num_change)
                norm_num_change[i_strategy].append(i_num_change/len(i_dose))
                num_cycleinDPM[i_strategy].append(i_num_cycleinDPM)
                norm_num_cycleinDPM[i_strategy].append(i_num_cycleinDPM/len(i_dose))
                σ1dose[i_strategy].append(i_σ1dose)
                σ2dose[i_strategy].append(i_σ2dose)
                norm_σ1dose[i_strategy].append(i_σ1dose/len(i_dose))
                norm_σ2dose[i_strategy].append(i_σ2dose/len(i_dose))

                assert (len(num_change[i_strategy]) == len(num_cycleinDPM[i_strategy]) == len(σ1dose[i_strategy])
                        == len(σ2dose[i_strategy]) == i+1)
            if i_strategy != 'best cycle':
                assert i_σ1dose/len(i_dose) + i_σ2dose/len(i_dose) == 1
            if i_strategy == 'strategy0':
                assert i_num_change <= 1
        if i_strategy != 'best cycle':
            assert (len(num_change[i_strategy]) == len(num_cycleinDPM[i_strategy]) == len(σ1dose[i_strategy])
                    == len(σ2dose[i_strategy]) == len(i_dose_strategy))

    color_total = 'b'
    color_ind = 'r'
    # color_numerical = 'g'
    strategyname_plot = ['strategy1', 'strategy1 cycle', 'strategy2', 'strategy2 cycle']
    data = [num_change, norm_num_change, σ1dose, σ2dose, norm_σ1dose, norm_σ2dose]
    data_cycle = [num_cycleinDPM, norm_num_cycleinDPM]
    bins = [np.arange(0, Numstep + 1), np.arange(0, 10, 1)/10, np.arange(0, Numstep + 1), np.arange(0, Numstep + 1),
            np.arange(0, 10, 1)/10,  np.arange(0, 10, 1)/10]
    bins_cycle = [np.arange(0, Numstep + 1), np.arange(0, 10, 1)/10]
    xlim = [(0, Numstep+1), (0, 1), (0, Numstep+1), (0, Numstep+1), (0, 1), (0, 1)]
    xlim_cycle = [(0, Numstep+1), (0, 1)]
    title = ['numdrugchange', 'normdrugchange', 'drug1dose', 'drug2dose', 'normdrug1dose', 'normdrug2dose']
    title_cycle = ['numcycle', 'normnumcycle']
    label = ['numerical best', 'sig best']
    ind = ['ind_beststrategy', 'ind_sigbeststrategy']
    xlabel = 'values'
    for i_strategy in strategyname_plot:
        for i, i_data in enumerate(data):
            i_data = np.array(i_data[i_strategy])
            i_bins, i_xticks, i_xlim = bins[i], bins[i], xlim[i]
            i_title = f'{title[i]} in {i_strategy}'
            for j, j_ind in enumerate(ind):
                j_label = label[j]
                j_ind = info_stoptime[j_ind][i_strategy]
                DPM_reversiblemodel_analysis_dose_1(i_data, j_ind, i_bins, i_xticks, i_title, j_label, i_xlim)
                plt.xlim(xlim[i])
                filename = os.path.join(pathsave, f'{title[i]} {i_strategy} {j_label}.pdf')
                plt.savefig(filename, format='pdf', bbox_inches='tight')
                plt.close()

        # num_cycleinDPM and norm_num_cycleinDPM
        if i_strategy in ('strategy1 cycle', 'strategy2 cycle'):
            for i, i_data in enumerate(data_cycle):
                i_data = np.array(i_data[i_strategy])
                i_bins, i_xticks, i_xlim = bins_cycle[i], bins_cycle[i], xlim_cycle[i]
                i_title = f'{title_cycle[i]} in {i_strategy}'
                for j, j_ind in enumerate(ind):
                    j_label = label[j]
                    j_ind = info_stoptime[j_ind][i_strategy]
                    DPM_reversiblemodel_analysis_dose_1(i_data, j_ind, i_bins, i_xticks, i_title, j_label, i_xlim)
                    plt.xlim(xlim[i])
                    filename = os.path.join(pathsave, f'{title_cycle[i]} {i_strategy} {j_label}.pdf')
                    plt.savefig(filename, format='pdf', bbox_inches='tight')
                    plt.close()
    data = [norm_num_change, norm_σ1dose, norm_σ2dose]
    data_cycle = [norm_num_cycleinDPM]
    title = ['normdrugchange', 'normdrug1dose', 'normdrug2dose']
    title_cycle = ['normnumcycle']
    ylim = (0, duration_5year + 30)
    for i_strategy in strategyname_plot:
        i_stoptime = stoptime[i_strategy]
        i_stoptime = np.array([i_val if i_val < duration_5year else duration_5year for i_val in i_stoptime])
        for i, i_data in enumerate(data):
            i_data = np.array(i_data[i_strategy])
            i_title = f'{title[i]} in {i_strategy}'
            assert len(i_stoptime) == len(i_data)
            i_model = sm.OLS(i_stoptime, i_data).fit()
            DPM_reversiblemodel_plot_linearreg(i_data, i_stoptime, i_model, i_title, ylim)
            plt.close()
            for j, j_ind in enumerate(ind):
                j_label = label[j]
                j_ind = info_stoptime[j_ind][i_strategy]
                j_data = i_data[j_ind]
                j_stoptime = i_stoptime[j_ind]
                i_title = f'{j_label} {title[i]} in {i_strategy}'
                DPM_reversiblemodel_plot_linearreg(j_data, j_stoptime, i_model, i_title, ylim)
                plt.close()
        if i_strategy in ('strategy1 cycle', 'strategy2 cycle'):
            for i, i_data in enumerate(data_cycle):
                i_data = np.array(i_data[i_strategy])
                i_title = f'{title_cycle[i]} in {i_strategy}'
                assert len(i_stoptime) == len(i_data)
                i_model = sm.OLS(i_stoptime, i_data).fit()
                DPM_reversiblemodel_plot_linearreg(i_data, i_stoptime, i_model, i_title, ylim)
                plt.close()
                for j, j_ind in enumerate(ind):
                    j_label = label[j]
                    j_ind = info_stoptime[j_ind][i_strategy]
                    j_data = i_data[j_ind]
                    j_stoptime = i_stoptime[j_ind]
                    i_title = f'{j_label} {title[i]} in {i_strategy}'
                    DPM_reversiblemodel_plot_linearreg(j_data, j_stoptime, i_model, i_title, ylim)
                    plt.close()

    return dict(zip(['num change', 'norm num change', 'num cycleinDPM', 'norm num cycleinDPM', 'σ1dose', 'norm σ1dose',
                     'σ2dose', 'norm σ2dose', 'bestcycle_period', 'bestcycle_firstdrug'],
                    [num_change, norm_num_change, num_cycleinDPM, norm_num_cycleinDPM, σ1dose, norm_σ1dose, σ2dose,
                     norm_σ2dose, bestcycle_period, bestcycle_firstdrug]))


def DPM_reversiblemodel_analysis_HZ(data_ref, data, Simduration):
    if (len(data_ref) != 0) & (len(data) != 0):
        treat = np.concatenate((np.zeros(len(data_ref)), np.ones(len(data))))
        val = data_ref + data
        E = [1 if i_val <= Simduration else 0 for i_val in data_ref]
        E.extend([1 if i_val <= Simduration else 0 for i_val in data])
        E = np.array(E)
        d = {'val': val, 'E': E, 'treat': treat}
        df = pd.DataFrame(data=d)
        cph = CoxPHFitter()
        cph.fit(df, duration_col='val', event_col='E')
        hz_ratio = cph.hazard_ratios_.values[0]
    else:
        hz_ratio = None
    return hz_ratio


def DPM_analysis_pairwise_logrank_test(stoptime, Strategy_name, Simduration):
    G, T_, E = tuple(), tuple(), tuple()
    flag_empty = False
    for i_strategy in Strategy_name:
        i_stop = stoptime[i_strategy]
        if not i_stop:
            flag_empty = True
            break
        G = G + tuple(len(i_stop) * [str(i_strategy)])
        E = E + tuple([1 if i_val <= Simduration else 0 for i_val in i_stop])
        T_ = T_ + tuple(i_stop)
    p = statistics.pairwise_logrank_test(T_, G, E) if not flag_empty else None
    p_out = dict(zip(p.name, p.p_value)) if not flag_empty else None
    return p_out


def DPM_reversiblemodel_analysis_stoptime(stoptime, para, dosage, strategy_name, simduration, mutonly=False, reveronly=False):
    def DPM_reversiblemodel_analysis_stoptime_1(stoptime_df_, para_df_):
        survival_median = stoptime_df_.median()[strategy_name].astype(int)
        survival_mean = stoptime_df_.mean()[strategy_name].apply(np.ceil).astype(int)
        num_survival_5y = stoptime_df_[strategy_name].apply(lambda x: x >= simduration).apply(sum).to_dict()
        percentage_survival_5y = (stoptime_df_[strategy_name].apply(lambda x: x >= simduration).apply(sum)/stoptime_df_.shape[0] * 100)
        percentage_survival_5y = percentage_survival_5y.round(2).to_dict()

        result = dict(zip(strategy_name, [None] * len(strategy_name)))

        (ind_beststrategy, improve_beststrategy, num_beststrategy, percentage_beststrategy,
         ind_sigbeststrategy, improve_sigbeststrategy, num_sigbeststrategy, percentage_sigbeststrategy,
         ind_betterstrategy, improve_betterstrategy, num_betterstrategy, percentage_betterstrategy,
         ind_sigbetterstrategy, improve_sigbetterstrategy, num_sigbetterstrategy, percentage_sigbetterstrategy
         ) = (deepcopy(result), deepcopy(result), deepcopy(result), deepcopy(result),
              deepcopy(result), deepcopy(result), deepcopy(result), deepcopy(result),
              deepcopy(result), deepcopy(result), deepcopy(result), deepcopy(result),
              deepcopy(result), deepcopy(result), deepcopy(result), deepcopy(result))

        for i_strategy in strategy_name:
            sub_strategy_name = [j_strategy for j_strategy in strategy_name if j_strategy is not i_strategy]
            max_sub_strategy = stoptime_df_[sub_strategy_name].max(axis=1)
            i_diff = stoptime_df_[i_strategy] - max_sub_strategy
            # No. of cases strategy numerically better than all the others
            i_ind = pd.Series(stoptime_df_[i_strategy] > max_sub_strategy)
            ind_beststrategy, improve_beststrategy, num_beststrategy, percentage_beststrategy \
                = DPM_reversiblemodel_analysis_stoptime_2(i_ind, i_diff, i_strategy, ind_beststrategy, improve_beststrategy, num_beststrategy,
                                                          percentage_beststrategy)

            # No. of cases strategy significantly better than all the others.
            i_ind = (stoptime_df_[i_strategy] > max_sub_strategy + 8 * 7) & (stoptime_df_[i_strategy] > 1.25 * max_sub_strategy)
            ind_sigbeststrategy, improve_sigbeststrategy, num_sigbeststrategy, percentage_sigbeststrategy \
                = DPM_reversiblemodel_analysis_stoptime_2(i_ind, i_diff, i_strategy, ind_sigbeststrategy, improve_sigbeststrategy,
                                                          num_sigbeststrategy, percentage_sigbeststrategy)

            i_ind_betterstrategy, i_improve_betterstrategy, i_num_betterstrategy, i_percentage_betterstrategy \
                = deepcopy(result), deepcopy(result), deepcopy(result), deepcopy(result)

            i_ind_betterstrategy[i_strategy], i_improve_betterstrategy[i_strategy], i_num_betterstrategy[i_strategy], \
                i_percentage_betterstrategy[i_strategy] = 'N.A.', 'N.A.', 'N.A.', 'N.A.'

            i_ind_sigbetterstrategy, i_improve_sigbetterstrategy, i_num_sigbetterstrategy, i_percentage_sigbetterstrategy \
                = deepcopy(result), deepcopy(result), deepcopy(result), deepcopy(result)

            i_ind_sigbetterstrategy[i_strategy], i_improve_sigbetterstrategy[i_strategy], i_num_sigbetterstrategy[i_strategy], \
                i_percentage_sigbetterstrategy[i_strategy] = 'N.A.', 'N.A.', 'N.A.', 'N.A.'

            for i_substrategy in sub_strategy_name:
                # No. of cases strategy numerically better than all the others
                i_ind = pd.Series(stoptime_df_[i_strategy] > stoptime_df_[i_substrategy])
                i_diff = stoptime_df_[i_strategy] - stoptime_df_[i_substrategy]

                i_ind_betterstrategy, i_improve_betterstrategy, i_num_betterstrategy, i_percentage_betterstrategy \
                    = DPM_reversiblemodel_analysis_stoptime_2(i_ind, i_diff, i_substrategy, i_ind_betterstrategy, i_improve_betterstrategy,
                                                              i_num_betterstrategy, i_percentage_betterstrategy)

                # No. of cases strategy significantly better.
                i_ind = (stoptime_df_[i_strategy] > stoptime_df_[i_substrategy] + 8 * 7) & \
                        (stoptime_df_[i_strategy] > 1.25 * stoptime_df_[i_substrategy])
                i_ind_sigbetterstrategy, i_improve_sigbetterstrategy, i_num_sigbetterstrategy, i_percentage_sigbetterstrategy \
                    = DPM_reversiblemodel_analysis_stoptime_2(i_ind, i_diff, i_substrategy, i_ind_sigbetterstrategy, i_improve_sigbetterstrategy,
                                                              i_num_sigbetterstrategy, i_percentage_sigbetterstrategy)

            ind_betterstrategy[i_strategy], improve_betterstrategy[i_strategy], num_betterstrategy[i_strategy], \
                percentage_betterstrategy[i_strategy] = i_ind_betterstrategy, i_improve_betterstrategy, i_num_betterstrategy, \
                i_percentage_betterstrategy

            ind_sigbetterstrategy[i_strategy], improve_sigbetterstrategy[i_strategy], num_sigbetterstrategy[i_strategy], \
                percentage_sigbetterstrategy[i_strategy] = i_ind_sigbetterstrategy, i_improve_sigbetterstrategy, i_num_sigbetterstrategy, \
                i_percentage_sigbetterstrategy

        info_stoptime_keys = ['stoptime_df', 'para_df', 'survival_median', 'survival_mean', 'num_survival_5y', 'percentage_survival_5y',
                              'ind_beststrategy', 'improve_beststrategy', 'num_beststrategy', 'percentage_beststrategy',
                              'ind_sigbeststrategy', 'improve_sigbeststrategy', 'num_sigbeststrategy',
                              'percentage_sigbeststrategy', 'ind_betterstrategy', 'improve_betterstrategy',
                              'num_betterstrategy', 'percentage_betterstrategy', 'ind_sigbetterstrategy',
                              'improve_sigbetterstrategy', 'num_sigbetterstrategy', 'percentage_sigbetterstrategy']
        info_stoptime_ = [stoptime_df_, para_df_, survival_median, survival_mean, num_survival_5y, percentage_survival_5y,
                          ind_beststrategy, improve_beststrategy, num_beststrategy, percentage_beststrategy,
                          ind_sigbeststrategy, improve_sigbeststrategy, num_sigbeststrategy, percentage_sigbeststrategy,
                          ind_betterstrategy, improve_betterstrategy, num_betterstrategy, percentage_betterstrategy,
                          ind_sigbetterstrategy, improve_sigbetterstrategy, num_sigbetterstrategy,
                          percentage_sigbetterstrategy]
        return dict(zip(info_stoptime_keys, info_stoptime_))

    def DPM_reversiblemodel_analysis_stoptime_2(ind, val, strategy_, ind_, improve_, num_, percentage_):
        ind_true = [i for i, x in enumerate(ind.to_list()) if x]
        val_sel = val.iloc[ind_true].values
        ind_[strategy_] = ind_true
        improve_[strategy_] = val_sel
        num_[strategy_] = len(ind_true)
        percentage_[strategy_] = len(ind_true)/len(ind)*100
        percentage_[strategy_] = round(percentage_[strategy_], 2)
        return ind_, improve_, num_, percentage_

    def DPM_reversiblemodel_analysis_stoptime_3(namebetter_, nameworse_):
        idx = info['ind_sigbetterstrategy'][namebetter_][nameworse_]
        idx_max = np.argmax(np.array(stoptime[namebetter_])[idx] - np.array(stoptime[nameworse_])[idx])
        idx = idx[idx_max]
        para_sel = {key: val[idx] for (key, val) in para.items()}
        return para_sel

    def DPM_reversiblemodel_analysis_stoptime_4(namebetter_, nameworse_):
        # select best name_better_ compared to all name_worse
        max_nameworse = stoptime_df[nameworse_].max(axis=1)
        diff = stoptime_df[namebetter_] - max_nameworse
        idx = pd.Index(diff).sort_values(ascending=False, return_indexer=True)[1]
        idx_choose = 0
        if not mutonly and not reveronly:
            for i in idx:
                i_dosage_ = dosage[namebetter_][i]
                assert dosage['paramID'][i] == para['paramID'][i]
                i_dosage_ = i_dosage_.split(';')
                if i_dosage_.count('3w sigma2 1th') >= 2 or i_dosage_.count('3w sigma2 1th') >= 2:
                    idx_choose = i
                    break
            para_sel = para_df.iloc[idx_choose].to_dict()
        else:
            para_sel = para_df.iloc[idx[0]].to_dict()
        return para_sel

    stoptime_df = pd.DataFrame.from_dict(stoptime)
    stoptime_df = stoptime_df.set_index('paramID', drop=False)
    para_df = pd.DataFrame.from_dict(para)
    para_df = para_df.set_index('paramID', drop=False)

    info = DPM_reversiblemodel_analysis_stoptime_1(stoptime_df, para_df)

    strategybetter, strategyworse = 'best cycle', 'mono half'
    _ = DPM_reversiblemodel_analysis_stoptime_3(strategybetter, strategyworse)

    if reveronly:
        strategysuperior = 'best cycle'
        strategy_name_compare = ['mono half', 'strategy0']
        strategy_name_run = ['mono 1', 'mono 2', 'best cycle', 'mono half', 'strategy0']
    elif mutonly:
        strategysuperior = 'strategy2'
        strategy_name_compare = ['strategy0', 'strategy1']
        strategy_name_run = ['strategy2', 'strategy0', 'strategy1']
    else:
        strategysuperior = 'strategy2 cycle'
        strategy_name_compare = ['best cycle', 'strategy1', 'strategy1 cycle', 'strategy2']
        strategy_name_run = ['best cycle', 'strategy1 cycle', 'strategy2', 'strategy2 cycle']
    para_ex = DPM_reversiblemodel_analysis_stoptime_4(strategysuperior, strategy_name_compare)
    return info, para_ex, strategy_name_run


def DPM_reversiblemodel_analysis_Nflow(Nflow, strategy_name, para, pathsave):
    def DPM_reversiblemodel_analysis_Nflow_1(i_Nflow, i_ylim, i_title):
        i_Nflow_pro, i_Nflow_tran, i_Nflow_muta = i_Nflow[cellflow_pro, :], i_Nflow[cellflow_tran, :], i_Nflow[cellflow_muta, :]

        # Nflow col
        titlestr = i_strategy + ' pro ' + i_title
        legstr = [cellflow[i] for i in cellflow_pro]
        DPM_reversiblemodel_plot_stackplot(i_Nflow_pro, titlestr, legstr, i_ylim)
        plt.savefig(os.path.join(i_pathsave, titlestr + '.pdf'), format='pdf', bbox_inches='tight')
        plt.close('all')

        titlestr = i_strategy + ' tran ' + i_title
        legstr = [cellflow[i] for i in cellflow_tran]
        DPM_reversiblemodel_plot_stackplot(i_Nflow_tran, titlestr, legstr, i_ylim)
        plt.savefig(os.path.join(i_pathsave, titlestr + '.pdf'), format='pdf', bbox_inches='tight')
        plt.close('all')

        titlestr = i_strategy + ' muta ' + i_title
        legstr = [cellflow[i] for i in cellflow_muta]
        DPM_reversiblemodel_plot_stackplot(i_Nflow_muta, titlestr, legstr, i_ylim)
        plt.savefig(os.path.join(i_pathsave, titlestr + '.pdf'), format='pdf', bbox_inches='tight')
        plt.close('all')
        return

    Nflowcol, Nflowcolextend, num = Nflow['valcol'], Nflow['valcolextend'], Nflow['num']
    for i_strategy in strategy_name:
        i_pathsave = os.path.join(pathsave, i_strategy)
        if not os.path.exists(i_pathsave):
            os.makedirs(i_pathsave)

        i_ind, = (num[i_strategy] != 0).nonzero()
        i_num = np.tile(num[i_strategy][i_ind].reshape(1, len(i_ind)), (len(cellflow), 1))
        i_Nflowcol = Nflowcol[i_strategy][:, i_ind]/i_num
        i_Nflowcolextend = Nflowcolextend[i_strategy]/len(para['paramID'])
        assert all([isclose(x, 1) for x in np.sum(i_Nflowcol, axis=0)])
        assert all([isclose(x, 1) for x in np.sum(i_Nflowcolextend, axis=0)])

        DPM_reversiblemodel_analysis_Nflow_1(i_Nflowcol, (0, 1), 'col')
        DPM_reversiblemodel_analysis_Nflow_1(i_Nflowcolextend, (0, 1), 'colextend')
    return


def DPM_reversiblemodel_analysis_sensitivity(stoptime, para, strategy_name, pathsave):
    def DPM_reversiblemodel_analysis_sensitivity_1(i_):
        i_par_ = para_df.iloc[i_]
        for i_key in par_sensitivity:
            idx = para_set[i_key].index(i_par_[i_key])
            if idx < len(para_set[i_key]) - 1:
                i_par_increase = deepcopy(i_par_)
                i_par_increase[i_key] = para_set[i_key][idx + 1]
                mask = pd.Series(True, index=para_df.index)
                for j_key in par_sensitivity:
                    j_flag = para_df[j_key] == i_par_increase[j_key]
                    mask = mask & j_flag
                    if mask.sum() == 0:
                        break
                if mask.sum() != 0:
                    assert mask.sum() == 1
                    i_idx = mask[mask].index.values[0]
                    i_stoptime, i_stoptime_sel = stoptime_df.iloc[i_idx, :], stoptime_df.iloc[i_, :]
                    for i_strategy_ in strategy_name:
                        i_sensitivity = (i_stoptime_sel.loc[i_strategy_] - i_stoptime[i_strategy_])/i_stoptime[i_strategy_]
                        sensitivity[i_key][i_strategy_].append(abs(i_sensitivity))

            i_sensitivity = sensitivity[i_key]
            i_num_strategy = [len(i_sensitivity[i_name]) for i_name in strategy_name]
            assert i_num_strategy.count(i_num_strategy[0]) == len(i_num_strategy)

        # pbar.update(1)
        info = ['{:.2e}'.format(i_/len(para_df))]
        num_ = [len(sensitivity[i_name][strategy_name[0]]) for i_name in sensitivity]
        info.extend(num_)
        print(info)
        return

    filename_sensitivity = os.path.join(pathsave, 'sensitivity.pckl')
    filename_sensitivity_para_df = os.path.join(pathsave, 'sensitivity_para_df.pckl')
    if not os.path.exists(filename_sensitivity):
        assert stoptime['paramID'] == para['paramID']
        para_set = dict(zip(par_sensitivity, [None]*len(par_sensitivity)))
        # sensitivity = dict(zip(strategy_name, [dict(zip(par_sensitivity, [[]]*len(par_sensitivity)))]*len(strategy_name)))
        result = dict(zip(strategy_name, [[-1] for _ in range(len(strategy_name))]))
        sensitivity = dict(zip(par_sensitivity, [deepcopy(result) for _ in range(len(par_sensitivity))]))
        stoptime_df = pd.DataFrame.from_dict(stoptime)
        if not os.path.exists(filename_sensitivity_para_df):
            para['S_01_2_S_01_1'] = [i/j for i, j in zip(para['S_01_2'], para['S_01_1'])]
            para['SR_00_1'] = [i/j for i, j in zip(para['S_00_1'], para['S_01_1'])]
            para['SR_00_2'] = [i / j for i, j in zip(para['S_00_2'], para['S_01_2'])]
            para['SR_11_1'] = [i / j for i, j in zip(para['S_11_1'], para['S_01_1'])]
            para['SR_11_2'] = [i / j for i, j in zip(para['S_11_2'], para['S_01_2'])]

            for i in par_sensitivity:
                i_par = list(set(para[i]))
                if i in ['SR_00_1', 'SR_00_2']:
                    i_par = ['{:.1e}'.format(x) for x in i_par]
                else:
                    i_par = ['{:.2e}'.format(x) for x in i_par]
                i_par = [float(x) for x in i_par]
                i_par = list(set(i_par))
                i_par.sort()
                para_set[i] = i_par

            para_df = pd.DataFrame.from_dict(para)
            para_df = para_df[['paramID'] + par_sensitivity]
            para_df.drop_duplicates(subset=para_df.columns.to_list()[1:], inplace=True)
            with tqdm(total=len(para_df), ncols=120) as pbar:
                for i in range(para_df.shape[0]):
                    assert para_df.loc[i, 'paramID'] == stoptime_df.loc[i, 'paramID']
                    for j in par_sensitivity:
                        ij_val = para_df.loc[i, j]
                        i_ind = [isclose(ij_val, i_val, rel_tol=1e-2) for i_val in para_set[j]]
                        i_ind = DPM_reversiblemodel_miscellaneous_indextrueinlist(i_ind)
                        assert len(i_ind) == 1
                        para_df.loc[i, j] = para_set[j][i_ind[0]]
                    pbar.update(1)
                with bz2.BZ2File(filename_sensitivity_para_df, 'wb') as f:
                    pickle.dump((para_df, para_set), f)
        else:
            with bz2.BZ2File(filename_sensitivity_para_df, 'rb') as f:
                para_df, para_set = pickle.load(f)
        num_cores = multiprocessing.cpu_count()
        # with tqdm(total=len(para_df), ncols=120) as pbar:
        Parallel(n_jobs=num_cores)(delayed(DPM_reversiblemodel_analysis_sensitivity_1)(i) for i in range(len(para_df)))
        # for i in range(len(para_df)):
        #     DPM_reversiblemodel_analysis_sensitivity_1(i)
        with bz2.BZ2File(filename_sensitivity, 'wb') as f:
            pickle.dump(sensitivity, f)
    else:
        with bz2.BZ2File(filename_sensitivity, 'rb') as f:
            sensitivity = pickle.load(f)

        result = dict(zip(strategy_name, [[-1] for _ in range(len(strategy_name))]))
        sensitivity_1 = dict(zip(par_sensitivity, [deepcopy(result) for _ in range(len(par_sensitivity))]))

        for i_parname in par_sensitivity:
            for i_strategy in strategy_name:
                sensitivity_1[i_parname][i_strategy] = sensitivity[i_strategy][i_parname]
    return


def DPM_reversiblemodel_analysis_sensitivity_αθμ(stoptime, para, filename_sensitivity_set, filename_sensitivity_para_df):
    def DPM_reversiblemodel_analysis_sensitivity_αθμ_1(i_, i_par_):
        para_np_cp = para_np[i_:, par_sensitivity_index_same]
        i_matching_rows, = np.equal(para_np_cp, i_par_[par_sensitivity_index_same]).all(axis=1).nonzero()
        if flag_assert:
            assert len(i_matching_rows) >= 1
        if len(i_matching_rows) > 1:
            i_matching_rows = i_matching_rows + i_
            i_str = str(tuple(np.sort(i_matching_rows))).replace(' ', '').replace('(', '').replace(')', '')
            # sensitivity_set[parname][i_] = i_str
            sensitivity_set.append(i_str)
        else:
            # sensitivity_set[parname][i_] = 'NaN'
            pass

        pbar.update(1)
        # if i_ % 1000 == 0:
        #     info = '{:.2e}'.format(i_/num_par)
        #     # num_ = [sum(1 for n in sensitivity_set[i_name] if n != str_format and n != 'NaN') for i_name in par_sensitivity2]
        #     # info.extend(num_)
        # print(i_)
        return

    if not os.path.exists(filename_sensitivity_set):
        assert stoptime['paramID'] == para['paramID']
        para_set = dict(zip(par_sensitivity, [None]*len(par_sensitivity)))
        stoptime_df = pd.DataFrame.from_dict(stoptime)
        if not os.path.exists(filename_sensitivity_para_df):
            para['S_01_2_S_01_1'] = [i/j for i, j in zip(para['S_01_2'], para['S_01_1'])]
            para['SR_00_1'] = [i/j for i, j in zip(para['S_00_1'], para['S_01_1'])]
            para['SR_00_2'] = [i/j for i, j in zip(para['S_00_2'], para['S_01_2'])]
            para['SR_11_1'] = [i/j for i, j in zip(para['S_11_1'], para['S_01_1'])]
            para['SR_11_2'] = [i/j for i, j in zip(para['S_11_2'], para['S_01_2'])]

            for i in par_sensitivity:
                i_par = list(set(para[i]))
                if i in ['SR_00_1', 'SR_00_2']:
                    i_par = ['{:.1e}'.format(x) for x in i_par]
                else:
                    i_par = ['{:.2e}'.format(x) for x in i_par]
                i_par = [float(x) for x in i_par]
                i_par = list(set(i_par))
                i_par.sort()
                para_set[i] = i_par

            para_set_num = {key: len(value) for key, value in para_set.items()}

            para_df = pd.DataFrame.from_dict(para)
            para_df = para_df[['paramID'] + par_sensitivity]
            para_df.drop_duplicates(subset=para_df.columns.to_list()[1:], inplace=True)
            with tqdm(total=len(para_df), ncols=120) as pbar:
                for i in range(para_df.shape[0]):
                    assert para_df.loc[i, 'paramID'] == stoptime_df.loc[i, 'paramID']
                    for j in par_sensitivity:
                        ij_val = para_df.loc[i, j]
                        i_ind = [isclose(ij_val, i_val, rel_tol=1e-2) for i_val in para_set[j]]
                        i_ind = DPM_reversiblemodel_miscellaneous_indextrueinlist(i_ind)
                        assert len(i_ind) == 1
                        para_df.loc[i, j] = para_set[j][i_ind[0]]
                    pbar.update(1)
                with bz2.BZ2File(filename_sensitivity_para_df, 'wb') as f:
                    pickle.dump((para_df, para_set, para_set_num), f)
        else:
            with bz2.BZ2File(filename_sensitivity_para_df, 'rb') as f:
                para_df, para_set, para_set_num = pickle.load(f)

        parname = 'αθμ'
        flag_assert = True
        num_cores = multiprocessing.cpu_count()
        para_np = para_df.iloc[:, 1:].to_numpy(dtype=int)
        num_par = len(para_df)
        # sensitivity_set = {parname: shared_memory.ShareableList([str_format] * num_par)}
        # Parallel(n_jobs=num_cores)(delayed(DPM_reversiblemodel_analysis_sensitivity2_1)(i, para_np[i, :])
        #                            for i in range(para_np.shape[0]))

        with tqdm(total=para_np.shape[0], ncols=120) as pbar:
            sensitivity_set = list()
            for i in range(para_np.shape[0]):
                DPM_reversiblemodel_analysis_sensitivity_αθμ_1(i, para_np[i, :])

        with bz2.BZ2File(filename_sensitivity_set, 'wb') as f:
            pickle.dump(sensitivity_set, f)
    return


def DPM_reversiblemodel_analysis_sensitivity_nonαθμ(stoptime, para, filename_sensitivity_set, filename_sensitivity_para_df):
    def DPM_reversiblemodel_analysis_sensitivity_nonαθμ_1(i_):
        i_par_ = para_np[i_, :]
        para_np_cp = para_np[i_:, :]
        for i_name in par_sensitivity_use2:
            i = ind_par_sensitivity_use[i_name]
            if flag_assert:
                assert i_par_[i] < para_set_num[i_name]
            i_pair_increase, i_pair_decrease = None, None
            if i_par_[i] < para_set_num[i_name] - 1:
                i_par_increase = deepcopy(i_par_)
                i_par_increase[i] = i_par_increase[i] + 1
                i_matching_rows, = np.equal(para_np_cp, i_par_increase).all(axis=1).nonzero()
                if i_matching_rows.any():
                    if flag_assert:
                        assert len(i_matching_rows) == 1
                    i_matching_rows = i_matching_rows[0] + i_
                    i_val = tuple(np.sort(np.array([i_, i_matching_rows])))
                    i_pair_increase = f'{i_val[0]},{i_val[1]}'
                    if flag_assert:
                        assert np.equal(np.delete(para_np[i_, :], i), np.delete(para_np[i_matching_rows, :], i)).all()
                        assert para_np[i_, i] != para_np[i_matching_rows, i]

            if i_par_[i] > 0:
                i_par_decrease = deepcopy(i_par_)
                i_par_decrease[i] = i_par_decrease[i] - 1
                i_matching_rows, = np.equal(para_np_cp, i_par_decrease).all(axis=1).nonzero()
                if i_matching_rows.any():
                    if flag_assert:
                        assert len(i_matching_rows) == 1
                    i_matching_rows = i_matching_rows[0] + i_
                    i_val = tuple(np.sort(np.array([i_, i_matching_rows])))
                    i_pair_decrease = f'{i_val[0]},{i_val[1]}'
                    if flag_assert:
                        assert np.equal(np.delete(para_np[i_, :], i), np.delete(para_np[i_matching_rows, :], i)).all()
                        assert para_np[i_, i] != para_np[i_matching_rows, i]

            if i_pair_increase is not None and i_pair_decrease is not None:
                sensitivity_set[i_name].append(i_pair_increase + ';' + i_pair_decrease)
            elif i_pair_increase is not None:
                sensitivity_set[i_name].append(i_pair_increase)
            elif i_pair_decrease is not None:
                sensitivity_set[i_name].append(i_pair_decrease)

        pbar.update(1)
        # if i_ % 10000 == 0:
        #     info = '{:.2e}'.format(i_/num_par)
        #     print(info)
        # print(i_)
        return

    if not os.path.exists(filename_sensitivity_set):
        assert stoptime['paramID'] == para['paramID']
        para_set = dict(zip(par_sensitivity, [None]*len(par_sensitivity)))
        stoptime_df = pd.DataFrame.from_dict(stoptime)
        if not os.path.exists(filename_sensitivity_para_df):
            para['S_01_2_S_01_1'] = [i/j for i, j in zip(para['S_01_2'], para['S_01_1'])]
            para['SR_00_1'] = [i/j for i, j in zip(para['S_00_1'], para['S_01_1'])]
            para['SR_00_2'] = [i/j for i, j in zip(para['S_00_2'], para['S_01_2'])]
            para['SR_11_1'] = [i/j for i, j in zip(para['S_11_1'], para['S_01_1'])]
            para['SR_11_2'] = [i/j for i, j in zip(para['S_11_2'], para['S_01_2'])]

            for i in par_sensitivity:
                i_par = list(set(para[i]))
                if i in ['SR_00_1', 'SR_00_2']:
                    i_par = ['{:.1e}'.format(x) for x in i_par]
                else:
                    i_par = ['{:.2e}'.format(x) for x in i_par]
                i_par = [float(x) for x in i_par]
                i_par = list(set(i_par))
                i_par.sort()
                para_set[i] = i_par

            para_set_num = {key: len(value) for key, value in para_set.items()}

            para_df = pd.DataFrame.from_dict(para)
            para_df = para_df[['paramID'] + par_sensitivity]
            para_df.drop_duplicates(subset=para_df.columns.to_list()[1:], inplace=True)
            with tqdm(total=len(para_df), ncols=120) as pbar:
                for i in range(para_df.shape[0]):
                    assert para_df.loc[i, 'paramID'] == stoptime_df.loc[i, 'paramID']
                    for j in par_sensitivity:
                        ij_val = para_df.loc[i, j]
                        i_ind = [isclose(ij_val, i_val, rel_tol=1e-2) for i_val in para_set[j]]
                        i_ind = DPM_reversiblemodel_miscellaneous_indextrueinlist(i_ind)
                        assert len(i_ind) == 1
                        para_df.loc[i, j] = i_ind[0]
                    pbar.update(1)
                with bz2.BZ2File(filename_sensitivity_para_df, 'wb') as f:
                    pickle.dump((para_df, para_set, para_set_num), f)
        else:
            with bz2.BZ2File(filename_sensitivity_para_df, 'rb') as f:
                para_df, para_set, para_set_num = pickle.load(f)

        flag_assert = True
        # para_df = para_df.iloc[:1500000, :]
        para_np = para_df.iloc[:, 1:].to_numpy(dtype=int)
        head_names = para_df.iloc[:, 1:].columns.to_list()
        ind_par_sensitivity_use = [head_names.index(element) for element in par_sensitivity_use2]
        ind_par_sensitivity_use = dict(zip(par_sensitivity_use2, ind_par_sensitivity_use))
        num_par = len(para_df)
        sensitivity_set = dict()
        str_format = ' ' * 30
        for i_key in par_sensitivity_use2:
            sensitivity_set[i_key] = list()

        num_cores = multiprocessing.cpu_count()
        with tqdm(total=len(para_df), ncols=120) as pbar:
            # Parallel(n_jobs=num_cores)(delayed(DPM_reversiblemodel_analysis_sensitivity4_1)(i) for i in range(len(para_df)))
            for i in range(num_par):
                DPM_reversiblemodel_analysis_sensitivity_nonαθμ_1(i)

        sensitivity_set_filter = dict(zip(par_sensitivity_use2, [list() for _ in range(len(par_sensitivity_use2))]))
        for i, i_key in enumerate(par_sensitivity_use2):
            i_sensitivity = sensitivity_set[i_key]
            with tqdm(total=len(i_sensitivity), ncols=120, desc=f'{i}, {i_key}') as pbar:
                for i_val in i_sensitivity:
                    if i_val != str_format and i_val not in sensitivity_set_filter[i_key]:
                        sensitivity_set_filter[i_key].append(i_val)
                    pbar.update(1)

        with bz2.BZ2File(filename_sensitivity_set, 'wb') as f:
            pickle.dump(sensitivity_set_filter, f)
    else:
        with bz2.BZ2File(filename_sensitivity_set, 'rb') as f:
            sensitivity_set = pickle.load(f)
    return


def DPM_reversiblemodel_analysis_sensitivityplot(stoptime,
                                                 strategy_name,
                                                 filename_sensitivity_para_df,
                                                 filename_sensitivity_set1,
                                                 filename_sensitivity_set2,
                                                 pathsave):
    def DPM_reversiblemodel_analysis_sensitivityplot_1(comb_, set_):
        comb_ = comb_.split(',')
        comb_ = [int(i_) for i_ in comb_]
        comb_ = list(itertools.combinations(comb_, 2))
        for i_combo in comb_:
            if tuple(sorted(i_combo)) not in set_:
                set_.append(i_combo)
        return set_

    filename_sensitivity = os.path.join(pathsave, 'sensitivity.pckl')
    filename_sensitivity_set_αθμ = os.path.join(pathsave, 'sensitivity_set_alphathetamu.pckl')
    filename_sensitivity_set_TSR1R2 = os.path.join(pathsave, 'sensitivity_set_TSR1R2.pckl')
    parname = par_sensitivity2 + par_sensitivity_use2
    result = dict(zip(strategy_name, [[] for _ in range(len(strategy_name))]))
    sensitivity = dict(zip(parname, [deepcopy(result) for _ in range(len(parname))]))
    with bz2.BZ2File(filename_sensitivity_para_df, 'rb') as f:
        para_df, para_set, para_set_num = pickle.load(f)

    if not os.path.exists(filename_sensitivity_set_αθμ):
        with bz2.BZ2File(filename_sensitivity_set1, 'rb') as f:
            sensitivity_αθμ = pickle.load(f)
        sensitivity_set_αθμ = list()
        with tqdm(total=len(sensitivity_αθμ), ncols=120, desc='αθμ') as pbar:
            for i_comb in sensitivity_set_αθμ:
                sensitivity_set_αθμ = DPM_reversiblemodel_analysis_sensitivityplot_1(i_comb, sensitivity_set_αθμ)
                pbar.update(1)
        with bz2.BZ2File(filename_sensitivity_set_αθμ, 'wb') as f:
            pickle.dump(sensitivity_set_αθμ, f)
    else:
        with bz2.BZ2File(filename_sensitivity_set_αθμ, 'rb') as f:
            sensitivity_set_αθμ = pickle.load(f)

    if not os.path.exists(filename_sensitivity_set_TSR1R2):
        with bz2.BZ2File(filename_sensitivity_set2, 'rb') as f:
            sensitivity_par = pickle.load(f)

        sensitivity_set_TSR1R2 = dict()
        for i_parname in par_sensitivity_use2:
            sensitivity_set_TSR1R2[i_parname] = sensitivity_par[i_parname]

        for i_name in par_sensitivity_use2:
            i_set = list(sensitivity_set_TSR1R2[i_name])
            i_set_filter = list()
            with tqdm(total=len(i_set), ncols=120, desc=f'{i_name}') as pbar:
                for i_comb in i_set:
                    if ';' not in i_comb:
                        i_set_filter = DPM_reversiblemodel_analysis_sensitivityplot_1(i_comb, i_set_filter)
                    else:
                        i_comb = i_comb.split(';')
                        for j_comb in i_comb:
                            i_set_filter = DPM_reversiblemodel_analysis_sensitivityplot_1(j_comb, i_set_filter)
                    pbar.update(1)
                sensitivity_set_TSR1R2[i_name] = i_set_filter

        with bz2.BZ2File(filename_sensitivity_set_TSR1R2, 'wb') as f:
            pickle.dump(sensitivity_set_TSR1R2, f)
    else:
        with bz2.BZ2File(filename_sensitivity_set_TSR1R2, 'rb') as f:
            sensitivity_set_TSR1R2 = pickle.load(f)

    if not os.path.exists(filename_sensitivity):
        sensitivity_set = dict()
        sensitivity_set[par_sensitivity2[0]] = sensitivity_set_αθμ
        for i_parname in par_sensitivity_use2:
            sensitivity_set[i_parname] = sensitivity_set_TSR1R2[i_parname]

        para_df = para_df.iloc[:, 1:]
        head_names = para_df.columns.to_list()
        numpar = list(range(len(head_names)))
        ind_par_sensitivity_use = [head_names.index(element) for element in par_sensitivity_use2]
        ind_par_sensitivity_use = dict(zip(par_sensitivity_use2, ind_par_sensitivity_use))
        ind_par_sensitivity_use['αθμ'] = [0, 1, 2, 3, 4, 5]
        for i_name in parname:
            i_set = sensitivity_set[i_name]
            i_index_diff = ind_par_sensitivity_use[i_name]
            if type(i_index_diff) == int:
                i_index_diff = list([i_index_diff])
            i_index_same = [x for x in numpar if x not in i_index_diff]
            with tqdm(total=len(i_set), ncols=120, desc=f'{i_name}') as pbar:
                for i_comb in i_set:
                    assert np.all(para_df.iloc[i_comb[0], i_index_same].values == para_df.iloc[i_comb[1], i_index_same].values)
                    assert np.isin(False, para_df.iloc[i_comb[0], i_index_diff].values == para_df.iloc[i_comb[1], i_index_diff].values)
                    pbar.update(1)

        for i_parname in parname:
            i_comb = sensitivity_set[i_parname]
            if len(i_comb) != 0:
                i_comb = list(map(list, zip(*i_comb)))
                for i_strategy in strategy_name:
                    i_stoptime = np.array(stoptime[i_strategy])
                    sensitivity[i_parname][i_strategy] = np.abs(i_stoptime[i_comb[0]] - i_stoptime[i_comb[1]])

        with bz2.BZ2File(filename_sensitivity, 'wb') as f:
            pickle.dump(sensitivity, f)
    else:
        with bz2.BZ2File(filename_sensitivity, 'rb') as f:
            sensitivity = pickle.load(f)
        df = pd.DataFrame(columns=['Strategy', 'para', 'val'])

        parname.remove('SR_00_1')
        for i_name in parname:
            for i_strategy in strategy_name[:-1]:
                i_sensitivity = sensitivity[i_name][i_strategy]
                i_len = i_sensitivity.shape[0]
                i_df = pd.DataFrame({'para': [i_name] * i_len, 'Strategy': [i_strategy] * i_len, 'val': i_sensitivity})
                df = pd.concat([df, i_df], ignore_index=True)

        ymin = -5
        ymax = 220
        # αθμ
        fig, ax = plt.subplots()
        boxplot(data=df[df['para'] == parname[0]], x='Strategy', y='val', ax=ax, showfliers=False, fill=False)
        fig.set_size_inches(12.5, 5)
        ax.set_ylim([ymin, ymax])
        plt.tight_layout()
        # plt.yscale("log")
        plt.savefig(os.path.join(pathsave, 'αθμ'+'.pdf'), format='pdf', bbox_inches='tight')
        plt.close()

        color = ["#DAA520", "#000000", "#0000FF", "#FF0000", "#008000", "#FF82AB", "#FF82AB", "#8A2BE2", "#8A2BE2", "#A9A9A9"]
        my_pal = dict(zip(strategy_name, color))

        # T1 and T2
        df_sel = df[df['para'].isin(['T1', 'T2'])]
        fig, ax = plt.subplots()
        ax_ = boxplot(data=df_sel, x='para', y='val', hue='Strategy', gap=0.15, ax=ax, palette=my_pal,  showfliers=False, fill=False)
        move_legend(ax_, 'upper left')
        fig.set_size_inches(12.5, 5)
        ax.set_ylim([ymin, 100])
        plt.tight_layout()
        # plt.yscale("log")
        plt.savefig(os.path.join(pathsave, 'T' + '.pdf'), format='pdf', bbox_inches='tight')
        plt.close()

        # S2/S1
        fig, ax = plt.subplots()
        boxplot(data=df[df['para'] == parname[3]], x='Strategy', y='val', ax=ax, showfliers=False, fill=False)
        fig.set_size_inches(12.5, 5)
        ax.set_ylim([ymin, ymax])
        plt.tight_layout()
        plt.savefig(os.path.join(pathsave, 'S_01_2ratio' + '.pdf'), format='pdf', bbox_inches='tight')
        plt.close()

        # SR
        df_sel = df[df['para'].isin(['SR_11_1', 'SR_11_2'])]  # 'SR_00_1'
        fig, ax = plt.subplots()
        ax_ = boxplot(data=df_sel, x='para', y='val', hue='Strategy', gap=0.15, ax=ax, palette=my_pal,  showfliers=False, fill=False)
        ax_.legend([], [], frameon=False)
        fig.set_size_inches(12.5, 5)
        ax.set_ylim([-1, 50])
        plt.tight_layout()
        plt.savefig(os.path.join(pathsave, 'SR' + '.pdf'), format='pdf', bbox_inches='tight')
        plt.close()

        # R1ratio and R2ratio
        df_sel = df[df['para'].isin(['R1ratio', 'R2ratio'])]  # 'SR_00_1'
        fig, ax = plt.subplots()
        ax_ = boxplot(data=df_sel, x='para', y='val', hue='Strategy', gap=0.15, ax=ax, palette=my_pal,  showfliers=False, fill=False)
        ax_.legend([], [], frameon=False)
        fig.set_size_inches(12.5, 5)
        ax.set_ylim([ymin, 300])
        plt.tight_layout()
        plt.savefig(os.path.join(pathsave, 'Rratio' + '.pdf'), format='pdf', bbox_inches='tight')
        plt.close()
    return


def DPM_reversiblemodel_analysis_sankeydiagram(para_source,
                                               para_name_source,
                                               para_target,
                                               category_source,
                                               category_target,
                                               name,
                                               pathsave_figure_sankeydiagram):
    def DPM_reversiblemodel_analysis_sankeydiagram1_1(i_, i_par, i_category):
        assert i_par['paramID'] == i_category['paramID']
        if name == 'mut to full':
            i_paramID = int(i_par['paramID'])
            i_index = paramID_source.index(i_paramID)
            i_par_source = para_source[i_index]
            i_flag = True
            for j_key in para_name_source:
                j_flag = isclose(i_par[j_key], i_par_source[j_key], rel_tol=1e-12)
                i_flag = i_flag and j_flag
            assert i_flag
            sankey_df.loc[category_source[i_index]['category'], i_category['category']] = \
                sankey_df.loc[category_source[i_index]['category'], i_category['category']] + 1
        elif name == 'rever to full':
            num_source = 0
            for j in range(len(para_source)):
                j_par = para_source[j]
                j_category = category_source[j]
                assert j_par['paramID'] == j_category['paramID']
                j_flag = True
                for k_key in para_name_source:
                    k_flag = isclose(i_par[k_key], j_par[k_key], rel_tol=1e-9)
                    j_flag = j_flag and k_flag
                if j_flag:
                    num_source = num_source + 1
                    sankey_df.loc[j_category['category'], i_category['category']] = \
                        sankey_df.loc[j_category['category'], i_category['category']] + 1
            assert num_source == 1 and sankey_df.values.sum() == i_ + 1
        pbar.update(1)
        return

    filename = os.path.join(pathsave_figure_sankeydiagram, name)
    filename = filename + '.pckl'
    if not os.path.isfile(filename):
        paramID_source = para_source['paramID']
        para_source = [dict(zip(para_source, x)) for x in zip(*para_source.values())]
        para_target = [dict(zip(para_target, x)) for x in zip(*para_target.values())]
        category_source = [dict(zip(category_source, x)) for x in zip(*category_source.values())]
        category_target = [dict(zip(category_target, x)) for x in zip(*category_target.values())]

        source_set = ["".join(seq) for seq in itertools.product("01", repeat=len(category_source[0]['category']))]
        target_set = ["".join(seq) for seq in itertools.product("01", repeat=4)]
        sankey_df = pd.DataFrame(0, index=pd.Index(source_set), columns=pd.Index(target_set))

        with tqdm(total=len(para_target), ncols=120, desc='running') as pbar:
            for i in range(len(para_target)):
                DPM_reversiblemodel_analysis_sankeydiagram1_1(i, para_target[i], category_target[i])
        assert sankey_df.values.sum() == len(para_target)
        with bz2.BZ2File(filename, 'wb') as f:
            pickle.dump(sankey_df, f)
    else:
        with bz2.BZ2File(filename, 'rb') as f:
            sankey_df = pickle.load(f)
    return sankey_df
