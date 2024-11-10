from DPM_reversiblemodel_constant import *
from DPM_reversiblemodel_model import \
    DPM_reversiblemodel_model_ODE, \
    DPM_reversiblemodel_model_simexpm, \
    DPM_reversiblemodel_model_prorate, \
    DPM_reversiblemodel_model_transitionrate, \
    DPM_reversiblemodel_fun_model_flow
from DPM_reversiblemodel_strategy import \
     DPM_reversiblemodel_strategy_simbydose, \
     DPM_reversiblemodel_strategy_bestcycle, \
     DPM_reversiblemodel_strategy_bestcyclesim, \
     DPM_reversiblemodel_strategy_0, \
     DPM_reversiblemodel_strategy_1, \
     DPM_reversiblemodel_strategy_2, \
     DPM_reversiblemodel_strategy_random
from DPM_reversiblemodel_miscellaneous import \
    DPM_reversiblemodel_stime, \
    DPM_reversiblemodel_cure, \
    DPM_reversiblemodel_miscellaneous_writehead, \
    DPM_reversiblemodel_miscellaneous_cycletreat2dose, \
    DPM_reversiblemodel_miscellaneous_dictlist, \
    DPM_reversiblemodel_miscellaneous_treat_mono_or_half
from DPM_reversiblemodel_readsave import DPM_reversiblemodel_readsave_save_csv, \
    DPM_reversiblemodel_readsave_stoptime_csv, \
    DPM_reversiblemodel_readsave_dosage_csv
from DPM_reversiblemodel_analysis import DPM_reversiblemodel_analysis_KM, \
    DPM_reversiblemodel_analysis_hz, \
    DPM_reversiblemodel_analysis_stoptime, \
    DPM_reversiblemodel_analysis_dose, \
    DPM_reversiblemodel_analysis_Nflow, \
    DPM_reversiblemodel_analysis_sensitivity, \
    DPM_reversiblemodel_analysis_sensitivity_αθμ, \
    DPM_reversiblemodel_analysis_sensitivity_nonαθμ, \
    DPM_reversiblemodel_analysis_sensitivityplot, \
    DPM_reversiblemodel_analysis_sankeydiagram
from DPM_reversiblemodel_plot import DPM_reversiblemodel_plot_km_multi, \
    DPM_reversiblemodel_plot_df2pdf, \
    DPM_reversiblemodel_plot_1strategy, \
    DPM_reversiblemodel_plot_multi, \
    DPM_reversiblemodel_plot_sankeydiagram


def DPM_reversiblemodel_fun_events(t, y, *_):
    if t is not None:
        event = np.sum(y, axis=0) - Limit_mortality
    else:
        event = 1
    return event


def DPM_reversiblemodel_fun_stopX(Nt, N_nonmut, N_mut, t, maxthreshold):
    if np.any(Nt >= maxthreshold):
        pos = np.argmax(Nt >= maxthreshold)
        Nt, N_nonmut, N_mut, t = Nt[:pos+1], N_nonmut[:pos+1], N_mut[:pos+1], t[:pos+1]
    return Nt, N_nonmut, N_mut, t


def DPM_reversiblemodel_fun_steadystate(pars, mutation, eachtimestep=False):
    t = np.arange(0, duration_steady + Simtimestep, Simtimestep)

    tsim = t if eachtimestep else [0, duration_steady]
    argsODE = (σ_notreat, mutation, pars)
    N = DPM_reversiblemodel_model_simexpm(argsODE, tsim, x0_init, LSsim=False)[0]

    ratio_0101 = N[ind_0101, -1] / np.sum(N[:, -1], 0)
    ratio_0001 = N[ind_0001, -1] / np.sum(N[:, -1], 0)
    ratio_0100 = N[ind_0100, -1] / np.sum(N[:, -1], 0)
    ratio_0000 = N[ind_0000, -1] / np.sum(N[:, -1], 0)
    if pars['μ1'] != 0 or pars['μ2'] != 0:
        assert ratio_0001 != 0
        assert ratio_0100 != 0
        assert ratio_0000 != 0

    return N, t, {'ratio_0101': ratio_0101, 'ratio_0001': ratio_0001, 'ratio_0100': ratio_0100, 'ratio_0000': ratio_0000}


def DPM_reversiblemodel_fun_initx(n_ratio, pars):
    n_nonmut = n_total * (1-pars['R1ratio']-pars['R2ratio'])

    n0101 = n_nonmut * n_ratio['ratio_0101']
    n0001 = n_nonmut * n_ratio['ratio_0001']
    n0100 = n_nonmut * n_ratio['ratio_0100']
    n0000 = n_nonmut * n_ratio['ratio_0000']

    n1101 = n_total * pars['R1ratio'] * n_ratio['ratio_0101']/(n_ratio['ratio_0101'] + n_ratio['ratio_0100'])
    n1100 = n_total * pars['R1ratio'] * n_ratio['ratio_0100']/(n_ratio['ratio_0101'] + n_ratio['ratio_0100'])

    n0111 = n_total * pars['R2ratio'] * n_ratio['ratio_0101']/(n_ratio['ratio_0101'] + n_ratio['ratio_0001'])
    n0011 = n_total * pars['R2ratio'] * n_ratio['ratio_0001']/(n_ratio['ratio_0101'] + n_ratio['ratio_0001'])

    n1111 = 0.
    return [n0101, n0001, n0100, n0000, n0111, n0011, n1101, n1100, n1111]


def DPM_reversiblemodel_fun_sim1par(strategy_run, par, pathsave, mutonly=False, reveronly=False):
    N, tmax, duration = [None]*len(strategy_run), 0, duration_5year
    par = pd.DataFrame.from_dict(par, orient='index').T
    par.rename(columns=dict(zip(pd_colname, par_key)), errors='raise', inplace=True)
    par = par.to_dict(orient='records')[0]
    mutation = False
    n_ratio = DPM_reversiblemodel_fun_steadystate(par, mutation)[2]
    x0 = DPM_reversiblemodel_fun_initx(n_ratio, par)
    mutation, LSsim, t = True, True, np.arange(0, duration + Simtimestep, Simtimestep)
    tdiff = np.append(0, np.diff(t))
    # Mono σ1
    argsODE = (σ1_full, mutation, par)
    N_mono_drug1, t_mono_drug1, *_ = DPM_reversiblemodel_model_simexpm(argsODE, tdiff, x0)
    d_mono_drug1 = np.tile(np.array([σ1_full]).T, (1, N_mono_drug1.shape[1] - 1))
    # Mono σ2
    argsODE = (σ2_full, mutation, par)
    N_mono_drug2, t_mono_drug2, *_ = DPM_reversiblemodel_model_simexpm(argsODE, tdiff, x0)
    d_mono_drug2 = np.tile(np.array([σ2_full]).T, (1, N_mono_drug2.shape[1] - 1))
    # Mono half
    argsODE = (σ_half, mutation, par)
    N_mono_half, t_mono_half, *_ = DPM_reversiblemodel_model_simexpm(argsODE, tdiff, x0)
    d_mono_half = np.tile(np.array([σ_half]).T, (1, N_mono_half.shape[1]-1))
    # Best cycle
    N_bestcycle, d_bestcycle, σ_bestcycle = DPM_reversiblemodel_strategy_bestcyclesim(x0, par, duration_5year, mutation)
    # Strategy 0
    N_strategy0, d_strategy0 = DPM_reversiblemodel_strategy_0(x0, par, duration, mutation)
    # Strategy 1
    N_strategy1, d_strategy1 = DPM_reversiblemodel_strategy_1(x0, par, duration, mutation)
    # Strategy 1 cycle
    N_strategy1_cycle, d_strategy1_cycle = DPM_reversiblemodel_strategy_1(x0, par, duration, mutation, use_cycle=True)
    # Strategy 2
    N_strategy2, d_strategy2 = DPM_reversiblemodel_strategy_2(x0, par, duration, mutation)
    # strategy2 cycle
    N_strategy2_cycle, d_strategy2_cycle = DPM_reversiblemodel_strategy_2(x0, par, duration, mutation, use_cycle=True)

    if 'mono 1' in strategy_run:
        pos = strategy_run.index('mono 1')
        N[pos] = N_mono_drug1
        tmax = N_mono_drug1.shape[1] if N_mono_drug1.shape[1] > tmax else tmax
    if 'mono 2' in strategy_run:
        pos = strategy_run.index('mono 2')
        N[pos] = N_mono_drug2
        tmax = N_mono_drug2.shape[1] if N_mono_drug2.shape[1] > tmax else tmax
    if 'mono half' in strategy_run:
        pos = strategy_run.index('mono half')
        N[pos] = N_mono_half
        tmax = N_mono_half.shape[1] if N_mono_half.shape[1] > tmax else tmax
    if 'best cycle' in strategy_run:
        pos = strategy_run.index('best cycle')
        N[pos] = N_bestcycle
        tmax = N_bestcycle.shape[1] if N_bestcycle.shape[1] > tmax else tmax
    if 'strategy0' in strategy_run:
        pos = strategy_run.index('strategy0')
        N[pos] = N_strategy0
        tmax = N_strategy0.shape[1] if N_strategy0.shape[1] > tmax else tmax
    if 'strategy1' in strategy_run:
        pos = strategy_run.index('strategy1')
        N[pos] = N_strategy1
        tmax = N_strategy1.shape[1] if N_strategy1.shape[1] > tmax else tmax
    if 'strategy1 cycle' in strategy_run:
        pos = strategy_run.index('strategy1 cycle')
        N[pos] = N_strategy1_cycle
        tmax = N_strategy1_cycle.shape[1] if N_strategy1_cycle.shape[1] > tmax else tmax
    if 'strategy2' in strategy_run:
        pos = strategy_run.index('strategy2')
        N[pos] = N_strategy2
        tmax = N_strategy2.shape[1] if N_strategy2.shape[1] > tmax else tmax
    if 'strategy2 cycle' in strategy_run:
        pos = strategy_run.index('strategy2 cycle')
        N[pos] = N_strategy2_cycle
        tmax = N_strategy2_cycle.shape[1] if N_strategy2_cycle.shape[1] > tmax else tmax

    if not os.path.exists(pathsave):
        os.makedirs(pathsave)
    strategy_name = strategy_run[0]
    # plot
    # mono 1
    if 'mono 1' in strategy_run:
        title = strategy_name + ': mono 1'
        DPM_reversiblemodel_plot_1strategy(N_mono_drug1, d_mono_drug1, tmax, title, pathsave, mutonly, reveronly)
    # mono 2
    if 'mono 2' in strategy_run:
        title = strategy_name + ': mono 2'
        DPM_reversiblemodel_plot_1strategy(N_mono_drug2, d_mono_drug2, tmax, title, pathsave, mutonly, reveronly)
    # cycle
    if 'mono half' in strategy_run:
        title = strategy_name + ': mono half'
        DPM_reversiblemodel_plot_1strategy(N_mono_half, d_mono_half, tmax, title, pathsave, mutonly, reveronly)
    if 'best cycle' in strategy_run:
        title = strategy_name + ': cycle ' + σ_bestcycle
        DPM_reversiblemodel_plot_1strategy(N_bestcycle, d_bestcycle, tmax, title, pathsave, mutonly, reveronly)
    if 'strategy0' in strategy_run:
        title = strategy_name + ': strategy0'
        DPM_reversiblemodel_plot_1strategy(N_strategy0, d_strategy0, tmax, title, pathsave, mutonly, reveronly)
    if 'strategy1' in strategy_run:
        title = strategy_name + ': strategy1'
        DPM_reversiblemodel_plot_1strategy(N_strategy1, d_strategy1, tmax, title, pathsave, mutonly, reveronly)
    if 'strategy1 cycle' in strategy_run:
        title = strategy_name + ': strategy1cycle'
        DPM_reversiblemodel_plot_1strategy(N_strategy1_cycle, d_strategy1_cycle, tmax, title, pathsave, mutonly, reveronly)
    if 'strategy2' in strategy_run:
        title = strategy_name + ': strategy2'
        DPM_reversiblemodel_plot_1strategy(N_strategy2, d_strategy2, tmax, title, pathsave, mutonly, reveronly)
    if 'strategy2 cycle' in strategy_run:
        title = strategy_name + ': strategy2cycle'
        DPM_reversiblemodel_plot_1strategy(N_strategy2_cycle, d_strategy2_cycle, tmax, title, pathsave, mutonly, reveronly)
    # plot multiple
    # total
    ind = ind_nonmut + ind_mut
    title = strategy_name + ': total'
    DPM_reversiblemodel_plot_multi(N, ind, strategy_run, title)
    plt.savefig(os.path.join(pathsave, 'total.pdf'), format='pdf', bbox_inches='tight')
    plt.close()
    # n1111
    if not reveronly:
        ind = ind_multi_mut
        title = strategy_name + ': multi_mut'
        DPM_reversiblemodel_plot_multi(N, ind, strategy_run, title)
        plt.savefig(os.path.join(pathsave, 'multi_mut.pdf'), format='pdf', bbox_inches='tight')
        plt.close()
    return


def DPM_reversiblemodel_fun_sim(x0, pars, σ_all, duration_all, mutation, LSsim, last=False):
    if len(σ_all) != len(duration_all):
        raise Exception("Length of σ and duration is not the same.")
    t = np.zeros([0], dtype=int)
    σ = np.zeros([Num_drug, 0], dtype='float')
    tnow = 0
    for i, i_σ in enumerate(σ_all):
        i_duration = duration_all[i]
        i_t = np.arange(tnow, tnow + i_duration + Simtimestep, Simtimestep)
        i_σ = np.tile(np.array([i_σ]).T, (1, i_t.shape[0]-1))
        t = np.append(t, i_t + tnow) if tnow == 0 else np.append(t, i_t[1:])
        σ = np.append(σ, i_σ, 1)
        tnow = i_t[-1]

    if LSsim:
        event = DPM_reversiblemodel_fun_events
        event.direction, event.terminal = 1, True
    else:
        event = None
    if last:
        t_eval = t[[0, -1]]
    else:
        t_eval = t
    argsODE = (σ, mutation, pars)
    sim = solve_ivp(DPM_reversiblemodel_model_ODE, t[[0, -1]],
                    list(x0),
                    method=ODE_method,
                    rtol=SCIPY_INTEGRATE_SOLVE_IVP['RTOL'],
                    atol=SCIPY_INTEGRATE_SOLVE_IVP['ATOL'],
                    t_eval=t_eval,
                    events=event,
                    args=argsODE)

    N_total = np.sum(sim.y, axis=0)
    N_nonmut = np.sum(sim.y[ind_nonmut, :], axis=0)
    N_mut = np.sum(sim.y[ind_mut, :], axis=0)

    return sim.y, t, N_total, N_nonmut, N_mut, sim


def DPM_reversiblemodel_fun_generatepar_fromPNAS(pathload='./pnas',
                                                 pathsave_mut='./reversible_parcsv_mut',
                                                 filename_pattern='_para',
                                                 joint=False):
    def DPM_reversiblemodel_fun_generatepar_fromPNAS_1(i_file_):
        df = pd.read_csv(i_file_)
        par_csv = df.to_dict(orient='records')
        return par_csv

    def DPM_reversiblemodel_fun_generatepar_fromPNAS_2(par_):
        keys = ['paramID']+pd_colname
        par = dict(zip(keys, [0]*len(keys)))
        par['paramID'] = int(par_['paramID'])
        assert par_['g0_S'] == par_['g0_R1'] == par_['g0_R2'] == par_['g0_R12']
        totalnum_ = par_['Spop'] + par_['R1pop'] + par_['R2pop'] + par_['R12pop']
        par['g'] = par_['g0_S']
        par['R1ratio'] = par_['R1pop']/totalnum_
        par['R2ratio'] = par_['R2pop']/totalnum_
        par['T1'], par['T2'] = par_['T.R1..S.'], par_['T.R2..S.']
        assert par_['T.R12..R1.'] == par_['T.R2..S.'] and par_['T.R12..R2.'] == par_['T.R1..S.']
        par['S_01_1'] = par_['Sa.S.D1.']
        par['S_01_2'] = par_['Sa.S.D2.']
        par['S_11_1'] = par_['Sa.R1.D1.']
        par['S_11_2'] = par_['Sa.R2.D2.']
        par['S_00_1'] = 0
        par['S_00_2'] = 0
        assert par['S_01_1'] >= par['S_01_2']
        assert par_['Sa.S.D2.'] == par_['Sa.R1.D2.'] and par_['Sa.S.D1.'] == par_['Sa.R2.D1.']
        assert par_['Sa.R12.D1.'] == par_['Sa.R1.D1.'] and par_['Sa.R12.D2.'] == par_['Sa.R2.D2.']
        return par

    def DPM_reversiblemodel_fun_generatepar_fromPNAS_3(par_i, rate_rever_):
        αθκ_i, μ_i = par_i[0], 10**par_i[1]

        Storr_σ0_i, Storr_σhalf_i, Storr_σfull_i, rrtoS_σ0_i, rrtoS_σhalf_i, rrtoS_σfull_i \
            = DPM_reversiblemodel_model_transitionrate(*αθκ_i, μ_i)

        '''criteria: no drug, rate to sensitive state must higher than rate to reversible resistant state, 
           bigger than largest muatation rate and smaller than 0.1
           full dose, rate to reversible resistant state is bigger than largest muatation rate and smaller than 0.1'''
        # 1 >= rrtoS_σ0_i >= Storr_σfull_i and rrtoS_σ0_i >= 1e-2 and Storr_σhalf_i >= 5e-4 and 5e-2 >= rrtoS_σ0_i >= 5e-3
        flag_sel = (rrtoS_σ0_i > Storr_σ0_i and
                    1e-1 >= rrtoS_σ0_i >= 1e-3 and
                    1e-1 >= Storr_σfull_i >= 1e-3)
        if flag_sel:
            par_α_θ_μ.append(αθκ_i + (μ_i,))
            rate_rever_i = pd.DataFrame.from_dict(dict(zip(columns, [[αθκ_i[0]], [αθκ_i[1]], [αθκ_i[2]],
                                                                     [μ_i], [Storr_σ0_i], [rrtoS_σ0_i],
                                                                     [Storr_σfull_i], [rrtoS_σfull_i]])), orient='columns')
            rate_rever_ = pd.concat([rate_rever_, rate_rever_i])
        pbar.update(1)
        return rate_rever_

    def DPM_reversiblemodel_fun_generatepar_fromPNAS_4(S_00_i, par_i):
        SR_00_1_i, SR_00_2_i = S_00_i[0], S_00_i[0]
        g_i, S_01_1_i, S_01_2_i, S_11_1_i, S_11_2_i = par_i['g'], par_i['S_01_1'], par_i['S_01_2'], par_i['S_11_1'], par_i['S_11_2']
        S_00_1_i, S_00_2_i = par_i['S_01_1'] * SR_00_1_i, par_i['S_01_2'] * SR_00_2_i

        rate_i = DPM_reversiblemodel_model_prorate(g_i, S_01_1_i, S_01_2_i, S_00_1_i, S_00_2_i, S_11_1_i, S_11_2_i)

        assert S_01_1_i > S_00_1_i and S_01_1_i > S_11_1_i
        assert S_01_2_i > S_00_2_i and S_01_2_i > S_11_2_i
        assert S_01_1_i >= S_01_2_i
        assert rate_i['IR1 σ1 full'] >= 0 and rate_i['IR2 σ2 full'] >= 0

        flag1 = S_00_1_i > S_11_1_i and S_00_2_i > S_11_2_i
        flag2 = par_i['S_11_1']/par_i['S_01_1'] < 0.1 and par_i['S_11_2']/par_i['S_01_2'] < 0.1
        flag3 = par_i['R1ratio'] > 0 and par_i['R2ratio'] > 0
        flag4 = par_i['R1ratio'] < 0.8 and par_i['R2ratio'] < 0.8
        flag5 = par_i['S_01_1']/i_par['g'] > 0.9 and par_i['S_01_2']/i_par['g'] > 0.9
        flag6 = par_i['g'] < 0.29

        if all([flag1, flag2, flag3, flag4, flag5, flag6]):
            par_i['S_00_1'], par_i['S_00_2'] = S_00_1_i, S_00_2_i
            i_paramID = str(par_i['paramID']) + '.' + str(paramID_count[str(par_i['paramID'])])
            paramID_count[str(par_i['paramID'])] += 1
            par_i['paramID'] = i_paramID
            assert par_i['S_00_1'] > par_i['S_11_1']
            assert par_i['S_00_2'] > par_i['S_11_2']
            PAR.append(par_i)
        pbar.update(1)
        return

    def DPM_reversiblemodel_fun_generatepar_fromPNAS_5(_, filename):
        ind = filename.split('_')[2]
        with bz2.BZ2File(os.path.join(pathsave_parraw, filename), 'rb') as f:
            par_no_α_θ_μ = pickle.load(f)
        par = []
        block_count_ = 1

        totalnum_ = len(par_α_θ_μ) * len(par_no_α_θ_μ)
        totalcombination_ = itertools.product(par_α_θ_μ, par_no_α_θ_μ)

        with tqdm(total=totalnum_, ncols=150, desc='generating par') as pbar_i:
            for par_i in totalcombination_:
                i_par_α_θ_μ_drug1, i_par_α_θ_μ_drug2, i_par_no_α_θ_μ = par_i[0], par_i[0], par_i[1]
                i_par_no_α_θ_μ['alpha1'] = i_par_α_θ_μ_drug1[0]
                i_par_no_α_θ_μ['theta1'] = i_par_α_θ_μ_drug1[1]
                i_par_no_α_θ_μ['kappa1'] = i_par_α_θ_μ_drug1[2]
                i_par_no_α_θ_μ['mu1'] = i_par_α_θ_μ_drug1[3]

                i_par_no_α_θ_μ['alpha2'] = i_par_α_θ_μ_drug2[0]
                i_par_no_α_θ_μ['theta2'] = i_par_α_θ_μ_drug2[1]
                i_par_no_α_θ_μ['kappa2'] = i_par_α_θ_μ_drug2[2]
                i_par_no_α_θ_μ['mu2'] = i_par_α_θ_μ_drug2[3]

                assert i_par_no_α_θ_μ['alpha1'] == i_par_no_α_θ_μ['alpha2']
                assert i_par_no_α_θ_μ['theta1'] == i_par_no_α_θ_μ['theta2']
                assert i_par_no_α_θ_μ['kappa1'] == i_par_no_α_θ_μ['kappa2']
                assert i_par_no_α_θ_μ['mu1'] == i_par_no_α_θ_μ['mu2']

                for i_, i_name in enumerate(pd_colname[:8]):
                    i_par_no_α_θ_μ[par_key[i_]] = i_par_no_α_θ_μ.pop(i_name)

                flag_i = DPM_reversiblemodel_fun_generatepar_fromPNAS_6(i_par_no_α_θ_μ)

                if all(flag_i):
                    assert 'paramID' in i_par_no_α_θ_μ.keys()
                    assert i_par_no_α_θ_μ['S_00_1'] > i_par_no_α_θ_μ['S_11_1']
                    assert i_par_no_α_θ_μ['S_00_2'] > i_par_no_α_θ_μ['S_11_2']
                    i_SR_00_1 = i_par_no_α_θ_μ['S_00_1']/i_par_no_α_θ_μ['S_01_1']
                    i_SR_00_2 = i_par_no_α_θ_μ['S_00_2']/i_par_no_α_θ_μ['S_01_2']
                    assert isclose(i_SR_00_1, i_SR_00_2, rel_tol=1e-9)
                    par.append(i_par_no_α_θ_μ)

                if len(par) == par_save_block_size:
                    filename = os.path.join(pathsave_par, f'par_reversiblemodel_{ind}_' + block_format.format(block_count_) +
                                            '_{:d}'.format(len(par)) + '.pckl')
                    with bz2.BZ2File(filename, 'wb') as f:
                        pickle.dump(par, f)
                    par = []
                    block_count_ += 1

                pbar_i.update(1)

            if len(par) != 0:
                filename = os.path.join(pathsave_par, f'par_reversiblemodel_{ind}_' + block_format.format(block_count_) +
                                        '_{:d}'.format(len(par)) + '.pckl')

                with bz2.BZ2File(filename, 'wb') as f:
                    pickle.dump(par, f)
        return

    def DPM_reversiblemodel_fun_generatepar_fromPNAS_6(par_i):
        α1_i, θ1_i, κ1_i, μ1_i = par_i['α1'], par_i['θ1'], par_i['κ1'], par_i['μ1']
        α2_i, θ2_i, κ2_i, μ2_i = par_i['α2'], par_i['θ2'], par_i['κ2'], par_i['μ2']
        T1_i = par_i['T1']
        T2_i = par_i['T2']

        _, Storr_σ1half_i, Storr_σ1full_i, _, _, _ = DPM_reversiblemodel_model_transitionrate(α1_i, θ1_i, κ1_i, μ1_i)
        _, Storr_σ2half_i, Storr_σ2full_i, _, _, _ = DPM_reversiblemodel_model_transitionrate(α2_i, θ2_i, κ2_i, μ2_i)

        assert Storr_σ1full_i > Storr_σ1half_i and Storr_σ2full_i > Storr_σ2half_i
        if T1_i > min([Storr_σ1full_i, Storr_σ1half_i]) or T2_i > min([Storr_σ2full_i, Storr_σ2half_i]):
            return False, False, False

        mutation = False
        _, _, n_ratio = DPM_reversiblemodel_fun_steadystate(par_i, mutation)
        x0 = DPM_reversiblemodel_fun_initx(n_ratio, par_i)

        mutation, LSsim = True, True
        t = np.arange(0, duration_5year + Simtimestep, Simtimestep)
        tdiff = np.append(0, np.diff(t))

        ''' σ1 full dose '''
        argsODE = (σ1_full, mutation, par_i)
        N_σ1_full, *_ = DPM_reversiblemodel_model_simexpm(argsODE, tdiff, x0)
        flag_mono_σ1_cure = DPM_reversiblemodel_cure(N_σ1_full)
        survi_time_mono_σ1 = DPM_reversiblemodel_stime(N_σ1_full)
        flag_mono_σ1 = True if not flag_mono_σ1_cure and survi_time_mono_σ1 < duration_5year else False
        if not flag_mono_σ1:
            return False, False, False

        '''σ2 full dose'''
        argsODE = (σ2_full, mutation, par_i)
        N_σ2_full, *_ = DPM_reversiblemodel_model_simexpm(argsODE, tdiff, x0)
        flag_mono_σ2_cure = DPM_reversiblemodel_cure(N_σ2_full)
        survi_time_mono_σ2 = DPM_reversiblemodel_stime(N_σ2_full)
        flag_mono_σ2 = True if not flag_mono_σ2_cure and survi_time_mono_σ2 < duration_5year else False
        if not flag_mono_σ2:
            return False, False, False

        '''σ1=0.5 and σ2=0.5'''
        argsODE = (σ_half, mutation, par_i)
        N_σ_half, *_ = DPM_reversiblemodel_model_simexpm(argsODE, tdiff, x0)
        flag_half_σ_cure = DPM_reversiblemodel_cure(N_σ_half)
        survi_time_half_σ = DPM_reversiblemodel_stime(N_σ_half)
        flag_half_σ = True if not flag_half_σ_cure and survi_time_half_σ < duration_5year else False
        return flag_mono_σ1, flag_mono_σ2, flag_half_σ

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%---begin---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
    if not os.path.exists(pathsave_mut):
        os.makedirs(pathsave_mut)
    if not os.path.exists(pathsave_par):
        os.makedirs(pathsave_par)
    if not joint:
        file_list = os.listdir(pathload)
        file_list = [filename for filename in file_list if re.search(filename_pattern, filename) is not None]
        paramID, block_count, PAR = 0, 1, []
        try:
            file_list.sort(key=lambda x: int(x.split('_')[-3]))
        except (ValueError, IndexError):
            pass
        file_list = [os.path.join(pathload, filename) for filename in file_list]
        with tqdm(total=len(file_list), ncols=150, desc='generating par') as pbar:
            for i_file in file_list:
                i_par_csv = DPM_reversiblemodel_fun_generatepar_fromPNAS_1(i_file)
                for i_par in i_par_csv:
                    i_par = DPM_reversiblemodel_fun_generatepar_fromPNAS_2(i_par)
                    PAR.append(i_par)
                pbar.update(1)

        df_PAR = pd.DataFrame.from_records(PAR)
        df_PAR.drop_duplicates(subset=df_PAR.columns.to_list()[1:], inplace=True)

        while len(df_PAR) > 0:
            i_df_PAR = df_PAR.iloc[:par_save_block_size, :]
            filename = os.path.join(pathsave_mut, f'par_reversiblemodel_para_' + block_format.format(block_count) +
                                    '_{:d}'.format(len(i_df_PAR)) + '.csv')
            i_df_PAR.to_csv(filename, index=False)
            block_count += 1
            df_PAR = df_PAR.iloc[par_save_block_size:, :]
    elif joint:
        if not os.path.exists(pathsave_parraw):
            os.makedirs(pathsave_parraw)
        file_list_raw = os.listdir(pathsave_parraw)
        file_list_raw = [filename for filename in file_list_raw if re.search(r'^rawpar_reversiblemodel', filename) is not None]

        file_list = os.listdir(pathsave_mut)
        file_list = [filename for filename in file_list if re.search(filename_pattern, filename) is not None]
        try:
            file_list.sort(key=lambda x: int(x.split('_')[-2]))
        except (ValueError, IndexError):
            pass
        file_list = [os.path.join(pathsave_mut, filename) for filename in file_list]

        α, μ = np.arange(α_min, α_max + α_step, α_step), 1
        # α = α_default
        par_αθκ = []
        max_Storr_σfull, min_Storr_σfull, max_rrtoS_σ0, min_rrtoS_σ0 = -1, np.inf, -1, np.inf
        for i_α in α:
            θ = np.arange(i_α, θ_max + θ_step, θ_step_nomut)
            θ = θ[θ <= θ_max]
            for i_θ in θ:
                i_par = (i_α, i_θ, κ)
                _, _, i_Storr_σfull, i_rrtoS_σ0, *_ = DPM_reversiblemodel_model_transitionrate(*i_par, μ)

                par_αθκ.append(i_par)
                max_Storr_σfull = max([i_Storr_σfull, max_Storr_σfull])
                min_Storr_σfull = min([i_Storr_σfull, min_Storr_σfull])

                max_rrtoS_σ0 = max([i_rrtoS_σ0, max_rrtoS_σ0])
                min_rrtoS_σ0 = min([i_rrtoS_σ0, min_rrtoS_σ0])

        val_max, val_min = max([max_Storr_σfull, max_rrtoS_σ0]), min([min_Storr_σfull, min_rrtoS_σ0])
        μ = np.arange(np.floor(-np.log10(val_max)) - 2, np.ceil(-np.log10(val_min)) + 2, 1)

        totalcombination = itertools.product(par_αθκ, μ)
        totalnum = len(par_αθκ) * len(μ)

        columns = ['α', 'θ', 'κ', 'μ', 's->rr σ0', 'rr->s σ0', 's->rr σ1', 'rr->s σ1']
        par_α_θ_μ, rate_rever = [], pd.DataFrame(columns=columns)

        with tqdm(total=totalnum, ncols=150, desc=f'generating par') as pbar:
            for i_par in totalcombination:
                rate_rever = DPM_reversiblemodel_fun_generatepar_fromPNAS_3(i_par, rate_rever)

        ind_select = [5, 7, 22, 24]
        par_α_θ_μ = [par_α_θ_μ[i] for i in ind_select]

        if not file_list_raw:
            PAR_pnas = []
            for i_file in file_list:
                i_par_csv = DPM_reversiblemodel_fun_generatepar_fromPNAS_1(i_file)
                PAR_pnas.extend(i_par_csv)

            paramID = []
            [paramID.append(i_par['paramID']) for i_par in PAR_pnas]
            assert len(set(paramID)) == len(paramID)
            paramID_str = [str(i_paramID) for i_paramID in paramID]
            paramID_count = dict(zip(paramID_str, [0]*len(paramID_str)))

            SR_00 = SR
            totalcombination = list(itertools.product(SR_00))
            totalnum = len(PAR_pnas) * len(SR)
            PAR = []
            with tqdm(total=totalnum, ncols=150, desc='generating par') as pbar:
                for i_par in PAR_pnas:
                    for i_SR_00 in totalcombination:
                        i_par_add = deepcopy(i_par)
                        DPM_reversiblemodel_fun_generatepar_fromPNAS_4(i_SR_00, i_par_add)
            totalnum = len(PAR)
            num_cores = 101  # multiprocessing.cpu_count()
            num_perchunk = int(np.ceil(totalnum/num_cores))
            num_inchunk = [num_perchunk] * num_cores
            num_inchunk[-1] = num_inchunk[-1] - (num_perchunk * num_cores - totalnum)

            totalcombination_chunked = []
            num = 0
            while num < totalnum:
                chunk_i = PAR[:num_perchunk]
                totalcombination_chunked.append(list(chunk_i))
                num = num + num_perchunk
                PAR = PAR[num_perchunk:]

            for i, i_chunk in enumerate(totalcombination_chunked):
                i_filename = os.path.join(pathsave_parraw, f'rawpar_reversiblemodel_{i}_' + '{:d}'.format(len(i_chunk)) + '.pckl')
                with bz2.BZ2File(i_filename, 'wb') as file:
                    pickle.dump(i_chunk, file)
        else:
            file_list_raw.sort(key=lambda x: int(x.split('_')[2]))

            # for i, i_file in enumerate(file_list_raw):
            #     DPM_reversiblemodel_fun_generatepar_fromPNAS_5(i, i_file)

            num_cores = multiprocessing.cpu_count()
            Parallel(n_jobs=num_cores)(delayed(DPM_reversiblemodel_fun_generatepar_fromPNAS_5)
                                       (i, i_file) for i, i_file in enumerate(file_list_raw))
    return


def DPM_reversiblemodel_fun_generatepar_fromfulltorever(pathload='./reversible_parcsv/',
                                                        pathsave_rever='./reversible_parcsv_rever',
                                                        filename_pattern='_para'):
    def DPM_reversiblemodel_fun_generatepar_fromfulltorever_1(i_file_):
        df = pd.read_csv(i_file_)
        par_csv = df.to_dict(orient='records')
        return par_csv

    def DPM_reversiblemodel_fun_generatepar_fromfulltorever_2(par_):
        par_['R1ratio'] = 0
        par_['R2ratio'] = 0
        par_['T1'] = 0
        par_['T2'] = 0
        par_['S_11_1'] = 1
        par_['S_11_2'] = 1
        return par_

    # %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%---begin---%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% #
    if not os.path.exists(pathsave_rever):
        os.makedirs(pathsave_rever)

    file_list = os.listdir(pathload)
    file_list = [filename for filename in file_list if re.search(filename_pattern, filename) is not None]
    paramID, block_count, PAR, PARset = 0, 1, [], set()

    try:
        file_list.sort(key=lambda x: int(x.split('_')[-2]))
    except (ValueError, IndexError):
        pass
    file_list = [os.path.join(pathload, filename) for filename in file_list]
    with tqdm(total=len(file_list), ncols=150, desc='generating par') as pbar:
        for i_file in file_list:
            i_par_csv = DPM_reversiblemodel_fun_generatepar_fromfulltorever_1(i_file)
            for i_par in i_par_csv:
                i_par = DPM_reversiblemodel_fun_generatepar_fromfulltorever_2(i_par)
                i_par_ = deepcopy(i_par)
                del i_par_['paramID']
                i_par_tuple = tuple(i_par_.values())

                if i_par_tuple not in PARset:
                    i_par['paramID'] = paramID
                    PARset.add(i_par_tuple)
                    PAR.append(i_par)
                    paramID += 1

                if len(PAR) == par_save_block_size:
                    i_df_par = pd.DataFrame.from_records(PAR)
                    filename = os.path.join(pathsave_rever, f'par_reversiblemodel_para_' + block_format.format(block_count) +
                                            '_{:d}'.format(len(i_df_par)) + '.csv')
                    i_df_par.to_csv(filename, index=False)
                    block_count += 1
                    PAR = []

            pbar.update(1)

    if len(PAR) > 0:
        i_df_par = pd.DataFrame.from_records(PAR)
        filename = os.path.join(pathsave_rever, f'par_reversiblemodel_para_' + block_format.format(block_count) +
                                '_{:d}'.format(len(i_df_par)) + '.csv')
        i_df_par.to_csv(filename, index=False)
        block_count += 1

    return


def DPM_reversiblemodel_fun_assignpar(par):
    α1, θ1, κ1, μ1, T1 = par[0]
    α2, θ2, κ2, μ2, T2 = par[1]
    g, S_01_1, S_01_2, S_00_1, S_00_2, S_11_1, S_11_2 = par[2]
    R1_ratio, R2_ratio = par[3]
    par = [α1, θ1, κ1, μ1, T1,
           α2, θ2, κ2, μ2, T2,
           g, S_01_1, S_01_2, S_00_1, S_00_2, S_11_1, S_11_2,
           R1_ratio, R2_ratio]
    return dict(zip(par_key, par))


def DPM_reversiblemodel_fun_generateparcsv(pathload='./reversible_par',
                                           pathsave='./reversible_parcsv',
                                           filename_pattern='.pckl'):
    if not os.path.exists(pathsave):
        os.makedirs(pathsave)

    file_list = os.listdir(pathload)
    file_list = [filename for filename in file_list if re.search(filename_pattern, filename) is not None]
    file_list.sort(key=lambda x: int(x.split('_')[2])*10+int(x.split('_')[3]))

    ind_chunk = list(set([int(i_file.split('_')[2]) for i_file in file_list]))
    ind_chunk.sort()
    paramID, block_count, PAR, PARset = 0, 1, [], set()
    with tqdm(total=len(file_list), ncols=150, desc='generating par') as pbar:
        for i_file in file_list:
            with bz2.BZ2File(os.path.join(pathload, i_file), 'rb') as f:
                par = pickle.load(f)
                for i, i_par in enumerate(par):
                    i_par['paramID'] = str(i_par['paramID'])
                    i_par_ = deepcopy(i_par)
                    del i_par_['paramID']
                    i_par_tuple = tuple(i_par_.values())
                    if i_par_tuple not in PARset:
                        i_SR_00_1 = i_par['S_00_1']/i_par['S_01_1']
                        i_SR_00_2 = i_par['S_00_2']/i_par['S_01_2']
                        i_α1, i_α2 = i_par['α1'], i_par['α2']
                        i_θ1, i_θ2 = i_par['θ1'], i_par['θ2']
                        i_κ1, i_κ2 = i_par['κ1'], i_par['κ2']
                        i_μ1, i_μ2 = i_par['μ1'], i_par['μ2']

                        assert isclose(i_SR_00_1, i_SR_00_2, rel_tol=1e-9)
                        assert i_α1 == i_α2
                        assert i_θ1 == i_θ2
                        assert i_κ1 == i_κ2
                        assert i_μ1 == i_μ2

                        PARset.add(i_par_tuple)
                        PAR.append(i_par)

                    assert i_par['S_00_1'] > i_par['S_11_1']
                    assert i_par['S_00_2'] > i_par['S_11_2']

                    if len(PAR) == par_save_block_size:
                        i_df_par = pd.DataFrame.from_records(PAR)
                        i_df_par = i_df_par.rename(columns=dict(zip(par_key, pd_colname)), errors='raise', inplace=False)
                        filename = os.path.join(pathsave, f'par_reversiblemodel_para_' + block_format.format(block_count) +
                                                '_{:d}'.format(len(i_df_par)) + '.csv')
                        i_df_par.to_csv(filename, index=False)
                        block_count += 1
                        PAR = []

                    assert len(PAR) < par_save_block_size

                pbar.update(1)
    if len(PAR) > 0:
        i_df_par = pd.DataFrame.from_records(PAR)
        i_df_par = i_df_par.rename(columns=dict(zip(par_key_all, pd_colname)), errors='raise', inplace=False)
        filename = os.path.join(pathsave, f'par_reversiblemodel_para_' + block_format.format(block_count) +
                                '_{:d}'.format(len(i_df_par)) + '.csv')
        i_df_par.to_csv(filename, index=False)
        block_count += 1
    return


def DPM_reversiblemodel_fun_run_csvfolder(pathload='./reversible_parcsv',
                                          pathsave='./reversible_csvresult',
                                          filename_pattern='par_reversiblemodel_para',
                                          Strategy_name=strategy_default[:-7],
                                          save_filename_param=True,
                                          save_filename_stopt=True,
                                          save_filename_dosage=True,
                                          mutonly=False,
                                          reveronly=False,
                                          use_parallel=False):
    if reveronly:
        [Strategy_name.remove(i_strategy) for i_strategy in ['strategy2', 'strategy2 cycle'] if i_strategy in Strategy_name]
    if not os.path.exists(pathsave):
        os.makedirs(pathsave)
    file_list_all = os.listdir(pathload)
    file_list = [filename for filename in file_list_all if re.search(filename_pattern, filename) is not None]
    file_list.sort(key=lambda x: int(x.split('_')[3]))
    file_list = [os.path.join(pathload, filename) for filename in file_list]

    if use_parallel and len(file_list) > 1:
        num_cores = multiprocessing.cpu_count()
        Parallel(n_jobs=num_cores)(delayed(DPM_reversiblemodel_fun_run_csv)
                                   (filename=i_filename,
                                    Strategy_name=Strategy_name,
                                    save_filename_param=save_filename_param,
                                    save_filename_stopt=save_filename_stopt,
                                    save_filename_dosage=save_filename_dosage,
                                    mutonly=mutonly,
                                    pathsave=pathsave)for i, i_filename in enumerate(file_list))
        print('Finished.')
    else:
        count_file = 1
        for i, i_filename in enumerate(file_list):
            print('Processing the ' + str(count_file) + 'th' + ' file of ' + str(len(file_list)) + ' total files...')
            DPM_reversiblemodel_fun_run_csv(filename=i_filename,
                                            Strategy_name=Strategy_name,
                                            save_filename_param=save_filename_param,
                                            save_filename_stopt=save_filename_stopt,
                                            save_filename_dosage=save_filename_dosage,
                                            mutonly=mutonly,
                                            pathsave=pathsave)
            count_file += 1
        print('Finished.')


def DPM_reversiblemodel_fun_run_csv(filename=None,
                                    Strategy_name=None,
                                    save_filename_param=None,
                                    save_filename_stopt=None,
                                    save_filename_dosage=None,
                                    mutonly=False,
                                    pathsave=None):
    df = pd.read_csv(filename)
    df.rename(columns=dict(zip(pd_colname, par_key)), errors='raise', inplace=True)
    par_csv = df.to_dict(orient='records')

    timenow = datetime.now()
    timenow = timenow.strftime('%Y%m%d')
    filename = filename.replace('.csv', '')
    if sys.platform in ['darwin', 'linux']:
        filename = filename[filename.rfind('/') + 1:]
    else:
        filename = filename[filename.rfind('\\') + 1:]

    filename_param = filename + '_result_para_' + timenow + '.csv' if save_filename_param else ''
    filename_stopt = filename + '_result_stopt_' + timenow + '.csv' if save_filename_stopt else ''
    filename_dosage = filename + '_result_dosage_' + timenow + '.csv' if save_filename_dosage else ''
    filename_Nflow = filename + '_result_Nflow_' + timenow + '.pckl'
    filename_pop = filename + '_result_pop_' + timenow + '.pckl'

    pathsave = os.path.abspath(pathsave)
    if not os.path.exists(pathsave):
        os.makedirs(pathsave)
    filename_param = os.path.join(pathsave, filename_param) if save_filename_param else ''
    filename_stopt = os.path.join(pathsave, filename_stopt) if save_filename_stopt else ''
    filename_dosage = os.path.join(pathsave, filename_dosage) if save_filename_dosage else ''
    filename_Nflow = os.path.join(pathsave, filename_Nflow)

    if not mutonly:
        save_filename_param and DPM_reversiblemodel_miscellaneous_writehead(filename_param, Heading_param_full)
    else:
        save_filename_param and DPM_reversiblemodel_miscellaneous_writehead(filename_param, Heading_param_mutonly)
    save_filename_stopt and DPM_reversiblemodel_miscellaneous_writehead(filename_stopt, Heading_stopt)
    save_filename_dosage and DPM_reversiblemodel_miscellaneous_writehead(filename_dosage, Heading_dosage)

    N_flowcol = dict(zip(Strategy_name, [np.zeros((len(cellflow), duration_5year))]*len(Strategy_name)))
    N_flowcolextend = dict(zip(Strategy_name, [np.zeros((len(cellflow), duration_5year))] * len(Strategy_name)))
    N_flow_num = dict(zip(Strategy_name, [np.zeros(duration_5year)]*len(Strategy_name)))
    with tqdm(total=len(par_csv), ncols=150, desc='Runing simulation for the parameters from .csv file input') as pbar:
        for i in range(len(par_csv)):
            i_par = par_csv[i]
            i_result = DPM_reversiblemodel_fun_run_simulation(i_par, Strategy_name, mutonly=mutonly)
            '''calculate Fig. 5'''
            for j, j_strategy in enumerate(strategy_default):
                if j_strategy in Strategy_name:
                    result_j = i_result[j]
                    if j_strategy == 'mono drug1':
                        d_j = np.tile(np.array([σ1_full], dtype=float).T, (1, result_j.shape[1] - 1))
                        N_j = result_j
                    elif j_strategy == 'mono drug2':
                        d_j = np.tile(np.array([σ2_full], dtype=float).T, (1, result_j.shape[1] - 1))
                        N_j = result_j
                    elif j_strategy == 'mono half':
                        d_j = np.tile(np.array([σ_half], dtype=float).T, (1, result_j.shape[1] - 1))
                        N_j = result_j
                    elif j_strategy == 'best cycle':
                        d_j = result_j['d'][0]
                        N_j = result_j['N']
                        d_j = DPM_reversiblemodel_miscellaneous_cycletreat2dose(d_j, math.ceil(N_j.shape[1]/Stepsize))
                        d_j = d_j[:, :N_j.shape[1]-1]
                    else:
                        d_j = result_j['d']
                        N_j = result_j['N']

                    _, i_strategy_flowcol, i_strategy_flowcolextend = DPM_reversiblemodel_fun_model_flow(i_par, N_j, d_j)

                    N_flowcol[j_strategy] = N_flowcol[j_strategy] + i_strategy_flowcol
                    N_flowcolextend[j_strategy] = N_flowcolextend[j_strategy] + i_strategy_flowcolextend
                    N_flow_num[j_strategy] = N_flow_num[j_strategy] + (np.amax(i_strategy_flowcol, axis=0) != 0)*1

            DPM_reversiblemodel_readsave_save_csv(i_par,
                                                  i_result,
                                                  filename_param,
                                                  filename_stopt,
                                                  filename_dosage,
                                                  filename_pop)
            pbar.update(1)

    with bz2.BZ2File(filename_Nflow, 'wb') as f:
        pickle.dump((N_flowcol, N_flowcolextend, N_flow_num), f)
    return


def DPM_reversiblemodel_fun_run_simulation_cycle(par):
    mutation = False
    n_ratio = DPM_reversiblemodel_fun_steadystate(par, mutation)[2]
    x0 = DPM_reversiblemodel_fun_initx(n_ratio, par)

    mutation, LSsim, t = True, True, np.arange(0, duration_5year + Simtimestep, Simtimestep)
    σ1_duration = np.arange(daysperweek, max_singledrugduration * daysperweek + daysperweek, daysperweek)
    σ2_duration = np.arange(daysperweek, max_singledrugduration * daysperweek + daysperweek, daysperweek)
    N_cycle_σ1first = -np.ones((len(σ1_duration), len(σ2_duration)))
    N_cycle_σ2first = -np.ones((len(σ1_duration), len(σ2_duration)))

    for i, i_σ1_duration in enumerate(σ1_duration):
        print(i)
        for j, j_σ2_duration in enumerate(σ2_duration):
            i_sim = DPM_reversible_run_cycle(i_σ1_duration, j_σ2_duration, t, x0, mutation, par)
            N_cycle_σ1first[i, j] = i_sim[0]
            N_cycle_σ2first[i, j] = i_sim[1]
    return N_cycle_σ1first, N_cycle_σ2first


def DPM_reversible_run_cycle(σ1_duration, σ2_duration, t, x0, mutation, pars):
    cycle_duration = σ1_duration + σ2_duration
    σ1_ = np.row_stack((np.ones((1, σ1_duration)), np.zeros((1, σ1_duration))))
    σ2_ = np.row_stack((np.zeros((1, σ2_duration)), np.ones((1, σ2_duration))))
    cycle_σ1first = np.concatenate((σ1_, σ2_), axis=1)
    cycle_σ2first = np.concatenate((σ2_, σ1_), axis=1)
    num_cycle = math.ceil((len(t) - 1) / cycle_duration)

    cycle_σ1first = np.tile(cycle_σ1first, num_cycle)[:, 0:len(t) - 1]
    cycle_σ2first = np.tile(cycle_σ2first, num_cycle)[:, 0:len(t) - 1]

    _, N_total_σ1first, *_ = DPM_reversiblemodel_strategy_simbydose(x0, cycle_σ1first, pars, mutation)
    _, N_total_σ2first, *_ = DPM_reversiblemodel_strategy_simbydose(x0, cycle_σ2first, pars, mutation)

    return N_total_σ1first[-1], N_total_σ2first[-1]


def DPM_reversiblemodel_fun_run_simulation(par, strategy_name, mutonly=False):
    strategy_name = [strategy_name] if type(strategy_name) is not list else strategy_name
    strategy_name = [i.lower() for i in strategy_name]
    mono_drug1, mono_drug2, mono_half, best_cycle, strategy0, strategy1, strategy1_cycle, strategy2, strategy2_cycle, strategy_random = \
        [None] * len(strategy_default)

    survi_time_nocycle, survi_time = [], []
    cure_nocycle, cure = [], []
    N_nocycle, N, d_nocycle, d = [], [], [], []

    mutation = False
    n_ratio = DPM_reversiblemodel_fun_steadystate(par, mutation)[2]
    x0 = DPM_reversiblemodel_fun_initx(n_ratio, par)
    if mutonly:
        assert np.all(np.array(x0)[ind_rever] == 0)

    mutation, LSsim, t = True, True, np.arange(0, duration_5year + Simtimestep, Simtimestep)
    tdiff = np.append(0, np.diff(t))
    ''' σ1 full dose '''
    if 'mono drug1' in strategy_name:
        argsODE = (σ1_full, mutation, par)
        mono_drug1, *_ = DPM_reversiblemodel_model_simexpm(argsODE, tdiff, x0)

        survi_time_mono_drug1 = DPM_reversiblemodel_stime(mono_drug1)
        d_mono_drug1 = np.tile(np.array([σ1_full], dtype=float).T, (1, mono_drug1.shape[1]-1))
        flag_mono_drug1_cure = DPM_reversiblemodel_cure(mono_drug1)

        if mutonly:
            assert np.all(mono_drug1[ind_rever, :] == 0)

        survi_time_nocycle.append(survi_time_mono_drug1), survi_time.append(survi_time_mono_drug1)
        cure_nocycle.append(flag_mono_drug1_cure), cure.append(flag_mono_drug1_cure)
        d_nocycle.append(d_mono_drug1), d.append(d_mono_drug1)
        N_nocycle.append(mono_drug1), N.append(mono_drug1)

    '''σ2 full dose'''
    if 'mono drug2' in strategy_name:
        argsODE = (σ2_full, mutation, par)
        mono_drug2, *_ = DPM_reversiblemodel_model_simexpm(argsODE, tdiff, x0)

        survi_time_mono_drug2 = DPM_reversiblemodel_stime(mono_drug2)
        d_mono_drug2 = np.tile(np.array([σ2_full], dtype=float).T, (1, mono_drug2.shape[1]-1))
        flag_mono_drug2_cure = DPM_reversiblemodel_cure(mono_drug2)

        if mutonly:
            assert np.all(mono_drug2[ind_rever, :] == 0)

        survi_time_nocycle.append(survi_time_mono_drug2), survi_time.append(survi_time_mono_drug2)
        cure_nocycle.append(flag_mono_drug2_cure), cure.append(flag_mono_drug2_cure)
        d_nocycle.append(d_mono_drug2), d.append(d_mono_drug2)
        N_nocycle.append(mono_drug2), N.append(mono_drug2)

    '''σ1=0.5 and σ2=0.5'''
    if 'mono half' in strategy_name:
        argsODE = (σ_half, mutation, par)
        mono_half, *_ = DPM_reversiblemodel_model_simexpm(argsODE, tdiff, x0)

        survi_time_mono_half = DPM_reversiblemodel_stime(mono_half)
        d_σ_half = np.tile(np.array([σ_half], dtype=float).T, (1, mono_half.shape[1]-1))
        flag_mono_half_cure = DPM_reversiblemodel_cure(mono_half)

        if mutonly:
            assert np.all(mono_half[ind_rever, :] == 0)

        survi_time_nocycle.append(survi_time_mono_half), survi_time.append(survi_time_mono_half)
        cure_nocycle.append(flag_mono_half_cure), cure.append(flag_mono_half_cure)
        d_nocycle.append(d_σ_half), d.append(d_σ_half)
        N_nocycle.append(mono_half), N.append(mono_half)

    duration = duration_5year
    '''best cycle'''
    if 'best cycle' in strategy_name:
        N_bestcycle, d_bestcycle = DPM_reversiblemodel_strategy_bestcycle(x0, par, duration, mutation)
        best_cycle = {'N': N_bestcycle, 'd': d_bestcycle}
        if mutonly:
            assert np.all(N_bestcycle[ind_rever, :] == 0)

        survi_time_best_cycle = DPM_reversiblemodel_stime(N_bestcycle)
        flag_best_cycle_cure = DPM_reversiblemodel_cure(N_bestcycle)

        survi_time.append(survi_time_best_cycle)
        cure.append(flag_best_cycle_cure)
        d.append(d_bestcycle)
        N.append(N_bestcycle)

    '''strategy0'''
    if 'strategy0' in strategy_name:
        N_strategy0, d_strategy0 = DPM_reversiblemodel_strategy_0(x0, par, duration, mutation)
        strategy0 = {'N': N_strategy0, 'd': d_strategy0}
        if mutonly:
            assert np.all(N_strategy0[ind_rever, :] == 0)

        survi_time_strategy0 = DPM_reversiblemodel_stime(N_strategy0)
        flag_strategy0_cure = DPM_reversiblemodel_cure(N_strategy0)

        survi_time_nocycle.append(survi_time_strategy0), survi_time.append(survi_time_strategy0)
        cure_nocycle.append(flag_strategy0_cure), cure.append(flag_strategy0_cure)
        d_nocycle.append(d_strategy0), d.append(d_strategy0)
        N_nocycle.append(N_strategy0), N.append(N_strategy0)

    '''strategy1'''
    if 'strategy1' in strategy_name:
        N_strategy1, d_strategy1 = DPM_reversiblemodel_strategy_1(x0, par, duration, mutation)
        strategy1 = {'N': N_strategy1, 'd': d_strategy1}
        if mutonly:
            assert np.all(N_strategy1[ind_rever, :] == 0)

        survi_time_strategy1 = DPM_reversiblemodel_stime(N_strategy1)
        flag_strategy1_cure = DPM_reversiblemodel_cure(N_strategy1)

        survi_time_nocycle.append(survi_time_strategy1), survi_time.append(survi_time_strategy1)
        cure_nocycle.append(flag_strategy1_cure), cure.append(flag_strategy1_cure)
        d_nocycle.append(d_strategy1), d.append(d_strategy1)
        N_nocycle.append(N_strategy1), N.append(N_strategy1)

    '''strategy1 cycle'''
    if 'strategy1 cycle' in strategy_name:
        N_strategy1_cycle, d_strategy1_cycle = DPM_reversiblemodel_strategy_1(x0, par, duration, mutation, use_cycle=True)
        strategy1_cycle = {'N': N_strategy1_cycle, 'd': d_strategy1_cycle}
        if mutonly:
            assert np.all(N_strategy1_cycle[ind_rever, :] == 0)

        survi_time_strategy1_cycle = DPM_reversiblemodel_stime(N_strategy1_cycle)
        flag_strategy1_cycle_cure = DPM_reversiblemodel_cure(N_strategy1_cycle)

        survi_time.append(survi_time_strategy1_cycle)
        cure.append(flag_strategy1_cycle_cure)
        d.append(d_strategy1_cycle)
        N.append(N_strategy1_cycle)

    '''strategy2'''
    if 'strategy2' in strategy_name:
        N_strategy2, d_strategy2 = DPM_reversiblemodel_strategy_2(x0, par, duration, mutation)
        strategy2 = {'N': N_strategy2, 'd': d_strategy2}
        if mutonly:
            assert np.all(N_strategy2[ind_rever, :] == 0)

        survi_time_strategy2 = DPM_reversiblemodel_stime(N_strategy2)
        flag_strategy2_cure = DPM_reversiblemodel_cure(N_strategy2)

        survi_time_nocycle.append(survi_time_strategy2), survi_time.append(survi_time_strategy2)
        cure_nocycle.append(flag_strategy2_cure), cure.append(flag_strategy2_cure)
        d_nocycle.append(d_strategy2), d.append(d_strategy2)
        N_nocycle.append(N_strategy2), N.append(N_strategy2)

    '''strategy2 cycle'''
    if 'strategy2 cycle' in strategy_name:
        N_strategy2_cycle, d_strategy2_cycle = DPM_reversiblemodel_strategy_2(x0, par, duration, mutation, use_cycle=True)
        strategy2_cycle = {'N': N_strategy2_cycle, 'd': d_strategy2_cycle}
        if mutonly:
            assert np.all(N_strategy2_cycle[ind_rever, :] == 0)

        survi_time_strategy2_cycle = DPM_reversiblemodel_stime(N_strategy2_cycle)
        flag_strategy2_cycle_cure = DPM_reversiblemodel_cure(N_strategy2_cycle)

        survi_time.append(survi_time_strategy2_cycle)
        cure.append(flag_strategy2_cycle_cure)
        d.append(d_strategy2_cycle)
        N.append(N_strategy2_cycle)

    '''random'''
    if 'random' in strategy_name:
        N_random, d_random = DPM_reversiblemodel_strategy_random(x0, par, duration, mutation)
        strategy_random = {'N': N_random, 'd': d_random}

    return mono_drug1, \
        mono_drug2, \
        mono_half, \
        best_cycle, \
        strategy0, \
        strategy1, \
        strategy1_cycle, \
        strategy2, \
        strategy2_cycle, \
        strategy_random


def DPM_reversiblemodel_fun_label(stoptime_df, strategy_name):
    stoptime_df['label'] = [-1]*stoptime_df.shape[0]
    with tqdm(total=stoptime_df.shape[0], ncols=150, desc='Calculating label...') as pbar:
        for i, i_row in stoptime_df.iterrows():
            i_row = i_row[strategy_name]
            if 'mono half' in i_row.index[i_row == i_row.max()]:
                stoptime_df['label'].iloc[i] = 0
            elif 'best cycle' in i_row.index[i_row == i_row.max()]:
                stoptime_df['label'].iloc[i] = 1
            elif 'strategy1' in i_row.index[i_row == i_row.max()]:
                stoptime_df['label'].iloc[i] = 2
            elif 'strategy2' in i_row.index[i_row == i_row.max()]:
                stoptime_df['label'].iloc[i] = 3
            elif 'strategy1 cycle' in i_row.index[i_row == i_row.max()]:
                stoptime_df['label'].iloc[i] = 4
            elif 'strategy2 cycle' in i_row.index[i_row == i_row.max()]:
                stoptime_df['label'].iloc[i] = 5
            else:
                raise ValueError
            pbar.update(1)
    return stoptime_df['label']


def DPM_reversiblemodel_fun_feature(para_df):
    for i_name in ['rate_Storr1_σhalf', 'rate_rr1toS_σhalf', 'rate_Storr2_σhalf', 'rate_rr2toS_σhalf']:
        para_df[i_name] = [-1] * para_df.shape[0]

    para_df.rename(columns={'rate_S1torr1_sigma0': 'rate_Storr1_sigma0', 'rate_rr1toS1_sigma0': 'rate_rr1toS_sigma0',
                            'rate_S2torr2_sigma0': 'rate_Storr2_sigma0', 'rate_rr2toS2_sigma0': 'rate_rr2toS_sigma0',
                            'rate_S1torr1_sigmafull': 'rate_Storr1_sigmafull', 'rate_rr1toS1_sigmafull': 'rate_rr1toS_sigmafull',
                            'rate_S2torr2_sigmafull': 'rate_Storr2_sigmafull', 'rate_rr2toS2_sigmafull': 'rate_rr2toS_sigmafull'},
                   inplace=True)

    with tqdm(total=para_df.shape[0], ncols=150, desc='Calculating feature...') as pbar:
        for i, i_row in para_df.iterrows():
            _, i_rate_Storr1_σhalf, _, _, i_rate_rr1toS_σhalf, _ = (
                DPM_reversiblemodel_model_transitionrate(i_row['alpha1'], i_row['theta1'], i_row['kappa1'], i_row['mu1']))

            _, i_rate_Storr2_σhalf, _, _, i_rate_rr2toS_σhalf, _ = (
                DPM_reversiblemodel_model_transitionrate(i_row['alpha2'], i_row['theta2'], i_row['kappa2'], i_row['mu2']))

            para_df['rate_Storr1_σhalf'].iloc[i] = i_rate_Storr1_σhalf
            para_df['rate_rr1toS_σhalf'].iloc[i] = i_rate_rr1toS_σhalf
            para_df['rate_Storr2_σhalf'].iloc[i] = i_rate_Storr2_σhalf
            para_df['rate_rr2toS_σhalf'].iloc[i] = i_rate_rr2toS_σhalf
            pbar.update(1)

    feature_name = ['alpha1', 'theta1', 'mu1', 'alpha2', 'theta2', 'mu2', 'g11', 'g10', 'g01', 'g00', 'S_01_1', 'S_01_2', 'S_00_1',
                    'S_00_2', 'S_11_1', 'S_11_2',  'T1', 'T2', 'R1ratio', 'R2ratio',
                    'rate_Storr1_sigma0', 'rate_rr1toS_sigma0', 'rate_Storr2_sigma0', 'rate_rr2toS_sigma0',
                    'rate_Storr1_σhalf', 'rate_rr1toS_σhalf', 'rate_Storr2_σhalf', 'rate_rr2toS_σhalf',
                    'rate_Storr1_sigmafull', 'rate_rr1toS_sigmafull', 'rate_Storr2_sigmafull', 'rate_rr2toS_sigmafull']
    para_df = para_df[feature_name]
    para_df_norm = para_df.copy()
    for i_feature_name in feature_name:
        para_df_norm[i_feature_name] = ((para_df_norm[i_feature_name]-para_df_norm[i_feature_name].min())
                                        / (para_df_norm[i_feature_name].max()-para_df_norm[i_feature_name].min()))
        assert para_df_norm[i_feature_name].min() == 0
        assert para_df_norm[i_feature_name].max() == 1
    return {'unnorm': para_df, 'norm': para_df_norm}


def DPM_reversiblemodel_fun_ayalysis(pathload='./reversible_csvresult/',
                                     pathload_mutonly='./reversible_csvresult_mut',
                                     pathload_reveronly='./reversible_csvresult_rever',
                                     strategy_name=strategy_default[:-7],
                                     mutonly=False,
                                     reveronly=False):

    def DPM_reversiblemodel_fun_analysis_1(i_stoptime_):
        i_filename_pattern = i_stoptime_.split('_stopt', 1)[0]
        i_file_para, i_file_dosage, i_file_pop, i_file_Nflow \
            = i_filename_pattern + '_para', i_filename_pattern + '_dosage', i_filename_pattern + '_pop', i_filename_pattern + '_Nflow'

        i_file_stoptime_ = os.path.join(pathload, i_stoptime_)
        i_stoptime = DPM_reversiblemodel_readsave_stoptime_csv(i_file_stoptime_, strategy_name)

        i_file_para = [i_file for i_file in file_para if i_file_para in i_file]
        i_file_dosage = [i_file for i_file in file_dosage if i_file_dosage in i_file]
        i_file_Nflow = [i_file for i_file in file_Nflow if i_file_Nflow in i_file]

        i_file_para = os.path.join(pathload, i_file_para[0])
        i_file_dosage = os.path.join(pathload, i_file_dosage[0])
        if (not mutonly) and (not reveronly):
            i_file_Nflow = os.path.join(pathload, i_file_Nflow[0])
        # i_file_pop = os.path.join(pathload, i_file_pop[0])

        i_para = pd.read_csv(i_file_para, usecols=['paramID'] + pd_colname).to_dict('list')
        i_dosage, i_dosage_idx = DPM_reversiblemodel_readsave_dosage_csv(i_file_dosage, strategy_name)

        i_num = len(i_para[list(i_para.keys())[0]])
        result = np.zeros((i_num,), dtype=float)
        Storr1_σ0, Storr2_σ0, rr1toS_σ0, rr2toS_σ0, Storr1_σfull, Storr2_σfull, rr1toS_σfull, rr2toS_σfull = \
            [deepcopy(result) for _ in range(8)]
        S_σ1_full, S_σ2_full, R1_σ1_full, R2_σ2_full, IR1_σ1_full, IR2_σ2_full = [deepcopy(result) for _ in range(6)]

        for i_ in range(i_num):
            Storr1_σ0[i_], _, Storr1_σfull[i_], rr1toS_σ0[i_], _, rr1toS_σfull[i_] = \
                DPM_reversiblemodel_model_transitionrate(i_para['alpha1'][i_],
                                                         i_para['theta1'][i_],
                                                         i_para['kappa1'][i_],
                                                         i_para['mu1'][i_])

            Storr2_σ0[i_], _, Storr2_σfull[i_], rr2toS_σ0[i_], _, rr2toS_σfull[i_] = \
                DPM_reversiblemodel_model_transitionrate(i_para['alpha2'][i_],
                                                         i_para['theta2'][i_],
                                                         i_para['kappa2'][i_],
                                                         i_para['mu2'][i_])

            i_rate = DPM_reversiblemodel_model_prorate(i_para['g'][i_],
                                                       i_para['S_01_1'][i_],
                                                       i_para['S_01_2'][i_],
                                                       i_para['S_00_1'][i_],
                                                       i_para['S_00_2'][i_],
                                                       i_para['S_11_1'][i_],
                                                       i_para['S_11_2'][i_])
            S_σ1_full[i_] = i_rate['S σ1 full']
            S_σ2_full[i_] = i_rate['S σ2 full']
            R1_σ1_full[i_] = i_rate['R1 σ1 full']
            R2_σ2_full[i_] = i_rate['R2 σ2 full']
            IR1_σ1_full[i_] = i_rate['IR1 σ1 full']
            IR2_σ2_full[i_] = i_rate['IR2 σ2 full']

        i_para['Storr1_σ0'] = Storr1_σ0
        i_para['rr1toS_σ0'] = rr1toS_σ0
        i_para['Storr2_σ0'] = Storr2_σ0
        i_para['rr2toS_σ0'] = rr2toS_σ0
        i_para['Storr1_σfull'] = rr1toS_σfull

        i_para['Storr1_σ0'] = Storr1_σ0
        i_para['rr1toS_σ0'] = rr1toS_σ0
        i_para['Storr2_σ0'] = Storr2_σ0
        i_para['rr2toS_σ0'] = rr2toS_σ0
        i_para['Storr1_σfull'] = Storr1_σfull
        i_para['rr1toS_σfull'] = rr1toS_σfull
        i_para['Storr2_σfull'] = Storr2_σfull
        i_para['rr2toS_σfull'] = rr2toS_σfull

        i_para['S σ1 full'] = S_σ1_full
        i_para['S σ2 full'] = S_σ2_full
        i_para['R1 σ1 full'] = R1_σ1_full
        i_para['R2 σ2 full'] = R2_σ2_full
        i_para['IR1 σ1 full'] = IR1_σ1_full
        i_para['IR2 σ2 full'] = IR2_σ2_full

        if mutonly:
            assert np.all(i_para['Storr1_σ0'] == 0) and \
                   np.all(i_para['rr1toS_σ0'] == 0) and \
                   np.all(i_para['Storr2_σ0'] == 0) and \
                   np.all(i_para['rr2toS_σ0'] == 0)

        for i_ in range(len(i_para[list(i_para.keys())[0]])):
            i_par = dict(zip(i_para.keys(), [i_para[i_key][i_] for i_key, _ in i_para.items()]))

            flag1_mutonly = i_par['S_11_1']/i_par['S_01_1'] < 0.1 and i_par['S_11_2']/i_par['S_01_2'] < 0.1
            flag2_mutonly = i_par['R1ratio'] > 0 and i_par['R2ratio'] > 0
            flag3_mutonly = i_par['S_01_1']/i_par['g'] > 0.9 and i_par['S_01_2']/i_par['g'] > 0.9
            flag4_mutonly = i_par['R1ratio'] < 0.8 and i_par['R2ratio'] < 0.8
            flag5_mutonly = i_par['g'] < 0.29
            flag_mutonly = flag1_mutonly and flag2_mutonly and flag3_mutonly and flag4_mutonly and flag5_mutonly

            if mutonly:
                flag = flag_mutonly
            else:
                flag = True
            if flag:
                for i_key in para.keys():
                    para[i_key].extend([i_par[i_key]])
                for i_key in stoptime.keys():
                    stoptime[i_key].extend([i_stoptime[i_key][i_]])
                for i_key in dosage.keys():
                    if i_key == 'mono drug1':
                        i_dose_mono_drug1 = DPM_reversiblemodel_miscellaneous_treat_mono_or_half(i_stoptime[i_key][i_], 'sigma1 full')
                        dosage[i_key].extend([i_dose_mono_drug1[:-1]])
                    elif i_key == 'mono drug2':
                        i_dose_mono_drug2 = DPM_reversiblemodel_miscellaneous_treat_mono_or_half(i_stoptime[i_key][i_], 'sigma2 full')
                        dosage[i_key].extend([i_dose_mono_drug2[:-1]])
                    elif i_key == 'mono half':
                        i_dose_mono_half = DPM_reversiblemodel_miscellaneous_treat_mono_or_half(i_stoptime[i_key][i_], 'half')
                        dosage[i_key].extend([i_dose_mono_half[:-1]])
                    else:
                        dosage[i_key].extend([i_dosage[i_key][i_]])
                assert i_par['paramID'] == i_dosage['paramID'][i_] == i_stoptime['paramID'][i_]
                for i_key in bool_dosage.keys():
                    if i_key == 'paramID':
                        bool_dosage[i_key].extend([i_dosage_idx[i_key][i_]])
                    else:
                        if 'strategy0' in strategy_name:
                            bool_dosage[i_key].extend([i_dosage_idx[i_key][i_]])
                        else:
                            bool_dosage[i_key].extend([None])

        if i_file_Nflow and (not mutonly) and (not reveronly):
            with bz2.BZ2File(i_file_Nflow, 'rb') as i_f:
                i_Nflowcol, i_Nflowcolextend, i_Nflownum = pickle.load(i_f)
            for i_strategy in strategy_name:
                Nflow['valcol'][i_strategy] = Nflow['valcol'][i_strategy] + i_Nflowcol[i_strategy]
                Nflow['valcolextend'][i_strategy] = Nflow['valcolextend'][i_strategy] + i_Nflowcolextend[i_strategy]
                Nflow['num'][i_strategy] = Nflow['num'][i_strategy] + i_Nflownum[i_strategy]

        print(round(len(para['paramID'])/10000, 2))
        return

    def DPM_reversiblemodel_fun_analysis_2(info_,
                                           hz_,
                                           strategy_name_,
                                           filename_stoptimetable_,
                                           filename_hzratio_):
        num_patient = len(info_['stoptime_df'])
        data = list()
        data.append(['Median survival,\nday'] + list(info_['survival_median']))
        data.append(['Mean survival,\nday'] + list(info_['survival_mean']))
        data.append(['No. of cases \nsurvival at 5y'] + list(info_['num_survival_5y'].values()))
        data.append(['Survival at\n5y, %'] + list(info_['percentage_survival_5y'].values()))
        data.append(['No. of cases strategy\nnumerically better\nthan all others'] + list(info_['num_beststrategy'].values()))
        data.append(['Cases strategy\nnumerically better\nthan all others, %'] + list(info_['percentage_beststrategy'].values()))
        data.append(['No. of cases strategy\nsignificantly better\nthan all others'] + list(info_['num_sigbeststrategy'].values()))
        data.append(['Cases strategy\nsignificantly better\nthan all others, %'] + list(info_['percentage_sigbeststrategy'].values()))
        col_names = [f'No. of patients\n: {num_patient}']
        if len(strategy_name_) == 9:
            colname = ['mono\n drug1',
                       'mono\ndrug2',
                       'mono\nhalf',
                       'best\ncycle',
                       'strategy0',
                       'strategy1',
                       'strategy1\ncycle',
                       'strategy2',
                       'strategy2\ncycle']
        elif len(strategy_name_) == 7:
            colname = ['mono\n drug1',
                       'mono\ndrug2',
                       'mono\nhalf',
                       'best\ncycle',
                       'strategy0',
                       'strategy1',
                       'strategy1\ncycle']
        else:
            colname = []

        col_names.extend(colname)
        print(tabulate(data, headers=col_names, tablefmt='rst', numalign='center'))
        data_df = pd.DataFrame(data, columns=col_names)
        DPM_reversiblemodel_plot_df2pdf(data_df, filename_stoptimetable_, (['white'], ['lightgray']))

        betterstrategy = [None] * len(strategy_name_)
        sigbetterstrategy = [None] * len(strategy_name_)
        hzratio = [None] * len(strategy_name_)
        improve_beststrategy, improve_sigbeststrategy = (pd.DataFrame(columns=['improve days', 'strategy']),
                                                         pd.DataFrame(columns=['improve days', 'strategy']))
        for ii, i_strategy in enumerate(strategy_name_):
            i_num_betterstrategy = list(info_['num_betterstrategy'][i_strategy].values())
            i_num_sigbetterstrategy = list(info_['num_sigbetterstrategy'][i_strategy].values())
            i_hzratio_strategy = hz_['hz'][i_strategy]
            i_hzratio_strategy[i_strategy] = 'N.A.'

            i_percentage_betterstrategy = list(info_['percentage_betterstrategy'][i_strategy].values())
            i_percentage_sigbetterstrategy = list(info_['percentage_sigbetterstrategy'][i_strategy].values())

            i_betterstrategy, i_sigbetterstrategy, i_hzratio = [colname[ii]], [colname[ii]], [colname[ii]]
            i_betterstrategy.extend([str(i_num_betterstrategy[j]) + '\n' + str(i_percentage_betterstrategy[j])
                                     for j in range(len(strategy_name_))])
            i_sigbetterstrategy.extend([str(i_num_sigbetterstrategy[j]) + '\n' + str(i_percentage_sigbetterstrategy[j])
                                        for j in range(len(strategy_name_))])
            i_hzratio.extend([str(round(i_hzratio_strategy[j], 3)) if type(i_hzratio_strategy[j]) != str else i_hzratio_strategy[j]
                              for j in strategy_name_])

            betterstrategy[ii], sigbetterstrategy[ii], hzratio[ii] = i_betterstrategy, i_sigbetterstrategy, i_hzratio

            i_improve_beststrategy = pd.DataFrame({'improve days': info_['improve_beststrategy'][i_strategy],
                                                   'strategy': [colname[ii]] * info_['improve_beststrategy'][i_strategy].shape[0]})
            i_improve_sigbeststrategy = pd.DataFrame({'improve days': info_['improve_sigbeststrategy'][i_strategy],
                                                      'strategy': [colname[ii]] * info_['improve_sigbeststrategy'][i_strategy].shape[0]})

            improve_beststrategy = pd.concat([improve_beststrategy, i_improve_beststrategy], ignore_index=True)
            improve_sigbeststrategy = pd.concat([improve_sigbeststrategy, i_improve_sigbeststrategy], ignore_index=True)

        hzratio_df = pd.DataFrame(hzratio, columns=col_names)
        DPM_reversiblemodel_plot_df2pdf(hzratio_df, filename_hzratio_, (['white'], ['lightgray']))

        improve_beststrategy[['improve days']] = improve_beststrategy[['improve days']].apply(pd.to_numeric)
        improve_sigbeststrategy[['improve days']] = improve_sigbeststrategy[['improve days']].apply(pd.to_numeric)

        return

    def DPM_reversiblemodel_fun_analysis_3(stoptime_, strategy_name_, filename_km_):
        km = {**{i_strategy: None for i_strategy in strategy_name_}}
        for i_strategy in strategy_name:
            i_stoptime_strategy = stoptime_[i_strategy]
            km[i_strategy] = DPM_reversiblemodel_analysis_KM(i_stoptime_strategy, duration_5year)
        par = {'duration': duration_5year, 'xtick step': 300, 'totalnum': len(stoptime_[strategy_name_[0]])}
        DPM_reversiblemodel_plot_km_multi(km, par, '')
        plt.savefig(filename_km_, format='pdf', bbox_inches='tight')
        plt.close()
        return

    def DPM_reversiblemodel_fun_analysis_4(stoptime_, strategy_name_, strategy_namepool_):
        stoptime_ = [dict(zip(stoptime_, x)) for x in zip(*stoptime_.values())]
        with tqdm(total=len(stoptime_), ncols=150, desc='running') as pbar:
            for i_stoptime in stoptime_:
                category['paramID'].extend([i_stoptime['paramID']])
                i_stoptime = {key: i_stoptime[key] for key in strategy_namepool_}
                i_stoptime_max = max(i_stoptime.values())
                i_strategy_namemax = [k for k, v in i_stoptime.items() if v == i_stoptime_max]
                if strategy_name_[0] in i_strategy_namemax and len(i_strategy_namemax) == 1:
                    category['category'].extend([0])
                elif strategy_name_[1] in i_strategy_namemax and len(i_strategy_namemax) == 1:
                    category['category'].extend([1])
                elif set(i_strategy_namemax) == set(strategy_name_):
                    category['category'].extend([2])
                elif strategy_name_[0] in i_strategy_namemax or strategy_name_[1] in i_strategy_namemax:
                    category['category'].extend([3])
                else:
                    category['category'].extend([4])
                pbar.update(1)
        return

    def DPM_reversiblemodel_fun_analysis_4_1(stoptime_, strategy_name_, strategy_namepool_):
        assert set(list(itertools.chain(*strategy_name_))).issubset(set(strategy_namepool_))
        stoptime_ = [dict(zip(stoptime_, x)) for x in zip(*stoptime_.values())]
        with tqdm(total=len(stoptime_), ncols=150, desc='running') as pbar:
            category['category'] = list()
            for ii, i_stoptime in enumerate(stoptime_):
                category['paramID'].extend([i_stoptime['paramID']])
                i_stoptime = {key: i_stoptime[key] for key in strategy_namepool_}
                i_stoptime_max = max(i_stoptime.values())
                i_strategy_namemax = [k for k, v in i_stoptime.items() if v == i_stoptime_max]
                i_category = 0
                for i_ in range(len(strategy_name_)):
                    if set(strategy_name_[i_]).intersection(set(i_strategy_namemax)):
                        i_category = i_category + 10 ** (len(strategy_name_)-i_-1)
                category['category'].extend([f'{i_category:0{len(strategy_name_)}}'])
                pbar.update(1)
            category_set = ["".join(seq) for seq in itertools.product("01", repeat=len(strategy_name_))]
            assert set(category['category']).issubset(set(category_set))
        return

    # --------begin------- #
    pathsave_figure = os.path.join(pathload, 'figure')
    if not os.path.exists(pathsave_figure):
        os.makedirs(pathsave_figure)

    pathsave_figure_improve_better = os.path.join(pathsave_figure, 'improve_better')
    if not os.path.exists(pathsave_figure_improve_better):
        os.makedirs(pathsave_figure_improve_better)

    pathsave_figure_improve_sigbetter = os.path.join(pathsave_figure, 'improve_sigbetter')
    if not os.path.exists(pathsave_figure_improve_sigbetter):
        os.makedirs(pathsave_figure_improve_sigbetter)

    pathsave_figure_paranalysis = os.path.join(pathsave_figure, 'paranalysis')
    if not os.path.exists(pathsave_figure_paranalysis):
        os.makedirs(pathsave_figure_paranalysis)

    pathsave_figure_stoptimeanalysis = os.path.join(pathsave_figure, 'stoptimeanalysis')
    if not os.path.exists(pathsave_figure_stoptimeanalysis):
        os.makedirs(pathsave_figure_stoptimeanalysis)

    pathsave_figure_doseanalysis = os.path.join(pathsave_figure, 'doseanalysis')
    if not os.path.exists(pathsave_figure_doseanalysis):
        os.makedirs(pathsave_figure_doseanalysis)

    pathsave_figure_popanalysis = os.path.join(pathsave_figure, 'popanalysis')
    if not os.path.exists(pathsave_figure_popanalysis):
        os.makedirs(pathsave_figure_popanalysis)

    pathsave_figure_Nflowanalysis = os.path.join(pathsave_figure, 'Nflowanalysis')
    if not os.path.exists(pathsave_figure_Nflowanalysis):
        os.makedirs(pathsave_figure_Nflowanalysis)

    pathsave_figure_sensitivity = os.path.join(pathsave_figure, 'sensitivity')
    if not os.path.exists(pathsave_figure_sensitivity):
        os.makedirs(pathsave_figure_sensitivity)

    pathsave_figure_sankeydiagram = os.path.join(pathsave_figure, 'sankeydiagram')
    if not os.path.exists(pathsave_figure_sankeydiagram):
        os.makedirs(pathsave_figure_sankeydiagram)

    filename_para = os.path.join(pathload, 'result_para.pckl')
    filename_stoptime = os.path.join(pathload, 'result_stoptime.pckl')
    filename_dosage = os.path.join(pathload, 'result_dosage.pckl')
    filename_Nflow = os.path.join(pathload, 'result_Nflow.pckl')
    filename_hz = os.path.join(pathload, 'hz.pckl')
    filename_category = os.path.join(pathload, 'category.pckl')

    filename_para_mutonly = os.path.join(pathload_mutonly, 'result_para.pckl')
    filename_category_mutonly = os.path.join(pathload_mutonly, 'category.pckl')
    filename_para_reveronly = os.path.join(pathload_reveronly, 'result_para.pckl')
    filename_category_reveronly = os.path.join(pathload_reveronly, 'category.pckl')

    filename_stoptimetable = os.path.join(pathsave_figure, 'stoptimetable.pdf')
    filename_hzratio = os.path.join(pathsave_figure, f'hzratio.pdf')
    filename_km = os.path.join(pathsave_figure, f'km.pdf')

    file_format = '.csv'
    file_list = os.listdir(pathload)
    partial_filename_para = ['_result_para', file_format]
    file_para = [filename for filename in file_list if all([x in filename for x in partial_filename_para])]
    partial_filename_stoptime = ['_result_stopt', file_format]
    file_stoptime = [filename for filename in file_list if all([x in filename for x in partial_filename_stoptime])]
    partial_filename_dosage = ['_result_dosage', file_format]
    file_dosage = [filename for filename in file_list if all([x in filename for x in partial_filename_dosage])]
    partial_filename_pop = ['_result_pop', file_format]
    file_pop = [filename for filename in file_list if all([x in filename for x in partial_filename_pop])]
    partial_filename_Nflow = ['_result_Nflow', '.pckl']
    file_Nflow = [filename for filename in file_list if all([x in filename for x in partial_filename_Nflow])]

    strategy_name_full = ['strategy1 cycle', 'strategy2 cycle']
    strategy_name_fullpool = ['mono drug1',
                              'mono drug2',
                              'mono half',
                              'best cycle',
                              'strategy0',
                              'strategy1',
                              'strategy1 cycle',
                              'strategy2',
                              'strategy2 cycle']
    strategy_name_mut = ['strategy1', 'strategy2']
    strategy_name_mutpool = ['mono drug1',
                             'mono drug2',
                             'mono half',
                             'best cycle',
                             'strategy0',
                             'strategy1',
                             'strategy1 cycle',
                             'strategy2',
                             'strategy2 cycle']
    strategy_name_rever = ['best cycle', 'strategy1 cycle']
    strategy_name_reverpool = ['mono drug1',
                               'mono drug2',
                               'mono half',
                               'best cycle',
                               'strategy0',
                               'strategy1',
                               'strategy1 cycle']

    # If there are no files found in the current directory, exit.
    assert len(file_stoptime) == len(file_para) == len(file_dosage)
    try:
        filename_i = file_stoptime[0].split('_')
        for i, i_val in enumerate(filename_i):
            try:
                int(i_val)
                break
            except ValueError:
                pass
        file_stoptime.sort(key=lambda x: int(x.split('_')[i]))
        file_dosage.sort(key=lambda x: int(x.split('_')[i]))
        file_pop.sort(key=lambda x: int(x.split('_')[i]))
        file_Nflow.sort(key=lambda x: int(x.split('_')[i]))
    except ValueError:
        pass

    # Read csv stopping time file
    stoptime = DPM_reversiblemodel_miscellaneous_dictlist(['paramID']+strategy_name)
    dosage = DPM_reversiblemodel_miscellaneous_dictlist(['paramID'] + strategy_name)
    bool_dosage = DPM_reversiblemodel_miscellaneous_dictlist(['paramID'] + ['bool_diff_dosage', 'bool_same_dosage'])

    para = DPM_reversiblemodel_miscellaneous_dictlist(['paramID'] + pd_colname + rate_key)
    category = DPM_reversiblemodel_miscellaneous_dictlist(['paramID'] + ['category'])
    Nflow = {'valcol': dict(zip(strategy_name, [np.zeros((len(cellflow), duration_5year))]*len(strategy_name))),
             'valcolextend': dict(zip(strategy_name, [np.zeros((len(cellflow), duration_5year))] * len(strategy_name))),
             'num': dict(zip(strategy_name, [np.zeros(duration_5year)]*len(strategy_name)))}

    if not os.path.isfile(filename_para) or \
            not os.path.isfile(filename_stoptime) or \
            not os.path.isfile(filename_dosage) or \
            not os.path.isfile(filename_category) or \
            ((not mutonly) and (not reveronly) and (not os.path.isfile(filename_Nflow))):
        for i, i_file_stoptime in enumerate(file_stoptime):
            DPM_reversiblemodel_fun_analysis_1(i_file_stoptime)
            print(i, len(file_stoptime), round(i/len(file_stoptime), 2))
        if mutonly:
            DPM_reversiblemodel_fun_analysis_4(stoptime, strategy_name_mut, strategy_name_mutpool)
            strategy_name_mut = [['best cycle'],
                                 ['strategy0'],
                                 ['strategy1', 'strategy1 cycle'],
                                 ['strategy2', 'strategy2 cycle']]
            DPM_reversiblemodel_fun_analysis_4_1(stoptime, strategy_name_mut, strategy_name_mutpool)
        elif reveronly:
            DPM_reversiblemodel_fun_analysis_4_1(stoptime, strategy_name_rever, strategy_name_reverpool)
            strategy_name_rever = [['best cycle'],
                                   ['strategy0'],
                                   ['strategy1', 'strategy1 cycle']]
            DPM_reversiblemodel_fun_analysis_4_1(stoptime, strategy_name_rever, strategy_name_reverpool)
        else:
            DPM_reversiblemodel_fun_analysis_4(stoptime, strategy_name_full, strategy_name_fullpool)
            # 16 by 16 table, list all possibilities
            strategy_name_full = [['best cycle'], ['strategy0'], ['strategy1', 'strategy1 cycle'], ['strategy2', 'strategy2 cycle']]
            DPM_reversiblemodel_fun_analysis_4_1(stoptime, strategy_name_full, strategy_name_fullpool)

        with bz2.BZ2File(filename_para, 'wb') as f:
            pickle.dump(para, f)
        with bz2.BZ2File(filename_stoptime, 'wb') as f:
            pickle.dump(stoptime, f)
        with bz2.BZ2File(filename_dosage, 'wb') as f:
            pickle.dump(dosage, f)
        with bz2.BZ2File(filename_category, 'wb') as f:
            pickle.dump(category, f)
        if (not mutonly) and (not reveronly):
            with bz2.BZ2File(filename_Nflow, 'wb') as f:
                pickle.dump(Nflow, f)
    else:
        with bz2.BZ2File(filename_para, 'rb') as f:
            para = pickle.load(f)
        with bz2.BZ2File(filename_stoptime, 'rb') as f:
            stoptime = pickle.load(f)
        with bz2.BZ2File(filename_dosage, 'rb') as f:
            dosage = pickle.load(f)
        with bz2.BZ2File(filename_category, 'rb') as f:
            category = pickle.load(f)
        if (not mutonly) and (not reveronly):
            with bz2.BZ2File(filename_Nflow, 'rb') as f:
                Nflow = pickle.load(f)

    # ------------------------ sankeydiagram --------------------------#
    if os.path.isfile(filename_para_mutonly) and os.path.isfile(filename_para_reveronly) \
            and os.path.isfile(filename_category_mutonly) and os.path.isfile(filename_category_reveronly) and \
            (not mutonly) and (not reveronly):
        with bz2.BZ2File(filename_para_mutonly, 'rb') as f:
            para_mutonly = pickle.load(f)
        with bz2.BZ2File(filename_para_reveronly, 'rb') as f:
            para_reveronly = pickle.load(f)
        with bz2.BZ2File(filename_category_mutonly, 'rb') as f:
            category_mutonly = pickle.load(f)
        with bz2.BZ2File(filename_category_reveronly, 'rb') as f:
            category_reveronly = pickle.load(f)

        # mut only
        para_name_source = ['g', 'R1ratio', 'R2ratio', 'T1', 'T2', 'S_01_1', 'S_01_2', 'S_11_1', 'S_11_2']
        name = 'mut to full'
        sankey_mut = DPM_reversiblemodel_analysis_sankeydiagram(para_mutonly,
                                                                para_name_source,
                                                                para,
                                                                category_mutonly,
                                                                category,
                                                                name,
                                                                pathsave_figure_sankeydiagram)
        # rever only
        para_name_source = ['g', 'alpha1', 'theta1', 'kappa1', 'mu1',
                            'alpha2', 'theta2', 'kappa2', 'mu2',
                            'S_01_1', 'S_01_2', 'S_00_1', 'S_00_2']
        name = 'rever to full'
        sankey_rever = DPM_reversiblemodel_analysis_sankeydiagram(para_reveronly,
                                                                  para_name_source,
                                                                  para,
                                                                  category_reveronly,
                                                                  category,
                                                                  name,
                                                                  pathsave_figure_sankeydiagram)

        name = 'mut to full'
        DPM_reversiblemodel_plot_sankeydiagram(sankey_mut, category_mutonly, name, pathsave_figure_sankeydiagram)

        name = 'rever to full'
        DPM_reversiblemodel_plot_sankeydiagram(sankey_rever, category_reveronly, name, pathsave_figure_sankeydiagram)

    if (not mutonly) and (not reveronly):
        DPM_reversiblemodel_analysis_sensitivity(stoptime, para, strategy_name, pathsave_figure_sensitivity)

        filename_sensitivity_para_df = os.path.join(pathsave_figure_sensitivity, 'sensitivity_para_df.pckl')
        filename_sensitivity_set1 = os.path.join(pathsave_figure_sensitivity, 'sensitivity_1.pckl')
        DPM_reversiblemodel_analysis_sensitivity_αθμ(stoptime, para, filename_sensitivity_set1, filename_sensitivity_para_df)

        filename_sensitivity_set2 = os.path.join(pathsave_figure_sensitivity, 'sensitivity_2.pckl')
        DPM_reversiblemodel_analysis_sensitivity_nonαθμ(stoptime, para, filename_sensitivity_set2, filename_sensitivity_para_df)

        if os.path.exists(filename_sensitivity_set1) and \
           os.path.exists(filename_sensitivity_set2) and \
           os.path.exists(filename_sensitivity_para_df):
            DPM_reversiblemodel_analysis_sensitivityplot(stoptime, strategy_name,
                                                         filename_sensitivity_para_df,
                                                         filename_sensitivity_set1,
                                                         filename_sensitivity_set2,
                                                         pathsave_figure_sensitivity)

    if (not mutonly) and (not reveronly):
        DPM_reversiblemodel_analysis_Nflow(Nflow, strategy_name, para, pathsave_figure_Nflowanalysis)

    info, para_ex, strategy_name_run = DPM_reversiblemodel_analysis_stoptime(stoptime,
                                                                             para,
                                                                             dosage,
                                                                             strategy_name,
                                                                             duration_5year,
                                                                             mutonly,
                                                                             reveronly)

    if not os.path.isfile(filename_hz):
        hz = DPM_reversiblemodel_analysis_hz(info['stoptime_df'].to_dict(orient='list'), strategy_name, duration_5year)
        with bz2.BZ2File(filename_hz, 'wb') as f:
            pickle.dump(hz, f)
    else:
        with bz2.BZ2File(filename_hz, 'rb') as f:
            hz = pickle.load(f)

    '''plot examples'''
    DPM_reversiblemodel_fun_sim1par(strategy_name_run, para_ex, pathsave_figure_popanalysis, mutonly, reveronly)

    '''plot and summary info'''
    DPM_reversiblemodel_fun_analysis_2(info,
                                       hz['hz'],
                                       strategy_name,
                                       filename_stoptimetable,
                                       filename_hzratio)

    '''plot km'''
    DPM_reversiblemodel_fun_analysis_3(info['stoptime_df'].to_dict(orient='list'), strategy_name, filename_km)

    '''analysis dose'''
    _ = DPM_reversiblemodel_analysis_dose(dosage, strategy_name, info, stoptime, pathsave_figure_doseanalysis)
    return
