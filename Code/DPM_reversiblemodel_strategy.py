from DPM_reversiblemodel_constant import *
from DPM_reversiblemodel_model import \
    DPM_reversiblemodel_model_simexpm, \
    DPM_reversiblemodel_model_Sa
from DPM_reversiblemodel_miscellaneous import \
    DPM_reversiblemodel_miscellaneous_selectFalse, \
    DPM_reversiblemodel_miscellaneous_selectTrue, \
    DPM_reversiblemodel_stime, \
    DPM_reversiblemodel_miscellaneous_cycletreat2dose


def DPM_reversiblemodel_strategy_bestcyclesim(x0, par, duration, mutation):
    max_st, N_bestcycle, σ_bestcycle = 0, None, None
    for i_numweek in weekpercycle:
        i_N_σ1_1th, i_t_σ1_1th, i_N_σ2_1th, i_t_σ2_1th = \
            (DPM_reversiblemodel_strategy_cyclesim(x0, par, i_numweek, duration, mutation))
        i_bestcycle = [i_N_σ1_1th, i_N_σ2_1th]
        i_st = [DPM_reversiblemodel_stime(i_N_σ1_1th), DPM_reversiblemodel_stime(i_N_σ2_1th)]
        idx, i_max_st = np.argmax(np.array(i_st)), np.max(np.array(i_st))
        if i_max_st > max_st:
            N_bestcycle = i_bestcycle[idx]
            σ_bestcycle = str(i_numweek) + [' σ1 1th', ' σ2 1th'][idx]
            max_st = i_max_st
        elif i_max_st == max_st:
            if int(σ_bestcycle.split(' ')[0]) < i_numweek:
                N_bestcycle = i_bestcycle[idx]
                σ_bestcycle = str(i_numweek) + [' σ1 1th', ' σ2 1th'][idx]
    d_bestcycle = DPM_reversiblemodel_miscellaneous_cycletreat2dose(σ_bestcycle, Numstep)
    d_bestcycle = d_bestcycle[:, :N_bestcycle.shape[1] - 1]
    return N_bestcycle, d_bestcycle, σ_bestcycle


def DPM_reversiblemodel_strategy_cyclesim(x0, par, numweek, duration, mutation):
    σ_cycle = [σ1_full, σ2_full]
    tnow = 0
    N_σ1_first, flag_σ1_first, x0_σ1_first = np.array([x0]).T, True, x0
    N_σ2_first, flag_σ2_first, x0_σ2_first = np.array([x0]).T, True, x0
    while (tnow < duration) and any((flag_σ1_first, flag_σ2_first)):
        if flag_σ1_first:
            i_N_σ1_first, i_σ1first_flag_mortality, i_σ1first_flag_cure = \
                    DPM_reversiblemodel_strategy_simexmp(x0_σ1_first, σ_cycle[0], par, numweek*7, mutation)
            x0_σ1_first, N_σ1_first = i_N_σ1_first[:, -1], np.append(N_σ1_first, i_N_σ1_first[:, 1:], 1)
            flag_σ1_first = False if (i_σ1first_flag_mortality or i_σ1first_flag_cure) else True

        if flag_σ2_first:
            i_N_σ2_first, i_σ2first_flag_mortality, i_σ2first_flag_cure = \
                DPM_reversiblemodel_strategy_simexmp(x0_σ2_first, σ_cycle[1], par, numweek*7, mutation)
            x0_σ2_first, N_σ2_first = i_N_σ2_first[:, -1], np.append(N_σ2_first, i_N_σ2_first[:, 1:], 1)
            flag_σ2_first = False if (i_σ2first_flag_mortality or i_σ2first_flag_cure) else True

        σ_cycle.reverse()
        tnow = tnow + numweek * 7
    t_σ1_first, t_σ2_first = np.arange(N_σ1_first.shape[1]), np.arange(N_σ2_first.shape[1])
    return N_σ1_first, t_σ1_first, N_σ2_first, t_σ2_first


def DPM_reversiblemodel_strategy_simbydose(x0, d, par, mutation):
    tdiff = np.append(0, np.diff(np.arange(0, 1 + Simtimestep, Simtimestep)))
    N, flag_sim, x0, flag_mortality, flag_cure = np.array([x0]).T, True, x0, False, False
    for i_d in d.T:
        i_argsODE = (tuple(i_d), mutation, par)
        i_N, *_, flag_mortality, flag_cure = DPM_reversiblemodel_model_simexpm(i_argsODE, tdiff, x0)
        x0, N = i_N[:, -1], np.append(N, i_N[:, 1:], 1)
        if flag_mortality or flag_cure:
            break
    return N, N.sum(axis=0), flag_mortality, flag_cure


def DPM_reversiblemodel_strategy_simexmp(x0, σ, par, duration, mutation):
    t = np.arange(0, duration + Simtimestep, Simtimestep)
    tdiff = np.append(0, np.diff(t))
    argsODE = (σ, mutation, par)
    i_N, *_, flag_mortality, flag_cure = DPM_reversiblemodel_model_simexpm(argsODE, tdiff, x0)
    return i_N, flag_mortality, flag_cure


def DPM_reversiblemodel_strategy_bestcycle(x0, par, duration, mutation):
    N, cycle_name, N_σ1first, N_σ2first, survi_σ1first, survi_σ2first, σ_cycle = ([], [], [], [], [], [],
                                                                                  [σ1_full, σ2_full])
    for i_weekpercycle in weekpercycle:
        tnow = 0
        i_N_σ1_first, flag_σ1_first, x0_σ1_first = np.array([x0]).T, True, x0
        i_N_σ2_first, flag_σ2_first, x0_σ2_first = np.array([x0]).T, True, x0
        while (tnow < duration) and any((flag_σ1_first, flag_σ2_first)):
            if flag_σ1_first:
                ii_N_σ1_first, ii_σ1first_flag_mortality, ii_σ1first_flag_cure = \
                    DPM_reversiblemodel_strategy_simexmp(x0_σ1_first, σ_cycle[0], par, i_weekpercycle*7, mutation)
                x0_σ1_first, i_N_σ1_first = ii_N_σ1_first[:, -1], np.append(i_N_σ1_first, ii_N_σ1_first[:, 1:], 1)
                flag_σ1_first = False if (ii_σ1first_flag_mortality or ii_σ1first_flag_cure) else True

            if flag_σ2_first:
                ii_N_σ2_first, ii_σ2first_flag_mortality, ii_σ2first_flag_cure = \
                    DPM_reversiblemodel_strategy_simexmp(x0_σ2_first, σ_cycle[1], par, i_weekpercycle*7, mutation)
                x0_σ2_first, i_N_σ2_first = ii_N_σ2_first[:, -1], np.append(i_N_σ2_first, ii_N_σ2_first[:, 1:], 1)
                flag_σ2_first = False if (ii_σ2first_flag_mortality or ii_σ2first_flag_cure) else True

            σ_cycle.reverse()
            tnow = tnow + i_weekpercycle * 7
        N_σ1first.append(i_N_σ1_first)
        N_σ2first.append(i_N_σ2_first)
        survi_σ1first.append(DPM_reversiblemodel_stime(i_N_σ1_first))
        survi_σ2first.append(DPM_reversiblemodel_stime(i_N_σ2_first))

    max_survi = max(survi_σ1first+survi_σ2first)

    N_all = N_σ1first + N_σ2first
    for i in range(len(strategy_cycle_name)):
        if DPM_reversiblemodel_stime(N_all[i]) >= max_survi:
            N.append(N_all[i])
            cycle_name.append(strategy_cycle_name[i])

    return N[0], cycle_name


def DPM_reversiblemodel_strategy_0(x0, par, Simduration, mutation):
    # Strategy 0: Current personalized medicine:
    nadir = sum(x0)
    tnow = 0.0
    N_strategy0_i, d_strategy0_i = np.array([x0]).T, np.zeros([Num_drug, 0], dtype=float)
    drug_used = np.zeros([Num_drug], dtype=int)
    treat_i = np.zeros([Num_drug, 1], dtype=int)

    Sa = DPM_reversiblemodel_model_Sa(par)

    # Index of the majority cell type.
    index_majority_cell = np.argmax(x0)
    index_first_use_drug = np.argmax(Sa[index_majority_cell, :])

    treat_i[index_first_use_drug] = 1
    drug_used[index_first_use_drug] = 1

    t = np.arange(0, Stepsize + Simtimestep, Simtimestep)
    tdiff = np.append(0, np.diff(t))
    while tnow < Simduration:
        d_i = np.tile(treat_i, (1, t.shape[0] - 1))
        argsODE = (tuple(treat_i.flatten()), mutation, par)
        N_i, *_, flag_mortality_i, flag_cure_i = DPM_reversiblemodel_model_simexpm(argsODE, tdiff, x0)

        N_strategy0_i = np.append(N_strategy0_i, N_i[:, 1:], 1)
        d_strategy0_i = np.append(d_strategy0_i, d_i[:, :N_i.shape[1]-1], 1)

        # If total cell population is >= Limit_mortality or all cell type < 1, stop.
        if flag_mortality_i or flag_cure_i:
            break

        # If total cell population is bigger than 2 times of nadir, total cell population is bigger than Limit_radiologicdetection
        # (tumor can be detected), there are drugs have not been used (each drug is used only once) or
        # (2) If the total cell population reemerges from a level below the detection,
        # and there are drugs have not been used (each drug is used only once), switch to the other unused drugs.
        if ((N_i[:, -1].sum() >= 2 * nadir and N_i[:, -1].sum() >= Limit_radiologicdetection) or
                (N_i[:, 0].sum() < Limit_radiologicdetection <= N_i[:, -1].sum())) and any(drug_used == 0):
            treat_i = np.zeros([Num_drug, 1], dtype=int)
            index_drug_not_used = [i for i, x in enumerate(list(drug_used)) if x == 0]
            treat_i[index_drug_not_used] = 1
            drug_used[index_drug_not_used] = 1
            nadir = N_i[:, -1].sum()

        # If the total cell population is smaller than the current nadir and bigger than the Limit_radiologicdetection (tumor can be detected),
        # update the nadir.
        if Limit_radiologicdetection <= N_i[:, -1].sum() < nadir:
            nadir = N_i[:, -1].sum()
        # Update x0 and tnow.
        x0 = N_i[:, -1]
        tnow += Stepsize

    return N_strategy0_i, d_strategy0_i


def DPM_reversiblemodel_strategy_cycle(x0, σ_cycle, duration_σ, mutation, par, duration=Stepsize, LSsim=True):
    tnow, N, flag_mortality, flag_cure = 0, np.array([x0]).T, False, False

    cycle_duration = 2 * duration_σ
    σ_1th = np.array([σ_cycle[0]]).T
    σ_2th = 1 - σ_1th
    σ_1th = np.repeat(σ_1th, duration_σ,  axis=1)
    σ_2th = np.repeat(σ_2th, duration_σ, axis=1)
    d = np.concatenate((σ_1th, σ_2th), axis=1)
    num_cycle = int(Stepsize/cycle_duration)
    assert Stepsize % cycle_duration == 0
    d = np.tile(d, num_cycle)

    while tnow < duration and (not flag_mortality or not LSsim) and not flag_cure:
        i_N, flag_mortality, flag_cure = DPM_reversiblemodel_strategy_simexmp(x0, σ_cycle[0], par, duration_σ, mutation)
        x0, N = i_N[:, -1], np.append(N, i_N[:, 1:], 1)
        σ_cycle.reverse()
        tnow = tnow + duration_σ

    return N, d, flag_mortality, flag_cure


def DPM_reversiblemodel_strategy_1(x0, par, Simduration, mutation, use_cycle=False):
    # Strategy 1: Minimize the total cell population.
    # In each Stepsize, select the d_i that minimizes the total cell population.
    def DPM_reversiblemodel_strategy_1_1(σ_cycle):
        N_cycle, d_cycle, flag_mortality_cycle, flag_cure_cycle = \
            DPM_reversiblemodel_strategy_cycle(x0, σ_cycle, i_weekpercycle*7, mutation, par)
        DPM_reversiblemodel_strategy_1_2(N_cycle, d_cycle, flag_mortality_cycle, flag_cure_cycle)
        return

    def DPM_reversiblemodel_strategy_1_2(N_1, d_1, flag_mortality_1, flag_cure_cycle_1):
        N_i_total.append(N_1)
        d_i_total.append(d_1)
        N_i_end_total.append(N_1[:, -1].sum())
        flag_mortality.append(flag_mortality_1)
        flag_cure.append(flag_cure_cycle_1)
        t_sim.append(N_1.shape[1] - 1)
        return

    tnow = 0.0
    N_strategy1_i,  d_strategy1_i = np.array([x0]).T, np.zeros([Num_drug, 0], dtype=float)

    t_i = np.arange(0, Stepsize + Simtimestep, Simtimestep)
    tdiff = np.append(0, np.diff(t_i))
    while tnow < Simduration:
        N_i_total, d_i_total, N_i_end_total, t_sim, flag_mortality, flag_cure = [], [], [], [], [], []
        for i_dose_combination in dose_combination:
            treat_i = np.array([i_dose_combination], dtype=float).T
            d_i = np.tile(treat_i, (1, t_i.shape[0] - 1))

            argsODE = (i_dose_combination, mutation, par)
            N_i, *_, flag_mortality_i, flag_cure_i = DPM_reversiblemodel_model_simexpm(argsODE, tdiff, x0)

            DPM_reversiblemodel_strategy_1_2(N_i, d_i, flag_mortality_i, flag_cure_i)
        if use_cycle:
            for i_weekpercycle in weekpercycle_DPM:
                DPM_reversiblemodel_strategy_1_1([σ1_full, σ2_full])
                DPM_reversiblemodel_strategy_1_1([σ2_full, σ1_full])

        if any(flag_cure):
            N_i_total, d_i_total, flag_cure = DPM_reversiblemodel_miscellaneous_selectTrue(flag_cure, N_i_total, d_i_total)
            index_select = 0
        else:
            if not all(flag_mortality):
                # apply the strategy on the dose without causing mortality
                N_i_total, d_i_total, N_i_end_total, flag_mortality = \
                    DPM_reversiblemodel_miscellaneous_selectFalse(flag_mortality, N_i_total, d_i_total, N_i_end_total)

                # find the d_i that minimizes total cell population.
                index_select = np.argmin(N_i_end_total)
            else:
                # This means simulation of all the dose combination reach mortality. Then select the longest survival time.
                index_select = np.argmax(t_sim)

        N_i_minimum, d_i_minimum = N_i_total[index_select], d_i_total[index_select]

        N_strategy1_i = np.append(N_strategy1_i, N_i_minimum[:, 1:], 1)
        d_strategy1_i = np.append(d_strategy1_i, d_i_minimum[:, :N_i_minimum.shape[1]-1], 1)

        # If total cell population is >= Limit_mortality or all cell type < 1, stop.
        if flag_mortality[index_select] or flag_cure[index_select]:
            break

        # Update x0 and tnow.
        x0 = N_i_minimum[:, -1]
        tnow += Stepsize
    return N_strategy1_i, d_strategy1_i


def DPM_reversiblemodel_strategy_2_strategy_3(x0, mutation, par, use_cycle):
    def DPM_reversiblemodel_strategy_2_strategy_3_1(σ_cycle):
        N_cycle, d_cycle, flag_mortality_cycle, flag_cure_cycle = \
            DPM_reversiblemodel_strategy_cycle(x0, σ_cycle, i_weekpercycle * 7, mutation, par)
        DPM_reversiblemodel_strategy_2_strategy_3_2(N_cycle, d_cycle, flag_mortality_cycle, flag_cure_cycle)
        return

    def DPM_reversiblemodel_strategy_2_strategy_3_2(N_1, d_1, flag_mortality_1, flag_cure_1):
        N_i_total.append(N_1)
        d_i_total.append(d_1)
        N_i_end_total.append(N_1[:, -1].sum())
        N_i_end_multi_mut.append(N_1[ind_multi_mut, -1].sum())
        flag_mortality.append(flag_mortality_1)
        flag_cure.append(flag_cure_1)
        t_sim.append(N_1.shape[1] - 1)
        return

    t_i = np.arange(0, Stepsize + Simtimestep, Simtimestep)
    tdiff = np.append(0, np.diff(t_i))
    N_i_total, d_i_total, N_i_end_total, N_i_end_multi_mut, t_sim, flag_mortality, flag_cure = [], [], [], [], [], [], []
    for i_dose_combination in dose_combination:
        treat_i = np.array([i_dose_combination], dtype=float).T
        d_i = np.tile(treat_i, (1, t_i.shape[0] - 1))

        argsODE = (i_dose_combination, mutation, par)
        N_i, *_, flag_mortality_i, flag_cure_i = DPM_reversiblemodel_model_simexpm(argsODE, tdiff, x0)

        N_i_total.append(N_i)
        d_i_total.append(d_i)
        N_i_end_total.append(N_i[:, -1].sum())
        N_i_end_multi_mut.append(N_i[ind_multi_mut, -1].sum())
        flag_mortality.append(flag_mortality_i)
        flag_cure.append(flag_cure_i)
        t_sim.append(N_i.shape[1] - 1)
    if use_cycle:
        for i_weekpercycle in weekpercycle_DPM:
            DPM_reversiblemodel_strategy_2_strategy_3_1([σ1_full, σ2_full])
            DPM_reversiblemodel_strategy_2_strategy_3_1([σ2_full, σ1_full])

    return N_i_total, d_i_total, N_i_end_total, N_i_end_multi_mut, t_sim, flag_mortality, flag_cure


def DPM_reversiblemodel_strategy_2(x0, par, Simduration, mutation, use_cycle=False):
    # Strategy 2: Minimize the risk of incurable cells developing unless there is an immediate threat of mortality.
    # Strategy 2.2: threshold is 1e11
    tnow = 0.0
    N_strategy2_i, d_strategy2_i = np.array([x0]).T, np.zeros([Num_drug, 0], dtype=float)

    while tnow < Simduration:
        N_i_total, d_i_total, N_i_end_total, N_i_end_multi_mut, t_sim, flag_mortality, flag_cure = \
            DPM_reversiblemodel_strategy_2_strategy_3(x0, mutation, par, use_cycle)

        if any(flag_cure):
            N_i_total, d_i_total, flag_cure = DPM_reversiblemodel_miscellaneous_selectTrue(flag_cure, N_i_total,
                                                                                           d_i_total)
            index_select = 0
        else:
            if not all(flag_mortality):
                # apply the strategy on the dose without causing mortality
                N_i_total, d_i_total, N_i_end_total, N_i_end_multi_mut, flag_mortality = \
                    DPM_reversiblemodel_miscellaneous_selectFalse(flag_mortality, N_i_total, d_i_total, N_i_end_total,
                                                                  N_i_end_multi_mut)

                # If the current total cell population does not exceed the threshold, minimize the
                # multiply-resistant population.
                if sum(x0) <= strategy2threshold:
                    index_select = np.argmin(N_i_end_multi_mut)
                # If the current total cell population exceeds the threshold, minimize the total population.
                else:
                    index_select = np.argmin(N_i_end_total)
            else:
                index_select = np.argmax(t_sim)

        N_i_minimum, d_i_minimum = N_i_total[index_select], d_i_total[index_select]

        N_strategy2_i = np.append(N_strategy2_i, N_i_minimum[:, 1:], 1)
        d_strategy2_i = np.append(d_strategy2_i, d_i_minimum[:, :N_i_minimum.shape[1]-1], 1)

        # If total cell population is >= Limit_mortality or all cell type < 1, stop.
        if flag_mortality[index_select] or flag_cure[index_select]:
            break

        # Update x0 and tnow.
        x0 = N_i_minimum[:, -1]
        tnow += Stepsize
    return N_strategy2_i, d_strategy2_i


def DPM_reversiblemodel_strategy_stopX(N, d, maxthreshold):
    flag_mortality = False
    if N[:, -1].sum() >= maxthreshold:
        pos = np.argmax(np.sum(N, axis=0) >= maxthreshold)
        N, d = N[:, :pos], d[:, :pos]
        flag_mortality = True
    return N, d, flag_mortality


def DPM_reversiblemodel_strategy_simstep(x0, parode, d, mutation):
    N_op, d_op = np.array([x0]).T, np.zeros([Num_drug, 0], dtype=float)
    t_i = np.arange(0, Stepsize + Simtimestep, Simtimestep)
    tdiff = np.append(0, np.diff(t_i))
    treat = np.array([d/2, 1-d/2])
    for i, i_treat in enumerate(treat.T):
        i_treat = np.array([i_treat], dtype=float).T
        d_i = np.tile(i_treat, (1, t_i.shape[0] - 1))
        argsODE = (tuple(i_treat.flatten()), mutation, parode)
        N_i, *_, flag_mortality_i, flag_cure_i = DPM_reversiblemodel_model_simexpm(argsODE, tdiff, x0)

        N_op = np.append(N_op, N_i[:, 1:], 1)
        d_op = np.append(d_op, d_i[:, :N_i.shape[1]-1], 1)

        # If total cell population is >= Limit_mortality or all cell type < 1, stop.
        if flag_mortality_i or flag_cure_i:
            break
        x0 = N_op[:, -1]
    return N_op, d_op


def DPM_reversiblemodel_strategy_simstep_byname(x0, parode, d, mutation):
    N_op, d_op = np.array([x0]).T, np.zeros([Num_drug, 0], dtype=float)
    t_i = np.arange(0, Stepsize + Simtimestep, Simtimestep)
    tdiff = np.append(0, np.diff(t_i))
    for i, (i_key, i_val) in enumerate(d.items()):
        if i_val in ['6w σ1 1th', '6w σ2 1th', '3w σ1 1th', '3w σ2 1th']:
            if i_val[0] == '6':
                i_weekpercycle = 6
            elif i_val[0] == '3':
                i_weekpercycle = 3
            else:
                raise ValueError('type of d is not correct.')
            if 'σ1' in i_val:
                σ_cycle = [σ1_full, σ2_full]
            elif 'σ2' in i_val:
                σ_cycle = [σ2_full, σ1_full]
            else:
                raise ValueError('type of d is not correct.')
            N_i, d_i, flag_mortality_i, flag_cure_i = \
                DPM_reversiblemodel_strategy_cycle(x0, σ_cycle, i_weekpercycle * 7, mutation, parode)
        elif i_val in ['σ1 full', 'σ2 full', 'half']:
            if i_val == 'σ1 full':
                i_treat = σ1_full
            elif i_val == 'σ2 full':
                i_treat = σ2_full
            elif i_val == 'half':
                i_treat = σ_half
            else:
                raise ValueError('type of d is not correct.')
            i_treat = np.array([i_treat], dtype=float).T
            d_i = np.tile(i_treat, (1, t_i.shape[0] - 1))
            argsODE = (tuple(i_treat.flatten()), mutation, parode)
            N_i, *_, flag_mortality_i, flag_cure_i = DPM_reversiblemodel_model_simexpm(argsODE, tdiff, x0)
        else:
            raise ValueError('type of d is not correct.')
        N_op = np.append(N_op, N_i[:, 1:], 1)
        d_op = np.append(d_op, d_i[:, :N_i.shape[1] - 1], 1)

        # If total cell population is >= Limit_mortality or all cell type < 1, stop.
        if flag_mortality_i or flag_cure_i:
            break
        x0 = N_op[:, -1]
    return N_op, d_op
