from DPM_reversiblemodel_lib import np, math, csv, compress, groupby, os, re, pd, tqdm
from DPM_reversiblemodel_constant import Stepsize, \
    duration_5year,\
    strategy_default, \
    σ1_full, \
    σ2_full, \
    σ_half, \
    Num_drug, \
    teatmentname, \
    Numstep, \
    pd_colname


def DPM_reversiblemodel_miscellaneous_dictlist(keys):
    val = dict(zip(keys, [None]*len(keys)))
    for i_key in val.keys():
        val[i_key] = []
    return val


def DPM_reversiblemodel_stime(N):
    if all(N[:, -1] < 1):
        return duration_5year
    else:
        return min([N.shape[1] - 1, duration_5year])


def DPM_reversiblemodel_cure(N):
    return all(N[:, -1] < 1)


def DPM_reversiblemodel_miscellaneous_selectFalse(flaglist, *args):
    ind = list(filter(lambda i: not flaglist[i], range(len(flaglist))))
    out = []
    for i_val in args:
        out.append([i_val[i] for i in ind])
    out.append([flaglist[i] for i in ind])
    return out


def DPM_reversiblemodel_miscellaneous_selectTrue(flaglist, *args):
    ind = np.where(flaglist)[0]
    out = []
    for i_val in args:
        out.append([i_val[i] for i in ind])
    out.append([flaglist[i] for i in ind])
    return out


def DPM_reversiblemodel_miscellaneous_selectdose(survi_time, d):
    max_survi = np.amax(survi_time)
    index_max = np.argwhere(survi_time == max_survi).flatten()
    d_out = []
    for i, i_d in enumerate(d):
        if i in index_max:
            if d_out:
                if isinstance(i_d, np.ndarray):
                    flag = [not np.array_equal(i_d, i_din) if isinstance(i_din, np.ndarray) else True for i_din in d_out]
                elif isinstance(i_d, list):
                    flag = [not i_d == i_din if isinstance(i_d, list) else True for i_din in d_out]
                else:
                    raise TypeError('type of d is not correct, either a np array or a list.')
                all(flag) and isinstance(i_d, np.ndarray) and d_out.append(i_d)
                all(flag) and isinstance(i_d, list) and d_out.extend(i_d)
            else:
                isinstance(i_d, np.ndarray) and d_out.append(i_d)
                isinstance(i_d, list) and d_out.extend(i_d)
    return d_out, max_survi


def DPM_reversiblemodel_miscellaneous_findtreat(d, survi_time):
    d_out = []
    if isinstance(d, np.ndarray):
        while d.any():
            d_i, d = d[:, :Stepsize], d[:, Stepsize:]
            d_i = DPM_reversiblemodel_miscellaneous_matchtreat(d_i)
            d_out.append(d_i)
    if isinstance(d, str):
        num_d = math.ceil(survi_time/Stepsize)
        if d in ['12w σ1 1th', '12w σ2 1th']:
            if d == '12w σ1 1th':
                d_out = ['σ1 full', 'σ2 full'] * math.floor(num_d/2) + num_d % 2 * ['σ1 full']
            elif d == '12w σ2 1th':
                d_out = ['σ2 full', 'σ1 full'] * math.floor(num_d/2) + num_d % 2 * ['σ2 full']
            else:
                raise ValueError('type of d is not correct.')
        else:
            d_out = [d] * num_d
    return d_out


def DPM_reversiblemodel_miscellaneous_matchtreat(d):
    d = d[0, :]
    num_change = sum(np.diff(d) != 0)
    if d[0] == 1:
        if num_change == 0:
            return 'σ1 full'
        elif num_change == 1:
            return '3w σ1 1th'
        elif num_change in [2, 3]:
            raise ValueError('type of drug 1 dose is not correct.')
            # return '3w σ1 1th'
        else:
            raise ValueError('type of drug 1 dose is not correct.')
    elif d[0] == 0.5:
        if num_change == 0:
            return 'half'
        else:
            raise ValueError('type of half dose is not correct.')
    elif d[0] == 0:
        if num_change == 0:
            return 'σ2 full'
        elif num_change == 1:
            return '3w σ2 1th'
        elif num_change in [2, 3]:
            raise ValueError('type of drug 1 dose is not correct.')
            # return '3w σ2 1th'
        else:
            raise ValueError('type of drug 2dose is not correct.')
    else:
        raise ValueError('type of d is not correct.')


def DPM_reversiblemodel_miscellaneous_dfdose2dose(df_dose, survi):
    assert df_dose['Strategy name'] in strategy_default
    if df_dose['Strategy name'] == 'mono drug1':
        d = np.tile(np.array([σ1_full], dtype=float).T, (1, survi))
    elif df_dose['Strategy name'] == 'mono drug2':
        d = np.tile(np.array([σ2_full], dtype=float).T, (1, survi))
    elif df_dose['Strategy name'] == 'mono half':
        d = np.tile(np.array([σ_half], dtype=float).T, (1, survi))
    elif df_dose['Strategy name'] == 'best cycle':
        d = df_dose.values[2:].astype(str)
        dall = []
        for i_d in d:
            if i_d != 'nan':
                i_d.replace('sigma', 'σ')
                i_d = DPM_reversiblemodel_miscellaneous_cycletreat2dose(i_d, Numstep)
                i_d = i_d[:, :survi]
                dall.append(i_d)
        return dall, [survi] * len(dall)
    else:
        d = df_dose.values[2:]
        d = DPM_reversiblemodel_miscellaneous_treat2dose(d)
    d = d[:, :survi]
    return d, survi


def DPM_reversiblemodel_miscellaneous_treat2dose(treat):
    d = np.zeros([Num_drug, 0], dtype=float)
    for i_treat in treat:
        if i_treat == '-1':
            break
        i_treat = i_treat.replace('sigma', 'σ')
        assert i_treat in teatmentname
        if i_treat == 'σ1 full':
            i_d = np.tile(np.array([σ1_full], dtype=float).T, (1, Stepsize))
            d = np.append(d, i_d, 1)
        elif i_treat == 'σ2 full':
            i_d = np.tile(np.array([σ2_full], dtype=float).T, (1, Stepsize))
            d = np.append(d, i_d, 1)
        elif i_treat == 'half':
            i_d = np.tile(np.array([σ_half], dtype=float).T, (1, Stepsize))
            d = np.append(d, i_d, 1)
        elif i_treat in ['3w σ1 1th', '3w σ2 1th']:
            d = DPM_reversiblemodel_miscellaneous_cycletreat2dose(i_treat, 1)
            # i_duration = int(i_treat[0]) * 7
            # σ_cycle = [σ1_full, σ2_full] if 'σ1' in i_treat else [σ2_full, σ1_full]
            # cycle_duration = 2 * i_duration
            # σ_1th = np.array([σ_cycle[0]]).T
            # σ_2th = 1 - σ_1th
            # σ_1th = np.repeat(σ_1th, i_duration, axis=1)
            # σ_2th = np.repeat(σ_2th, i_duration, axis=1)
            # i_d = np.concatenate((σ_1th, σ_2th), axis=1)
            # # num_cycle = int(Stepsize / cycle_duration)
            # # d = np.tile(d, num_cycle)
            # assert Stepsize % cycle_duration == 0
            # assert i_d.shape[1] == Stepsize
            # d = np.append(d, i_d, 1)
    return d


def DPM_reversiblemodel_miscellaneous_cycletreat2dose(i_treat, num_step):
    try:
        i_duration = int(i_treat[:2]) * 7
    except ValueError:
        i_duration = int(i_treat[0]) * 7
    σ_cycle = [σ1_full, σ2_full] if 'σ1' in i_treat else [σ2_full, σ1_full]
    cycle_duration = 2 * i_duration
    σ_1th = np.array([σ_cycle[0]]).T
    σ_2th = 1 - σ_1th
    σ_1th = np.repeat(σ_1th, i_duration, axis=1)
    σ_2th = np.repeat(σ_2th, i_duration, axis=1)
    i_d = np.concatenate((σ_1th, σ_2th), axis=1)
    num_cycle = int(Stepsize / cycle_duration)
    i_d = np.tile(i_d, num_cycle) if num_cycle != 0 else i_d
    assert Stepsize % cycle_duration == 0 or Stepsize % cycle_duration == Stepsize
    assert i_d.shape[1] == Stepsize or i_d.shape[1] == 2*Stepsize
    i_d = np.tile(i_d, num_step)
    return i_d


def DPM_reversiblemodel_miscellaneous_treat_mono_or_half(stoptime, treat):
    treat = treat + ';'
    treat = treat * math.ceil(stoptime/Stepsize)
    if math.ceil(stoptime/Stepsize) < math.ceil(duration_5year/Stepsize):
        num = math.ceil(duration_5year/Stepsize) - math.ceil(stoptime/Stepsize)
        treat = treat + '-1;' * num
    return treat


def DPM_reversiblemodel_miscellaneous_writehead(filename, head):
    with open(filename, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(head)


def DPM_reversiblemodel_miscellaneous_indextrueinlist(list_):
    return list(compress(range(len(list_)), list_))


def DPM_reversiblemodel_miscellaneous_allequal(data):
    g = groupby(data)
    return next(g, True) and not next(g, False)


def DPM_reversiblemodel_miscellaneous_d_changetime(d):
    d_diff = np.diff(d)
    d_diff_boolen = d_diff.astype(bool)
    if d_diff_boolen.ndim != 1:
        pos_changedose = np.sum(d_diff_boolen, axis=0).astype(bool)
    else:
        pos_changedose = d_diff_boolen
    pos_changedose, = np.where(pos_changedose)
    # Position where drug dose changes.
    pos_changedose += 1
    # At the beginning, index=0, drug dose always changes. It can be seen as the start of applying treatment, previous drug doses are all 0.
    pos_changedose = np.insert(pos_changedose, 0, 0)
    return pos_changedose


def DPM_reversiblemodel_miscellaneous_renameparamID():
    pathload_par = './reversible_parcsv/'
    pathload_result = './reversible_csvresult/'
    file_format = '.csv'
    file_list = os.listdir(pathload_result)
    partial_filename_para = ['_result_para', file_format]
    file_para = [filename for filename in file_list if all([x in filename for x in partial_filename_para])]
    partial_filename_stoptime = ['_result_stopt', file_format]
    file_stoptime = [filename for filename in file_list if all([x in filename for x in partial_filename_stoptime])]
    partial_filename_dosage = ['_result_dosage', file_format]
    file_dosage = [filename for filename in file_list if all([x in filename for x in partial_filename_dosage])]
    file_stoptime.sort(key=lambda x: int(x.split('_')[3]))
    file_para.sort(key=lambda x: int(x.split('_')[3]))
    file_dosage.sort(key=lambda x: int(x.split('_')[3]))

    file_list_all = os.listdir(pathload_par)
    filename_pattern = 'par_reversiblemodel_para'
    file_list = [filename for filename in file_list_all if re.search(filename_pattern, filename) is not None]
    file_list.sort(key=lambda x: int(x.split('_')[3]))

    assert len(file_para) == len(file_stoptime) == len(file_dosage) == len(file_list)
    with tqdm(total=len(file_list), ncols=150) as pbar:
        for i in range(len(file_list)):
            i_file_para = os.path.join(pathload_result, file_para[i])
            i_file_stoptime = os.path.join(pathload_result, file_stoptime[i])
            i_file_dosage = os.path.join(pathload_result, file_dosage[i])
            i_file_csv = os.path.join(pathload_par, file_list[i])

            i_para = pd.read_csv(i_file_para)
            i_stoptime = pd.read_csv(i_file_stoptime)
            i_dosage = pd.read_csv(i_file_dosage, low_memory=False)
            i_csv = pd.read_csv(i_file_csv)

            for j in range(len(i_para)):
                j_i_para = i_para.iloc[j, :]
                j_i_csv = i_csv.iloc[j, :]

                for k_name in pd_colname:
                    assert j_i_para[k_name] == j_i_csv[k_name]

                i_para.loc[j, 'paramID'] = i_csv['paramID'][j]
                i_stoptime.loc[j, 'paramID'] = i_csv['paramID'][j]

            i_para.to_csv(i_file_para, index=False)
            i_stoptime.to_csv(i_file_stoptime, index=False)
            i_dosage.to_csv(i_file_dosage, index=False)
            pbar.update(1)
    return
