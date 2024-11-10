from DPM_reversiblemodel_model import *
from DPM_reversiblemodel_miscellaneous import DPM_reversiblemodel_miscellaneous_d_changetime

color = ['black', 'gray', 'darkorange', 'red', 'darkgreen', 'deepskyblue', 'blueviolet', 'gold', 'c', 'm', 'y']


# def DPM_reversiblemodel_plot_rate_Svsrr(θ_step, α_step, κ):
#     def DPM_reversiblemodel_plot_rate_Svsrr_1(ax_1, E_1, label_1, σ1):
#         for i, i_α1 in enumerate(α):
#             E_1_α_i = E_1[str(i_α1)]
#             ax_1.plot(E_1_α_i['θ'],
#                       np.exp(-κ * (E_1_α_i['E+'] - E_1_α_i['E-'])),
#                       color=color[i],
#                       label=r'$e^{-κ\times(E^{+}-E^{-})}_{S->rr}, α=%s$' % (str(round(i_α1, 1))))
#             ax_1.plot(E_1_α_i['θ'],
#                       np.exp(-κ * (E_1_α_i['E-'] - E_1_α_i['E+'])),
#                       ls='--',
#                       color=color[i],
#                       label=r'$e^{-κ\times(E^{-}-E^{+})}_{rr->S}, α=%s$' % (str(round(i_α1, 1))))
#
#         ax_1.set_title(f'σ={σ1}, κ={κ}')
#         ax_1.set_xlabel('θ')
#         ax_1.set_xticks(np.arange(0, 1.1, .1))
#         ax_1.set_yscale("log", base=10)
#         ax_1.grid()
#         if label_1:
#             ax_σ_0.legend(frameon=False, bbox_to_anchor=(-0.1, 1), ncol=2)
#
#     θ = np.arange(0., 1. + θ_step, θ_step)
#     θ = θ[θ <= θ_max]
#     α = np.array([0.2, 0.5, 0.8])  # np.arange(0., 1., α_step)
#     α = α[α <= α_max]
#
#     # E+ under σ=0, 0.5, 1
#     Eplus_σ_0 = DPM_reversiblemodel_model_Eplus(0, θ)
#     Eplus_σ_half = DPM_reversiblemodel_model_Eplus(.5, θ)
#     Eplus_σ_full = DPM_reversiblemodel_model_Eplus(1., θ)
#
#     E_σ_0_α, E_σ_half_α, E_σ_full_α = dict(dict()), dict(dict()), dict(dict())
#     for i_α in α:
#         ind = np.flatnonzero(θ >= i_α)[0]
#         i_θ = θ[ind:]
#         # σ = 0
#         E_σ_0_α[str(i_α)] = {'θ': i_θ, 'E+': Eplus_σ_0[ind:], 'E-': DPM_reversiblemodel_model_Eminus(0., i_θ, i_α)}
#         # σ = 0.5
#         E_σ_half_α[str(i_α)] = {'θ': i_θ, 'E+': Eplus_σ_half[ind:], 'E-': DPM_reversiblemodel_model_Eminus(0.5, i_θ, i_α)}
#         # σ = 1
#         E_σ_full_α[str(i_α)] = {'θ': i_θ, 'E+': Eplus_σ_full[ind:], 'E-': DPM_reversiblemodel_model_Eminus(1., i_θ, i_α)}
#
#     plt.rcParams['font.size'] = 12
#     fig, (ax_σ_0, ax_σ_half, ax_σ_full) = plt.subplots(1, 3)
#
#     '''σ = 0'''
#     label = True
#     σ = 0
#     DPM_reversiblemodel_plot_rate_Svsrr_1(ax_σ_0, E_σ_0_α, label, σ)
#     '''σ = 0.5'''
#     label = False
#     σ = 0.5
#     DPM_reversiblemodel_plot_rate_Svsrr_1(ax_σ_half, E_σ_half_α, label, σ)
#     '''σ = 1.0'''
#     σ = 1.0
#     DPM_reversiblemodel_plot_rate_Svsrr_1(ax_σ_full, E_σ_full_α, label, σ)
#
#     fig.set_size_inches(18, 6)
#     plt.tight_layout(pad=1, w_pad=0.3, h_pad=1.1)
#     return


def DPM_reversiblemodel_plot_nonmut(t, n, title, leg_plot):
    n_ratio = np.zeros(n.shape)
    n_ratio[0, :] = n[0, :]/np.sum(n, 0)
    n_ratio[1, :] = n[1, :]/np.sum(n, 0)
    n_ratio[2, :] = n[2, :]/np.sum(n, 0)
    n_ratio[3, :] = n[3, :]/np.sum(n, 0)

    plt.rcParams.update({'font.size': 20})
    fig, (ax_n, ax_ratio) = plt.subplots(1, 2)

    ax_n.plot(t, n[0, :], color='g', lw=2, label='0101')
    ax_n.plot(t, n[1, :], color='k', lw=2, label='0001')
    ax_n.plot(t, n[2, :], color='b', lw=2, label='0100')
    ax_n.plot(t, n[3, :], color='r', lw=2, label='0000')
    ax_n.plot(t, np.sum(n, 0) * 1.5, color='m', lw=2, label='total')
    if leg_plot:
        leg = ax_n.legend(loc='best', ncol=1, frameon=1, fontsize=19)
        leg.draw_frame(False)
    ax_n.set_xlabel('Time[days]')
    ax_n.set_ylabel('Cell Number')
    ax_n.set_yscale('log', base=10)
    # ax_n.set_ylim([1, 1e13])
    ax_n.set_title(title)

    ax_ratio.plot(t, n_ratio[0, :], color='g', lw=2, label='$0101$')
    ax_ratio.plot(t, n_ratio[1, :], color='k', lw=2, label='$0001$')
    ax_ratio.plot(t, n_ratio[2, :], color='b', lw=2, label='$0100$')
    ax_ratio.plot(t, n_ratio[3, :], color='r', lw=2, label='$0000$')
    ax_ratio.set_xlabel('Time[days]')
    ax_ratio.set_ylabel('Cell type percentage')
    ax_ratio.set_ylim([0, 1])
    ax_ratio.set_title(title)

    fig.set_size_inches(17,  8)
    plt.tight_layout(pad=0.4, w_pad=0.4, h_pad=1.0)
    return


def DPM_reversiblemodel_plot_cellnumwithd(t, N, d, tmax, title, Ntype=None):
    assert Ntype in ('nonmut', 'mut')
    logbase = 10
    Legend_boxanchor = (0.5, 1.22)
    diffpts = DPM_reversiblemodel_miscellaneous_d_changetime(d)
    plt.rcParams.update({'font.size': 20})
    fig, ax = plt.subplots()
    if Ntype == 'nonmut':
        color_nonmut = ['k', 'g', 'b', 'r']
        Legend_colnum = 3
        Legend_order = np.array([0, 3, 1, 4, 2, 5])
        for i, i_idx in enumerate(ind_nonmut):
            ax.plot(t, N[i_idx, :], color=color_nonmut[i], label=celltype[i_idx])
    else:
        Legend_colnum = 4
        color_mut = ['k', 'g', 'b', 'r', 'm']
        Legend_order = np.array([0, 4, 1, 5, 2, 6, 3])
        for i, i_idx in enumerate(ind_mut):
            ax.plot(t, N[i_idx, :], color=color_mut[i], label=celltype[i_idx])

    ymaxval = max([math.ceil(math.log(np.max(np.sum(N, 0)), logbase)) + Ylim_add_plot,
                   math.ceil(math.log(Limit_mortality, logbase))+1 + Ylim_add_plot])
    ylimmin = Ylimmin_plot
    if ymaxval > math.log(Limit_mortality, logbase):
        plt.axhline(y=Limit_mortality, color='k', linestyle='dashed')
        plt.text(tmax*3/4, 2*Limit_mortality, 'Limit of mortality', fontsize=13)
    if ylimmin < 1:
        plt.axhline(y=1, color='k', linestyle='dashed')
        plt.text(tmax*3/4, 2*1, 'Limit of proliferation', fontsize=13)
    # plot dose
    diffpts = np.append(diffpts, t[-1])
    diffpts = np.append(diffpts, np.arange(0, math.ceil(t.max()/Stepsize))*Stepsize)
    diffpts = np.unique(diffpts)
    diffpts = diffpts.astype(int)
    diffpts.sort()
    for i_treat in range(len(diffpts)-1):
        d_i = d[:, diffpts[i_treat]:diffpts[i_treat+1]]
        d_i = np.unique(d_i, axis=1)
        d_i_sum = d_i.sum()
        y_basal = ymaxval - Ylim_add_plot
        for i_drug in range(Num_drug):
            if d_i[i_drug] != 0:
                plt.fill_between(t[diffpts[i_treat]:diffpts[i_treat+1]+1], logbase**y_basal,
                                 logbase**(y_basal+Ylim_add_plot*d_i[i_drug]/d_i_sum),
                                 color=color_drug_plot[i_drug])
                plt.hlines(y=logbase**(y_basal+Ylim_add_plot*d_i[i_drug]/d_i_sum), xmin=t[diffpts[i_treat]],
                           xmax=t[diffpts[i_treat+1]], color='black')
                # if i_treat != 0:
                plt.vlines(x=t[diffpts[i_treat]], ymin=logbase**(ymaxval-Ylim_add_plot), ymax=logbase**ymaxval, color='black')
                y_basal = y_basal + Ylim_add_plot*d_i[i_drug]/d_i_sum

    # Plot artificial label.
    y_basal = ymaxval-Ylim_add_plot
    d_lab = np.full(Num_drug, 1/Num_drug, dtype=float)
    t_lab = np.array(range(-int(t[-1]), -int(t[-1] / 2), 1))
    for i in range(Num_drug):
        plt.fill_between(t_lab, logbase**y_basal, logbase**(y_basal+Ylim_add_plot*d_lab[i]),
                         color=color_drug_plot[i], label=Lable_drug[i])
        y_basal = y_basal+logbase*d_lab[i]
    handles, labels = plt.gca().get_legend_handles_labels()
    leg = plt.legend([handles[idx] for idx in Legend_order], [labels[idx] for idx in Legend_order],
                     bbox_to_anchor=Legend_boxanchor, loc='upper center', ncol=Legend_colnum)
    leg.get_frame().set_linewidth(0)
    ax.set_xlabel('Time[days]')
    ax.set_ylabel('Cell Number')
    ax.set_yscale('log', base=logbase)
    ax.set_title(title)
    ax.set_xlim([Xlimmin_plot, tmax+Stepsize/4])
    # ax.set_xticks(diffpts[:-1])
    # ax.set_xticks(np.arange(0, math.ceil(t.max()/Stepsize))*Stepsize)
    xtick = np.arange(0, math.ceil(tmax/Stepsize))*Stepsize
    ax.set_xticks(xtick[::1])
    ax.grid()
    ax.set_ylim([ylimmin, logbase ** ymaxval])
    ylocmaj = mpl.ticker.LogLocator(base=logbase, numticks=len(range(-2, int(ymaxval), 2)))
    ax.yaxis.set_major_locator(ylocmaj)
    ylocmin = mpl.ticker.LogLocator(base=logbase, subs=np.arange(0, logbase)*1/logbase, numticks=100)
    ax.yaxis.set_minor_locator(ylocmin)
    ax.yaxis.set_minor_formatter(mpl.ticker.NullFormatter())
    fig.set_size_inches([14.5,  8.5])
    plt.tight_layout(pad=0.25, w_pad=0.25, h_pad=1.0)
    return


def DPM_reversiblemodel_plot_cellnumwithd_mutonly_reveronly(t, N, d, tmax, title, mutonly, reveronly):
    assert mutonly or reveronly is True
    logbase = 10
    Legend_boxanchor = (0.5, 1.22)
    diffpts = DPM_reversiblemodel_miscellaneous_d_changetime(d)
    plt.rcParams.update({'font.size': 20})
    fig, ax = plt.subplots()
    color_nonmut = ['k', 'g', 'b', 'r']
    Legend_colnum = 6
    Legend_order = np.array([0, 1, 2, 3, 4, 5])
    if mutonly:
        ind = ind_mutonly
    else:
        ind = ind_reveronly
    for i, i_idx in enumerate(ind):
        ax.plot(t, N[i_idx, :], color=color_nonmut[i], label=celltype[i_idx])

    ymaxval = max([math.ceil(math.log(np.max(np.sum(N, 0)), logbase)) + Ylim_add_plot,
                   math.ceil(math.log(Limit_mortality, logbase))+1 + Ylim_add_plot])
    ylimmin = Ylimmin_plot
    if ymaxval > math.log(Limit_mortality, logbase):
        plt.axhline(y=Limit_mortality, color='k', linestyle='dashed')
        plt.text(tmax*3/4, 2*Limit_mortality, 'Limit of mortality', fontsize=13)
    if ylimmin < 1:
        plt.axhline(y=1, color='k', linestyle='dashed')
        plt.text(tmax*3/4, 2*1, 'Limit of proliferation', fontsize=13)
    # plot dose
    diffpts = np.append(diffpts, t[-1])
    diffpts = np.append(diffpts, np.arange(0, math.ceil(t.max()/Stepsize))*Stepsize)
    diffpts = np.unique(diffpts)
    diffpts = diffpts.astype(int)
    diffpts.sort()
    for i_treat in range(len(diffpts)-1):
        d_i = d[:, diffpts[i_treat]:diffpts[i_treat+1]]
        d_i = np.unique(d_i, axis=1)
        d_i_sum = d_i.sum()
        y_basal = ymaxval - Ylim_add_plot
        for i_drug in range(Num_drug):
            if d_i[i_drug] != 0:
                plt.fill_between(t[diffpts[i_treat]:diffpts[i_treat+1]+1], logbase**y_basal,
                                 logbase**(y_basal+Ylim_add_plot*d_i[i_drug]/d_i_sum),
                                 color=color_drug_plot[i_drug])
                plt.hlines(y=logbase**(y_basal+Ylim_add_plot*d_i[i_drug]/d_i_sum), xmin=t[diffpts[i_treat]],
                           xmax=t[diffpts[i_treat+1]], color='black')
                # if i_treat != 0:
                plt.vlines(x=t[diffpts[i_treat]], ymin=logbase**(ymaxval-Ylim_add_plot), ymax=logbase**ymaxval, color='black')
                y_basal = y_basal + Ylim_add_plot*d_i[i_drug]/d_i_sum

    # Plot artificial label.
    y_basal = ymaxval-Ylim_add_plot
    d_lab = np.full(Num_drug, 1/Num_drug, dtype=float)
    t_lab = np.array(range(-int(t[-1]), -int(t[-1] / 2), 1))
    for i in range(Num_drug):
        plt.fill_between(t_lab, logbase**y_basal, logbase**(y_basal+Ylim_add_plot*d_lab[i]),
                         color=color_drug_plot[i], label=Lable_drug[i])
        y_basal = y_basal+logbase*d_lab[i]
    handles, labels = plt.gca().get_legend_handles_labels()
    leg = plt.legend([handles[idx] for idx in Legend_order], [labels[idx] for idx in Legend_order],
                     bbox_to_anchor=Legend_boxanchor, loc='upper center', ncol=Legend_colnum)
    leg.get_frame().set_linewidth(0)
    ax.set_xlabel('Time[days]')
    ax.set_ylabel('Cell Number')
    ax.set_yscale('log', base=logbase)
    ax.set_title(title)
    ax.set_xlim([Xlimmin_plot, tmax+Stepsize/4])
    # ax.set_xticks(diffpts[:-1])
    # ax.set_xticks(np.arange(0, math.ceil(t.max()/Stepsize))*Stepsize)
    xtick = np.arange(0, math.ceil(tmax/Stepsize))*Stepsize
    ax.set_xticks(xtick[::1])
    ax.grid()
    ax.set_ylim([ylimmin, logbase ** ymaxval])
    ylocmaj = mpl.ticker.LogLocator(base=logbase, numticks=len(range(-2, int(ymaxval)+2, 2)))
    ax.yaxis.set_major_locator(ylocmaj)
    ylocmin = mpl.ticker.LogLocator(base=logbase, subs=np.arange(0, logbase)*1/logbase, numticks=100)
    ax.yaxis.set_minor_locator(ylocmin)
    ax.set_yticks([10**-2, 10**1, 10**4, 10**7, 10**10, 10**13, 10**16])
    ax.set_ylim([10**-2, 10**16])
    ax.yaxis.set_minor_formatter(mpl.ticker.NullFormatter())
    plt.tight_layout(pad=0.25, w_pad=0.25, h_pad=1.0)
    fig.set_size_inches([14.5,  8.5])
    return


def DPM_reversiblemodel_plot_cellnumpercentagewithd(t, N, d, tmax, title, Ntype=None):
    assert Ntype in ('nonmut', 'mut')
    Legend_boxanchor = (0.5, 1.22)
    diffpts = DPM_reversiblemodel_miscellaneous_d_changetime(d)
    N_ratio = np.zeros(N.shape)
    for i in range(N.shape[0]):
        N_ratio[i, :] = N[i, :] / np.sum(N, 0) * 100
    plt.rcParams.update({'font.size': 20})
    fig, ax = plt.subplots()
    if Ntype == 'nonmut':
        Legend_colnum = 3
        Legend_order = np.array([0, 3, 1, 4, 2, 5])
        color_nonmut = ['k', 'g', 'b', 'r']
        for i, i_idx in enumerate(ind_nonmut):
            ax.plot(t, N_ratio[i_idx, :], color=color_nonmut[i], label=celltype[i_idx])
    else:
        Legend_colnum = 4
        color_mut = ['k', 'g', 'b', 'r', 'm']
        Legend_order = np.array([0, 4, 1, 5, 2, 6, 3])
        for i, i_idx in enumerate(ind_mut):
            ax.plot(t, N_ratio[i_idx, :], color=color_mut[i], label=celltype[i_idx])

    ylim_add_plot = 10
    ylimmin, ylimmax = 0, 100+ylim_add_plot
    # plot dose
    diffpts = np.append(diffpts, t[-1])
    diffpts = diffpts.astype(int)
    for i_treat in range(len(diffpts) - 1):
        d_i = d[:, diffpts[i_treat]:diffpts[i_treat+1]]
        d_i = np.unique(d_i, axis=1)
        d_i_sum = d_i.sum()
        y_basal = 101
        for i_drug in range(Num_drug):
            if d_i[i_drug] != 0:
                plt.fill_between(t[diffpts[i_treat]:diffpts[i_treat+1]+1], y_basal, y_basal+ylim_add_plot*d_i[i_drug]/d_i_sum,
                                 color=color_drug_plot[i_drug])
                plt.hlines(y=y_basal+ylim_add_plot*d_i[i_drug]/d_i_sum, xmin=t[diffpts[i_treat]],
                           xmax=t[diffpts[i_treat+1]], color='black')
                # if i_treat != 0:
                plt.vlines(x=t[diffpts[i_treat]], ymin=y_basal, ymax=ylimmax, color='black')
                y_basal = y_basal + ylim_add_plot*d_i[i_drug]/d_i_sum

    # Plot artificial label.
    y_basal = ylimmax-Ylim_add_plot
    d_lab = np.full(Num_drug, 1/Num_drug, dtype=float)
    t_lab = np.array(range(-int(t[-1]), -int(t[-1] / 2), 1))
    for i in range(Num_drug):
        plt.fill_between(t_lab, y_basal, y_basal+Ylim_add_plot*d_lab[i], color=color_drug_plot[i], label=Lable_drug[i])
        y_basal = y_basal+d_lab[i]
    handles, labels = plt.gca().get_legend_handles_labels()
    leg = plt.legend([handles[idx] for idx in Legend_order], [labels[idx] for idx in Legend_order],
                     bbox_to_anchor=Legend_boxanchor, loc='upper center', ncol=Legend_colnum)
    leg.get_frame().set_linewidth(0)
    ax.set_xlabel('Time[days]')
    ax.set_ylabel('Percentage of cell type')
    ax.set_title(title)
    ax.set_xlim([Xlimmin_plot, tmax+Stepsize/4])
    # ax.set_xticks(np.arange(0, math.ceil(t.max()/Stepsize))*Stepsize)
    ax.set_xticks(np.arange(0, math.ceil(tmax/Stepsize))*Stepsize)
    ax.grid()
    ax.set_ylim([ylimmin, ylimmax])
    plt.tight_layout(pad=0.4, w_pad=0.4, h_pad=1.0)
    fig.set_size_inches(13, 9)
    return


def DPM_reversiblemodel_plot_multi(N, ind, name, title):
    ##
    maxt = 0
    logbase = 10

    color_sel = ['k', 'g', 'b', 'r', 'm']
    plt.rcParams.update({'font.size': 20})
    fig, ax = plt.subplots()
    ax.set_yscale('log', base=logbase)
    for i in range(len(N)):
        if type(ind) == int:
            i_N = N[i][ind, :]
        else:
            i_N = N[i][ind, :].sum(axis=0)
        if i_N.ndim == 2:
            i_t = np.arange(0, i_N.shape[0])
        else:
            i_t = np.arange(0, len(i_N))
        ax.plot(i_t, i_N, color=color_sel[i], label=name[i])
        maxt = i_N.shape[0] if i_N.shape[0] > maxt else maxt

    ymaxval = math.ceil(math.log(Limit_mortality, logbase)) + 1 + Ylim_add_plot
    ylimmin = Ylimmin_plot
    if ymaxval > math.log(Limit_mortality, logbase):
        plt.axhline(y=Limit_mortality, color='k', linestyle='dashed')
        plt.text(np.max(maxt)*2/3, 2 * Limit_mortality, 'Limit of mortality', fontsize=13)
    if ylimmin < 1:
        plt.axhline(y=1, color='k', linestyle='dashed')
        plt.text(np.max(maxt)*2/3, 2 * 1, 'Limit of proliferation', fontsize=13)
    leg = plt.legend(loc='best')
    leg.get_frame().set_linewidth(0)
    ax.set_xlabel('Time[days]')
    ax.set_ylabel('Cell Number')

    ax.set_title(title)
    ax.set_xlim([Xlimmin_plot, maxt + Stepsize / 4])
    xticks = np.arange(0, math.ceil(maxt)/Stepsize)*Stepsize
    ax.set_xticks(xticks[::1])
    ax.grid()
    ax.set_ylim([ylimmin, logbase ** ymaxval])
    ylocmaj = mpl.ticker.LogLocator(base=logbase, numticks=len(range(-2, int(ymaxval), 2)))
    ax.yaxis.set_major_locator(ylocmaj)
    ylocmin = mpl.ticker.LogLocator(base=logbase, subs=np.arange(0, logbase) * 1 / logbase, numticks=100)
    ax.yaxis.set_minor_locator(ylocmin)
    ax.yaxis.set_minor_formatter(mpl.ticker.NullFormatter())
    plt.tight_layout(pad=0.25, w_pad=0.25, h_pad=1.0)
    fig.set_size_inches([14.5,  8.5])
    ##
    return


def DPM_reversiblemodel_plot_mut(t, n, title, leg_plot):
    n_ratio = np.zeros(n.shape)
    n_ratio[4, :] = n[4, :]/np.sum(n, 0)
    n_ratio[5, :] = n[5, :]/np.sum(n, 0)
    n_ratio[6, :] = n[6, :]/np.sum(n, 0)
    n_ratio[7, :] = n[7, :]/np.sum(n, 0)
    n_ratio[8, :] = n[8, :]/np.sum(n, 0)

    fig, (ax_n, ax_ratio) = plt.subplots(1, 2)

    ax_n.plot(t, n[4, :], color='g', lw=2, label='0111')
    ax_n.plot(t, n[5, :], color='k', lw=2, label='0011')
    ax_n.plot(t, n[6, :], color='b', lw=2, label='1101')
    ax_n.plot(t, n[7, :], color='r', lw=2, label='1100')
    ax_n.plot(t, n[8, :], color='c', lw=2, label='1111')

    if leg_plot:
        leg = ax_n.legend(loc='best', ncol=1, frameon=1, fontsize=19)
        leg.draw_frame(False)
    ax_n.set_xlabel('Time[days]')
    ax_n.set_ylabel('Cell Number')
    if np.max(n[4:, :]) > 0.:
        ax_n.set_yscale('log', base=10)

    ax_n.set_title(title)

    ax_ratio.plot(t, n_ratio[4, :], color='g', lw=2, label='0111')
    ax_ratio.plot(t, n_ratio[5, :], color='k', lw=2, label='0011')
    ax_ratio.plot(t, n_ratio[6, :], color='b', lw=2, label='1101')
    ax_ratio.plot(t, n_ratio[7, :], color='r', lw=2, label='1100')
    ax_ratio.plot(t, n_ratio[8, :], color='c', lw=2, label='1111')
    ax_ratio.set_xlabel('Time[days]')
    ax_ratio.set_ylabel('Cell subtype percentage')
    ax_ratio.set_title(title)

    fig.set_size_inches(17,  8)
    plt.tight_layout(pad=0.4, w_pad=0.4, h_pad=1.0)
    return


def DPM_reversiblemodel_plot_km_multi(km, par, titlestr):
    ##
    plt.rcParams['font.size'] = 14
    plt.figure()
    fig = plt.gcf()
    fig.set_size_inches(18, 10)
    ax = plt.axes()
    keys = km.keys()
    legh = []
    legstr = []
    color_ = ['#DAA520', 'k', 'b', 'r', 'g', '#FF82AB', '#FF82AB', 'blueviolet', 'blueviolet', 'c', 'c', 'm', 'm', '#808080', '#808080']
    linestyle = ['-', '-', '-', '-', '-', '-', '--', '-', '--', '-', '--', '-', '--', '-', '--', '-', '--']
    for i, i_key in enumerate(keys):
        i_km = km[i_key]
        l_1i, = plt.plot(i_km['t'], i_km['val'], color=color_[i], linestyle=linestyle[i])
        ax.fill_between(i_km['t'], i_km['interval_lower'], i_km['interval_upper'], color=color_[i], alpha=.3)
        l_2i, = plt.plot([], [], ' ')
        legh.extend([l_1i, l_2i])
        if i_km['median_survival'] == np.inf:
            i_km['median_survival'] = duration_5year
        legstr.extend([i_key, f"median survival time of {i_key} (days): {int(i_km['median_survival'])}"])

    plt.title(f"{int(par['totalnum'])} of patients")
    plt.xticks(list(range(0, par['duration']+par['xtick step'], par['xtick step'])))
    plt.yticks(np.arange(0, 1.2, .2))
    plt.xlabel('Survival time[days]')
    plt.ylabel('Fraction of surviving patients')
    plt.legend(legh, legstr, loc='upper right', frameon=False, ncol=2)
    plt.tight_layout()
    ##
    return


def DPM_reversiblemodel_plot_hist(fig, ax, data, bins=None, normed_area=False, normed_height=True, xticks=None, xticklabel=None,
                                  xlabel=None, ylim=None, title=None, legend=False, legend_loc='best', legend_fontsize='medium',
                                  legend_labels=None, facecolor=None, edgecolor=None, xlog=False, label=''):
    # Plot one 1D histogram
    # Default colors
    facecolor = facecolor if facecolor is not None else 'b'
    edgecolor = None if edgecolor is None else None
    xticklable = [str(i) for i in xticks] if xticklabel is None else xticklabel

    hist_kwargs = dict()
    hist_kwargs['x'] = data
    hist_kwargs['bins'] = bins
    hist_kwargs['histtype'] = 'bar'
    hist_kwargs['facecolor'] = facecolor
    hist_kwargs['edgecolor'] = edgecolor
    hist_kwargs['alpha'] = 0.6
    hist_kwargs['rwidth'] = 1
    hist_kwargs['align'] = 'right'

    # Calculate weights if normalizing bins by height
    if normed_height and not normed_area:
        hist_kwargs['weights'] = np.ones_like(hist_kwargs['x'], dtype=float)
        hist_kwargs['weights'] /= float(len(hist_kwargs['x']))
        hist_kwargs['weights'] *= 100

    # fig, ax = plt.subplots()
    fig.set_size_inches(14, 9)
    counts, bins, _ = ax.hist(**hist_kwargs, label=label)
    xlog and ax.set_xscale('log', base=10)
    ax.set_xticks(xticks)
    ax.get_xaxis().set_major_formatter(mpl.ticker.ScalarFormatter())
    ax.set_xticklabels(xticklable)
    ax.set_xlabel(xlabel)
    ax.set_ylabel('Percentage')
    ax.set_xlim((bins[0], bins[-1]))
    ax.set_ylim((0, 105))
    # plt.tight_layout()

    ylim is not None and ax.set_ylim(ylim)
    title is not None and ax.set_title(title)
    legend and legend_labels is not None and ax.set_legend(legend_labels, loc=legend_loc, prop={'size': legend_fontsize})
    return counts, bins


def DPM_reversiblemodel_plot_df2pdf(df, filename, color_, numpages=(1, 1), pagesize=(11, 8.5)):
    def DPM_reversiblemodel_plot_df2pdf_1(df_, pagesize_):
        # draw as table
        alternating_colors = [color_[0] * len(df_.columns)]
        alternating_colors = alternating_colors * len(df_)
        # plt.rc('font', size=30)
        fig_, ax = plt.subplots(figsize=pagesize)
        ax.axis('tight')
        ax.axis('off')
        the_table = ax.table(cellText=df_.values,
                             rowLabels=None,
                             colLabels=df_.columns,
                             rowColours=color_[1]*len(df_),
                             colColours=color_[1]*len(df_.columns),
                             cellColours=alternating_colors,
                             loc='center')
        the_table.auto_set_font_size(False)
        the_table.set_fontsize(7)   # 8
        the_table.scale(1, 2)       # 4
        return fig_

    with PdfPages(filename) as pdf:
        nh, nv = numpages
        rows_per_page = len(df)//nh
        cols_per_page = len(df.columns)//nv
        for i in range(0, nh):
            for j in range(0, nv):
                page = df.iloc[(i*rows_per_page):min((i+1)*rows_per_page, len(df)),
                               (j*cols_per_page):min((j+1)*cols_per_page, len(df.columns))]
                fig = DPM_reversiblemodel_plot_df2pdf_1(page, pagesize)
                if nh > 1 or nv > 1:
                    # Add part/page number at bottom-center of page
                    fig.text(0.5, 0.5/pagesize[0],
                             'Part-{}x{}: Page-{}'.format(i+1, j+1, i*nv + j + 1), ha='center', fontsize=9)
                pdf.savefig(fig, bbox_inches='tight')
                plt.close()
    return


def DPM_reversiblemodel_plot_linearreg(x, data, model, title, ylim, alpha=0.05, stepsize=0.1):
    x = x.flatten()
    # xx = np.arange(0, 10, 1)/10
    # pred = model.get_prediction(xx).summary_frame(alpha)

    plt.rcParams['font.size'] = 16
    fig, ax = plt.subplots()
    fig.set_size_inches(12, 9)
    l1 = ax.scatter(x, data, c='b')
    # l2, = ax.plot(x, list(model.fittedvalues), c='b')
    # l3 = ax.fill_between(xx, pred['mean_ci_lower'], pred['mean_ci_upper'], alpha=0.5)
    l4, = ax.plot([], [], ' ')
    l5, = ax.plot([], [], ' ')
    plt.title(title)
    # plt.xticks(xx)
    plt.xlabel('values')
    plt.ylabel('Survival time [days]')
    plt.ylim(ylim)
    # plt.legend([l1, l2, l3, l4, l5], ['exp', 'linear regression', f'{(1 - alpha)*100}% confidence interval',
    #                                   f'$r^2$={round(model.rsquared, 3)}, slope={round(model.params[0], 3)}'],
    #            loc='upper left', frameon=False)
    return


def DPM_reversiblemodel_plot_contour(data_ref, data, par):
    duration, binsize, name = par['duration'], par['binsize'], par['name']

    data_ref = [i_val if i_val <= duration else duration for i_val in data_ref]
    data = [i_val if i_val <= duration else duration for i_val in data]
    h, x, y = np.histogram2d(data_ref, data, bins=[np.arange(0, duration + binsize, binsize), np.arange(0, duration + binsize, binsize)],
                             range=[[0, duration + binsize], [0, duration + binsize]])
    h = h.T
    # h = h/np.sum(h)
    cmap = plt.colormaps['bwr']
    levels = mpl.ticker.MaxNLocator(nbins=math.ceil(np.log10(h.max()))).tick_values(0, np.log10(h.max()))
    norm = mpl.colors.BoundaryNorm(levels, ncolors=cmap.N, clip=True)

    # levels = np.array([2e-5, 1e-4, 1e-3, 1e-2, 1e-1])
    xc = (x[:-1] + x[1:]) / 2.0  # x-axis bin centers
    yc = (y[:-1] + y[1:]) / 2.0
    X, Y = np.meshgrid(xc, yc)
    plt.rcParams['font.size'] = 14
    plt.rcParams['figure.figsize'] = (13, 9)
    # plt.figure()
    # fig = plt.gcf()
    # ax = plt.axes()
    # cp = ax.contour(X, Y, h, levels=levels, cmap='bwr')  # , colors='black'
    # ax.set_facecolor('white')
    im = plt.pcolormesh(X, Y, np.log10(h, out=np.full_like(h, fill_value=np.nan), where=(h != 0)),
                        shading='nearest', edgecolors=None, linewidth=0.5, norm=norm, cmap=cmap)
    cbar = plt.colorbar(im)
    plt.xlim([0, duration])
    plt.ylim([0, duration])
    plt.xticks(list(range(0, duration + binsize, par['xtick step'])))
    plt.yticks(list(range(0, duration + binsize, par['xtick step'])))
    plt.title(f'{int(np.sum(h))} of patients')
    plt.xlabel(name[0] + ' [days]')
    plt.ylabel(name[1] + ' [days]')
    cbar.set_label('log10(# of patientss) ', rotation=270, labelpad=30)
    # plt.legend(loc='upper center', bbox_to_anchor=(0.5, 1.08), frameon=False, ncol=2)
    plt.tight_layout(pad=0.4, w_pad=0.4, h_pad=1.0)
    return


def DPM_reversiblemodel_plot_1strategy(N, d, tmax, title, pathsave, mutonly=False, reveronly=False):
    t = np.arange(0, N.shape[1])
    if mutonly or reveronly:
        pathsave_mutonly = os.path.join(pathsave, title.split(':')[1]+'.pdf')
        DPM_reversiblemodel_plot_cellnumwithd_mutonly_reveronly(t, N, d, tmax, title, mutonly, reveronly)
        plt.savefig(pathsave_mutonly, format='pdf', bbox_inches='tight')
        plt.close('all')
    else:
        # nonmut cell
        pathsave_nonmut = os.path.join(pathsave, title.split(':')[1]+' nonmut.pdf')
        DPM_reversiblemodel_plot_cellnumwithd(t, N, d, tmax, title, Ntype='nonmut')
        plt.savefig(pathsave_nonmut, format='pdf', bbox_inches='tight')
        plt.close('all')

        # pathsave_nonmut_ratio = os.path.join(pathsave, title.split(':')[1] + ' nonmut ratio.pdf')
        # DPM_reversiblemodel_plot_cellnumpercentagewithd(t, N, d, tmax, title, Ntype='nonmut')
        # plt.savefig(pathsave_nonmut_ratio, format='pdf', bbox_inches='tight')
        # plt.close('all')

        # mut cell
        pathsave_mut = os.path.join(pathsave, title.split(':')[1] + ' mut.pdf')
        DPM_reversiblemodel_plot_cellnumwithd(t, N, d, tmax, title, Ntype='mut')
        plt.savefig(pathsave_mut, format='pdf', bbox_inches='tight')
        plt.close('all')

        # pathsave_mut_ratio = os.path.join(pathsave, title.split(':')[1] + ' mut ratio.pdf')
        # DPM_reversiblemodel_plot_cellnumpercentagewithd(t, N, d, tmax, title, Ntype='mut')
        # plt.savefig(pathsave_mut_ratio, format='pdf', bbox_inches='tight')
        # plt.close('all')
    return


def DPM_reversiblemodel_plot_stackplot(data, title, legstr, ylim):
    plt.rcParams.update({'font.size': 15})
    fig, ax = plt.subplots()
    fig.set_size_inches(14, 8)
    p = plt.stackplot(np.arange(0, data.shape[1]), data, colors=color33)

    ax.set_xticks(range(0, duration_5year, 200))
    ax.set_ylabel('Percentage')
    ax.set_xlabel('[days]')
    ax.set_xlim((0, duration_5year))
    ax.set_ylim(ylim)
    ax.set_title(title)

    if data.shape[0] == 9:
        ncol = 5
        loc = 'upper left'
        bbox_to_anchor = (0, 1)
    elif data.shape[0] == 12:
        ncol = 4
        loc = 'upper left'
        bbox_to_anchor = (0, 1)
    else:
        ncol = None
        lco = None
        bbox_to_anchor = None

    leg = plt.legend(p, legstr, bbox_to_anchor=bbox_to_anchor, loc=loc, ncol=ncol, frameon=False)
    plt.tight_layout()
    return


def DPM_reversiblemodel_plot_sankeydiagram(sankey_df, category, name, pathsave_figure_sankeydiagram):
    def DPM_reversiblemodel_plot_sankeydiagram_1(sankey_df_, idx, name_):
        if isinstance(idx, str):
            new_val = sankey_df_.loc[idx]
        elif isinstance(idx, list):
            new_val = sankey_df_.loc[idx].sum(axis=0)
        else:
            raise Exception('idx type error')
        new_val.name = name_
        sankey_df_ = pd.concat([sankey_df_, pd.DataFrame(new_val).T], axis=0)
        return sankey_df_

    color_ = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#c3b530']
    source_name, color_all = None, None
    target_name = ['cyc', 's1', 's2', 's1,s2', 'other']
    target_comb = ["".join(seq) for seq in itertools.product("01", repeat=len(sankey_df.columns[0]))]
    target_other = list(set(target_comb).difference({'1000', '0010', '0001', '0011'}))
    sankey_df['cyc'] = sankey_df['1000']
    sankey_df['s1'] = sankey_df['0010']
    sankey_df['s2'] = sankey_df['0001']
    sankey_df['s1,s2'] = sankey_df['0011']
    sankey_df['other'] = sankey_df[target_other].sum(axis=1)
    sankey_df = sankey_df[target_name]
    if name == 'mut to full':
        source_name = ['s1', 's2', 's1,s2', 'cyc,s1,s2', 'cyc,s0,s1,s2', 'other']
        source_comb = ["".join(seq) for seq in itertools.product("01", repeat=len(sankey_df.index[0]))]
        source_other = list(set(source_comb).difference({'0010', '0001', '0011', '1011', '1111'}))

        sankey_df = DPM_reversiblemodel_plot_sankeydiagram_1(sankey_df, '0010', 's1')
        sankey_df = DPM_reversiblemodel_plot_sankeydiagram_1(sankey_df, '0001', 's2')
        sankey_df = DPM_reversiblemodel_plot_sankeydiagram_1(sankey_df, '0011', 's1,s2')
        sankey_df = DPM_reversiblemodel_plot_sankeydiagram_1(sankey_df, '1011', 'cyc,s1,s2')
        sankey_df = DPM_reversiblemodel_plot_sankeydiagram_1(sankey_df, '1111', 'cyc,s0,s1,s2')
        sankey_df = DPM_reversiblemodel_plot_sankeydiagram_1(sankey_df, source_other, 'other')
        sankey_df = sankey_df.loc[source_name]
        category_comb = category['category']
        category = list()
        for i in category_comb:
            if i == '0010':
                category.append('s1')
            elif i == '0001':
                category.append('s2')
            elif i == '0011':
                category.append('s1,s2')
            elif i == '1011':
                category.append('cyc,s1,s2')
            elif i == '1111':
                category.append('cyc,s0,s1,s2')
            else:
                category.append('other')
        color_all = color_ + color_[:-1]
    elif name == 'rever to full':
        source_name = ['cyc', 's1', 'cyc,s1', 'cyc,s0,s1']
        source_comb = ["".join(seq) for seq in itertools.product("01", repeat=len(sankey_df.index[0]))]

        sankey_df = DPM_reversiblemodel_plot_sankeydiagram_1(sankey_df, '100', 'cyc')
        sankey_df = DPM_reversiblemodel_plot_sankeydiagram_1(sankey_df, '001', 's1')
        sankey_df = DPM_reversiblemodel_plot_sankeydiagram_1(sankey_df, '101', 'cyc,s1')
        sankey_df = DPM_reversiblemodel_plot_sankeydiagram_1(sankey_df, '111', 'cyc,s0,s1')
        sankey_df = sankey_df.loc[source_name]
        category_comb = category['category']
        category = list()
        for i in category_comb:
            if i == '100':
                category.append('cyc')
            elif i == '001':
                category.append('s1')
            elif i == '101':
                category.append('cyc,s1')
            elif i == '111':
                category.append('cyc,s0,s1')
        color_all = color_[:-2] + color_[:-1]

    assert sankey_df.shape[0] == len(source_name) and sankey_df.shape[1] == len(target_name)
    assert len(color_all) == len(source_name) + len(target_name)
    x = [0.2] * len(source_name) + [0.8] * len(target_name)

    num_par = sankey_df.values.sum()
    num_partarget = sankey_df.sum(axis=0).values
    data = list()
    source, target, value, color_link, val = [], [], [], [], 0
    for i, i_name in enumerate(source_name):
        i_source = [source_name[i] + f', {round(category.count(i_name)/len(category) * 100, 2)}%']
        i_total = sankey_df.sum(axis=1)[i]
        i_valrow = 0
        for j in range(len(target_name)):
            j_target = sankey_df.iloc[i, j]
            j_val_row = j_target/i_total * 100
            j_val_col = j_target/num_partarget[j] * 100
            i_source.append(f'{j_target}\nrow {round(j_val_row, 2)}%\ncol {round(j_val_col, 2)}%')
            source.append(i)
            target.append(j + len(source_name))
            color_link.append(color_[i])
            value.append(j_target)
            if not math.isnan(j_val_row):
                i_valrow = i_valrow + j_val_row
        assert isclose(i_valrow, 100, abs_tol=1e-9)
        i_source.append(f'{i_total}, {round(i_total/num_par * 100, 2)}%')
        val = val + i_total/num_par
        data.append(i_source)
    assert isclose(val, 1, abs_tol=1e-9)
    lastline = ['']
    [lastline.append(f'{i}, {round(i/num_par * 100, 2)}%') for i in list(sankey_df.values.sum(axis=0))]
    lastline.append('')
    data.append(lastline)
    col_names = [f'No. of patients\n: {num_par}']
    col_names.extend(target_name)
    col_names.append('')
    print(tabulate(data, headers=col_names, tablefmt='rst', numalign='center'))
    data_df = pd.DataFrame(data, columns=col_names)
    filename_table = os.path.join(pathsave_figure_sankeydiagram, f'{name} num_sel.pdf')
    DPM_reversiblemodel_plot_df2pdf(data_df, filename_table, (['white'], ['lightgray']))

    y = []
    [y.append(int(i.split(',')[0])) for i in list(data_df.iloc[:, -1].values)[:-1]]
    y = list(np.cumsum(np.array(y)) / num_par)
    y_ = [0.001] + y[:-1]

    y = []
    [y.append(int(i.split(',')[0])) for i in list(data_df.iloc[-1, :].values[1:-1])]
    y = list(np.cumsum(np.array(y)) / num_par)
    y_ = y_ + [0.001] + y[:-1]

    fig = go.Figure(data=[go.Sankey(
        arrangement='snap',
        node=dict(
            pad=15,
            thickness=15,
            line=dict(color="black", width=0.5),
            label=source_name + target_name,
            x=x,
            y=y_,
            color=color_all
        ),
        link=dict(
            source=source,
            target=target,
            value=value,
            color=color_link
        ))])
    fig.update_layout(title_text="Sankey Diagram", font_size=13)
    fig.show()
    fig.update_layout(width=1920, height=1080)
    filename_sankeydiagram = os.path.join(pathsave_figure_sankeydiagram, f'{name}.pdf')
    fig.write_image(filename_sankeydiagram)
    return
