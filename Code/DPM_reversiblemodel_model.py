from DPM_reversiblemodel_constant import *


def DPM_reversiblemodel_model_transitionrate(α_i, θ_i, κ_i, μ_i):
    Storr_σ0 = μ_i * np.exp(-κ_i * (DPM_reversiblemodel_model_Eplus(σ0, θ_i) -
                                    DPM_reversiblemodel_model_Eminus(σ0,  θ_i, α_i)))
    Storr_σhalf = μ_i * np.exp(-κ_i * (DPM_reversiblemodel_model_Eplus(σhalf, θ_i) -
                                       DPM_reversiblemodel_model_Eminus(σhalf, θ_i, α_i)))
    Storr_σfull = μ_i * np.exp(-κ_i * (DPM_reversiblemodel_model_Eplus(σfull, θ_i) -
                                       DPM_reversiblemodel_model_Eminus(σfull, θ_i, α_i)))

    rrtoS_σ0 = μ_i * np.exp(-κ_i * -(DPM_reversiblemodel_model_Eplus(σ0, θ_i) -
                                     DPM_reversiblemodel_model_Eminus(σ0, θ_i, α_i)))
    rrtoS_σhalf = μ_i * np.exp(-κ_i * -(DPM_reversiblemodel_model_Eplus(σfull, θ_i) -
                                        DPM_reversiblemodel_model_Eminus(σfull, θ_i, α_i)))
    rrtoS_σfull = μ_i * np.exp(-κ_i * -(DPM_reversiblemodel_model_Eplus(σfull, θ_i) -
                                        DPM_reversiblemodel_model_Eminus(σfull, θ_i, α_i)))

    return Storr_σ0, Storr_σhalf, Storr_σfull, rrtoS_σ0, rrtoS_σhalf, rrtoS_σfull


def DPM_reversiblemodel_model_prorate(g, S_01_1, S_01_2, S_00_1, S_00_2, S_11_1, S_11_2):
    S_σ1_full = g - S_01_1
    S_σ2_full = g - S_01_2
    R1_σ1_full = g - S_00_1
    R2_σ2_full = g - S_00_2
    IR1_σ1_full = g - S_11_1
    IR2_σ2_full = g - S_11_2

    return dict(zip(['S σ1 full', 'S σ2 full', 'R1 σ1 full', 'R2 σ2 full', 'IR1 σ1 full', 'IR2 σ2 full'],
                    [S_σ1_full, S_σ2_full, R1_σ1_full, R2_σ2_full, IR1_σ1_full, IR2_σ2_full]))


def DPM_reversiblemodel_model_simexpm(argsODE, t, x0, LSsim=True):
    x0, N = np.array(x0), np.array([x0]).T
    A = DPM_reversiblemodel_model_matrix(*argsODE)
    flag_mortality, flag_cure = False, False
    for i_t in t[1:]:
        x0[x0 < 1] = 0
        i_N = expm(A * i_t) @ x0
        if sum(i_N) >= Limit_mortality and LSsim:
            flag_mortality = True
            break
        N = np.append(N, np.reshape(i_N, (-1, 1)), 1)
        if all(i_N < 1):
            flag_cure = True
            break
        x0 = i_N

    t = t[:N.shape[1]]
    t = np.cumsum(t)
    N_total = np.sum(N, axis=0)
    N_nonmut = np.sum(N[ind_nonmut, :], axis=0)
    N_mut = np.sum(N[ind_mut, :], axis=0)
    return N, t, N_total, N_nonmut, N_mut, flag_mortality, flag_cure


def DPM_reversiblemodel_model_Sa(par):
    Sa = np.array([[par['S_01_1'], par['S_01_2']],
                   [par['S_00_1'], par['S_01_2']],
                   [par['S_01_1'], par['S_00_2']],
                   [par['S_00_1'], par['S_00_2']],
                   [par['S_01_1'], par['S_11_2']],
                   [par['S_00_1'], par['S_11_2']],
                   [par['S_11_1'], par['S_01_2']],
                   [par['S_11_1'], par['S_00_2']],
                   [par['S_11_1'], par['S_11_2']]
                   ])
    return Sa


# Heaviside function
def DPM_reversiblemodel_model_Heaviside(x):
    return 1 * (x >= 0)


# Reduction and translocation factors to the expression
def DPM_reversiblemodel_model_A(σ, θ):
    return 1-σ*(1-θ)+1e-6


# Production function as a step-like function
def DPM_reversiblemodel_model_f(y, σ, θ, α):
    return DPM_reversiblemodel_model_A(σ, θ)*(α+(1-α)*DPM_reversiblemodel_model_Heaviside(y-θ))


# Potential function
def DPM_reversiblemodel_model_U(y, σ, θ, α):
    return -DPM_reversiblemodel_model_A(σ, θ)*(α+(1-α)*DPM_reversiblemodel_model_Heaviside(y-θ))*(y-θ)+(y**2-θ**2)/2


# E+, value unrelated to α, so set α=-1
def DPM_reversiblemodel_model_Eplus(σ, θ):
    return -DPM_reversiblemodel_model_U(DPM_reversiblemodel_model_A(σ, θ), σ, θ, -1)


# E-
def DPM_reversiblemodel_model_Eminus(σ, θ, α):
    return -DPM_reversiblemodel_model_U(α*DPM_reversiblemodel_model_A(σ, θ), σ, θ, α)


def DPM_reversiblemodel_model_ODE(t, x, σ, mutation, argsODE):
    σ_t = np.arange(1, σ.shape[1]+Simtimestep, Simtimestep)
    pos = np.where(t <= σ_t)[0][0]
    # print(t)

    σ_pos = σ[:, pos]
    # len(pos) == 0:
    # σ_pos = σ[:, -1]

    A = DPM_reversiblemodel_model_matrix(σ_pos, mutation, argsODE)
    return A @ x


def DPM_reversiblemodel_model_matrix(σ, mutation, argsODE):
    α1, θ1, κ1 = argsODE['α1'], argsODE['θ1'], argsODE['κ1']
    α2, θ2, κ2 = argsODE['α2'], argsODE['θ2'], argsODE['κ2']
    μ1, μ2 = argsODE['μ1'], argsODE['μ2']

    g = argsODE['g']
    S_01_1, S_00_1, S_11_1 = argsODE['S_01_1'], argsODE['S_00_1'], argsODE['S_11_1']
    S_01_2, S_00_2, S_11_2 = argsODE['S_01_2'], argsODE['S_00_2'], argsODE['S_11_2']

    T1 = argsODE['T1'] if mutation is True else 0.
    T2 = argsODE['T2'] if mutation is True else 0.

    σ1, σ2 = σ

    ΔE_S_rr1 = -κ1 * (DPM_reversiblemodel_model_Eplus(σ1, θ1) - DPM_reversiblemodel_model_Eminus(σ1, θ1, α1))
    ΔE_S_rr2 = -κ2 * (DPM_reversiblemodel_model_Eplus(σ2, θ2) - DPM_reversiblemodel_model_Eminus(σ2, θ2, α2))

    '''4 states, no mutation happens'''
    dn0101 = [g - μ1 * np.exp(ΔE_S_rr1) - μ2 * np.exp(ΔE_S_rr2) - S_01_1 * σ1 - S_01_2 * σ2,   # n0101
              μ1 * np.exp(-ΔE_S_rr1),  # n0001
              μ2 * np.exp(-ΔE_S_rr2),  # n0100
              0,                       # n0000
              0,                       # n0111
              0,                       # n0011
              0,                       # n1101
              0,                       # n1100
              0]                       # n1111

    dn0001 = [μ1 * np.exp(ΔE_S_rr1),   # n0101
              g - μ1 * np.exp(-ΔE_S_rr1) - μ2 * np.exp(ΔE_S_rr2) - S_00_1 * σ1 - S_01_2 * σ2,  # n0001
              0,                       # n0100
              μ2 * np.exp(-ΔE_S_rr2),  # n0000
              0,                       # n0111
              0,                       # n0011
              0,                       # n1101
              0,                       # n1100
              0]                       # n1111

    dn0100 = [μ2 * np.exp(ΔE_S_rr2),   # n0101
              0,                       # n0001
              g - μ1 * np.exp(ΔE_S_rr1) - μ2 * np.exp(-ΔE_S_rr2) - S_01_1 * σ1 - S_00_2 * σ2,  # n0100
              μ1 * np.exp(-ΔE_S_rr1),  # n0000
              0,                       # n0111
              0,                       # n0011
              0,                       # n1101
              0,                       # n1100
              0]                       # n1111

    dn0000 = [0,                       # n0101
              μ2 * np.exp(ΔE_S_rr2),   # n0001
              μ1 * np.exp(ΔE_S_rr1),   # n0100
              g - μ1 * np.exp(-ΔE_S_rr1) - μ2 * np.exp(-ΔE_S_rr2) - S_00_1 * σ1 - S_00_2 * σ2,  # n0000
              0,                       # n0111
              0,                       # n0011
              0,                       # n1101
              0,                       # n1100
              0]                       # n1111

    '''mutated states'''
    dn0111 = [g * T2,                  # n0101
              0,                       # n0001
              g * T2,                  # n0100
              0,                       # n0000
              g - μ1 * np.exp(ΔE_S_rr1) - S_01_1 * σ1 - S_11_2 * σ2,   # n0111
              μ1 * np.exp(-ΔE_S_rr1),  # n0011
              0,                       # n1101
              0,                       # n1100
              0]                       # n1111

    dn0011 = [0,                       # n0101
              g * T2,                  # n0001
              0,                       # n0100
              g * T2,                  # n0000
              μ1 * np.exp(ΔE_S_rr1),   # n0111
              g - μ1 * np.exp(-ΔE_S_rr1) - S_00_1 * σ1 - S_11_2 * σ2,  # n0011
              0,                       # n1101
              0,                       # n1100
              0]                       # n1111

    dn1101 = [g * T1,                  # n0101
              g * T1,                  # n0001
              0,                       # n0100
              0,                       # n0000
              0,                       # n0111
              0,                       # n0011
              g - μ2 * np.exp(ΔE_S_rr2) - S_11_1 * σ1 - S_01_2 * σ2,  # n1101
              μ2 * np.exp(-ΔE_S_rr2),  # n1100
              0]                       # n1111

    dn1100 = [0,                       # n0101
              0,                       # n0001
              g * T1,                  # n0100
              g * T1,                  # n0000
              0,                       # n0111
              0,                       # n0011
              μ2 * np.exp(ΔE_S_rr2),   # n1101
              g - μ2 * np.exp(-ΔE_S_rr2) - S_11_1 * σ1 - S_00_2 * σ2,   # n1100
              0]                       # n1111

    dn1111 = [0,                      # n0101
              0,                      # n0001
              0,                      # n0100
              0,                      # n0000
              g * T1,                 # n0111
              g * T1,                 # n0011
              g * T2,                 # n1101
              g * T2,                 # n1100
              g - S_11_1 * σ1 - S_11_2 * σ2]  # n1111
    # f = np.array([dn0101, dn0001, dn0100, dn0000, dn0111, dn0011, dn1101, dn1100, dn1111])
    return np.array([dn0101, dn0001, dn0100, dn0000, dn0111, dn0011, dn1101, dn1100, dn1111])


def DPM_reversiblemodel_fun_model_flow(par, N, d):
    d_σ1 = np.zeros([3, d.shape[1]])
    d_σ1[0, d[0, :] == 0] = 1
    d_σ1[1, d[1, :] == 0.5] = 1
    d_σ1[2, d[0, :] == 1] = 1
    assert np.all(np.sum(d_σ1, axis=0) == 1)

    d_σ2 = np.zeros([3, d.shape[1]])
    d_σ2[0, d[1, :] == 0] = 1
    d_σ2[1, d[1, :] == 0.5] = 1
    d_σ2[2, d[1, :] == 1] = 1
    assert np.all(np.sum(d_σ1, axis=0) == 1)

    σ1, σ2 = d

    α1, θ1, κ1 = par['α1'], par['θ1'], par['κ1']
    α2, θ2, κ2 = par['α2'], par['θ2'], par['κ2']
    μ1, μ2 = par['μ1'], par['μ2']

    ΔE_S_rr1_σ1_0 = -κ1 * (DPM_reversiblemodel_model_Eplus(σ0, θ1) - DPM_reversiblemodel_model_Eminus(σ0, θ1, α1))
    ΔE_S_rr1_σ1_half = -κ1 * (DPM_reversiblemodel_model_Eplus(σhalf, θ1) - DPM_reversiblemodel_model_Eminus(σhalf, θ1, α1))
    ΔE_S_rr1_σ1_full = -κ1 * (DPM_reversiblemodel_model_Eplus(σfull, θ1) - DPM_reversiblemodel_model_Eminus(σfull, θ1, α1))

    ΔE_S_rr2_σ2_0 = -κ2 * (DPM_reversiblemodel_model_Eplus(σ0, θ2) - DPM_reversiblemodel_model_Eminus(σ0, θ2, α2))
    ΔE_S_rr2_σ2_half = -κ2 * (DPM_reversiblemodel_model_Eplus(σhalf, θ2) - DPM_reversiblemodel_model_Eminus(σhalf, θ2, α2))
    ΔE_S_rr2_σ2_full = -κ2 * (DPM_reversiblemodel_model_Eplus(σfull, θ2) - DPM_reversiblemodel_model_Eminus(σfull, θ2, α2))

    pos_ΔE_S_rr1_σ1_0, neg_ΔE_S_rr1_σ1_0 = μ1 * np.exp(ΔE_S_rr1_σ1_0),  μ1 * np.exp(-ΔE_S_rr1_σ1_0)
    pos_ΔE_S_rr1_σ1_half, neg_ΔE_S_rr1_σ1_half = μ1 * np.exp(ΔE_S_rr1_σ1_half),  μ1 * np.exp(-ΔE_S_rr1_σ1_half)
    pos_ΔE_S_rr1_σ1_full, neg_ΔE_S_rr1_σ1_full = μ1 * np.exp(ΔE_S_rr1_σ1_full), μ1 * np.exp(-ΔE_S_rr1_σ1_full)

    pos_ΔE_S_rr2_σ2_0, neg_ΔE_S_rr2_σ2_0 = μ2 * np.exp(ΔE_S_rr2_σ2_0), μ2 * np.exp(-ΔE_S_rr2_σ2_0)
    pos_ΔE_S_rr2_σ2_half, neg_ΔE_S_rr2_σ2_half = μ2 * np.exp(ΔE_S_rr2_σ2_half), μ2 * np.exp(-ΔE_S_rr2_σ2_half)
    pos_ΔE_S_rr2_σ2_full, neg_ΔE_S_rr2_σ2_full = μ2 * np.exp(ΔE_S_rr2_σ2_full), μ2 * np.exp(-ΔE_S_rr2_σ2_full)

    pos_ΔE_σ1 = np.array([[pos_ΔE_S_rr1_σ1_0, pos_ΔE_S_rr1_σ1_half, pos_ΔE_S_rr1_σ1_full]])
    neg_ΔE_σ1 = np.array([[neg_ΔE_S_rr1_σ1_0, neg_ΔE_S_rr1_σ1_half, neg_ΔE_S_rr1_σ1_full]])

    pos_ΔE_σ2 = np.array([[pos_ΔE_S_rr2_σ2_0, pos_ΔE_S_rr2_σ2_half, pos_ΔE_S_rr2_σ2_full]])
    neg_ΔE_σ2 = np.array([[neg_ΔE_S_rr2_σ2_0, neg_ΔE_S_rr2_σ2_half, neg_ΔE_S_rr2_σ2_full]])

    pos_ΔE_σ1, neg_ΔE_σ1 = np.matmul(pos_ΔE_σ1, d_σ1).flatten(), np.matmul(neg_ΔE_σ1, d_σ1).flatten()
    pos_ΔE_σ2, neg_ΔE_σ2 = np.matmul(pos_ΔE_σ2, d_σ2).flatten(), np.matmul(neg_ΔE_σ2, d_σ2).flatten()

    g = par['g']
    S_01_1, S_00_1, S_11_1 = par['S_01_1'], par['S_00_1'], par['S_11_1']
    S_01_2, S_00_2, S_11_2 = par['S_01_2'], par['S_00_2'], par['S_11_2']

    T1, T2 = par['T1'], par['T2']

    n0101, n0001, n0100, n0000, n0111, n0011, n1101, n1100, n1111 = N[:, :-1]

    n0101_pro = np.absolute((g - S_01_1 * g * σ1 - S_01_1 * g * σ1) * n0101)  # 0
    n0001_pro = np.absolute((g - S_00_1 * g * σ1 - S_01_2 * g * σ2) * n0001)  # 1
    n0100_pro = np.absolute((g - S_01_1 * g * σ1 - S_00_2 * g * σ2) * n0100)  # 2
    n0000_pro = np.absolute((g - S_00_1 * g * σ1 - S_00_2 * g * σ2) * n0000)  # 3
    n0111_pro = np.absolute((g - S_01_1 * g * σ1 - S_11_2 * g * σ2) * n0111)  # 4
    n0011_pro = np.absolute((g - S_00_1 * g * σ1 - S_11_2 * g * σ2) * n0011)  # 5
    n1101_pro = np.absolute((g - S_11_1 * g * σ1 - S_01_2 * g * σ2) * n1101)  # 6
    n1100_pro = np.absolute((g - S_11_1 * g * σ1 - S_00_2 * g * σ2) * n1100)  # 7
    n1111_pro = np.absolute((g - S_11_1 * g * σ1 - S_11_2 * g * σ2) * n1111)  # 8

    n0101_to_n0001 = pos_ΔE_σ1 * n0101  # 9
    n0101_to_n0100 = pos_ΔE_σ2 * n0101  # 10
    n0001_to_n0101 = neg_ΔE_σ1 * n0001  # 11
    n0100_to_n0101 = neg_ΔE_σ2 * n0100  # 12
    n0001_to_n0000 = pos_ΔE_σ2 * n0001  # 13
    n0100_to_n0000 = pos_ΔE_σ1 * n0100  # 14
    n0000_to_n0001 = neg_ΔE_σ1 * n0000  # 15
    n0000_to_n0100 = neg_ΔE_σ2 * n0000  # 16
    n1101_to_n1100 = pos_ΔE_σ2 * n1101  # 17
    n1100_to_n1101 = neg_ΔE_σ2 * n1100  # 18
    n0111_to_n0011 = pos_ΔE_σ1 * n0111  # 19
    n0011_to_n0111 = neg_ΔE_σ1 * n0011  # 20

    n0101_to_n1101 = g * T1 * n0101  # 21
    n0101_to_n0111 = g * T2 * n0101  # 22
    n0001_to_n1101 = g * T1 * n0001  # 23
    n0001_to_n0011 = g * T2 * n0001  # 24
    n0100_to_n0111 = g * T2 * n0100  # 25
    n0100_to_n1100 = g * T1 * n0100  # 26
    n0000_to_n1100 = g * T1 * n0000  # 27
    n0000_to_n0011 = g * T2 * n0000  # 28
    n1101_to_n1111 = g * T2 * n1101  # 29
    n0111_to_n1111 = g * T1 * n0111  # 30
    n1100_to_n1111 = g * T2 * n1100  # 31
    n0011_to_n1111 = g * T1 * n0011  # 32

    N_flow = np.array([n0101_pro, n0001_pro, n0100_pro, n0000_pro, n0111_pro, n0011_pro, n1101_pro, n1100_pro,
                       n1111_pro, n0101_to_n0001, n0101_to_n0100, n0001_to_n0101, n0100_to_n0101, n0001_to_n0000,
                       n0100_to_n0000, n0000_to_n0001, n0000_to_n0100, n1101_to_n1100, n1100_to_n1101, n0111_to_n0011,
                       n0011_to_n0111, n0101_to_n1101, n0101_to_n0111, n0001_to_n1101, n0001_to_n0011, n0100_to_n0111,
                       n0100_to_n1100, n0000_to_n1100, n0000_to_n0011, n1101_to_n1111, n0111_to_n1111, n1100_to_n1111,
                       n0011_to_n1111])
    N_flow = N_flow[:, :duration_5year]

    val_sum = np.tile(np.sum(N_flow, axis=1).reshape(N_flow.shape[0], 1), (1, N_flow.shape[1]))
    N_flow_normrow = np.divide(N_flow, val_sum)
    assert all([isclose(x, 1) for x in np.sum(N_flow_normrow, axis=1)])

    val_sum = np.tile(np.sum(N_flow, axis=0).reshape(1, N_flow.shape[1]), (N_flow.shape[0], 1))
    N_flow_normcol = np.divide(N_flow, val_sum)
    assert all([isclose(x, 1) for x in np.sum(N_flow_normcol, axis=0)])

    N_flow_normcol0 = np.zeros((len(cellflow), duration_5year))
    N_flow_normcol0extend = np.zeros((len(cellflow), duration_5year))
    N_flow_normcol0[:, :N_flow_normcol.shape[1]] = N_flow_normcol
    N_flow_normcol0extend[:, :N_flow_normcol.shape[1]] = N_flow_normcol
    N_flow_normcol0extend[:, N_flow_normcol.shape[1]:] = np.repeat(N_flow_normcol[:, -1].reshape(-1, 1),
                                                                   duration_5year-N_flow_normcol.shape[1],
                                                                   axis=1)
    assert all([isclose(x, 1) for x in np.sum(N_flow_normcol0[:, :N_flow_normcol.shape[1]], axis=0)])
    assert all([isclose(x, 1) for x in np.sum(N_flow_normcol0extend[:, :N_flow_normcol.shape[1]], axis=0)])
    N_flow_normcol = N_flow_normcol0
    N_flow_normcolextend = N_flow_normcol0extend

    return N_flow, N_flow_normcol, N_flow_normcolextend
