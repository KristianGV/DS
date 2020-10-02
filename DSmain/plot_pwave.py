import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({
    "text.usetex": True
})

#dataset= np.genfromtxt('examples/aux/generic_wimp_oh2-planck-sigmav.dat')
# MAJORANA
i, m_p, sig1_p, sig1_p_errp, sig1_p_errm = np.loadtxt(
    'examples/aux/generic_wimp_oh2-planck-sigmav-pwave.dat', unpack=True, skiprows=7, usecols=range(0, 5))

i_s, m_s, sig1_s, sig1_s_errp, sig1_s_errm = np.loadtxt(
    'examples/aux/generic_wimp_oh2-planck-sigmav-swave.dat', unpack=True, skiprows=7, usecols=range(0, 5))

# DIRAC
i, m_p_d, sig1_p_d, sig1_p_errp_d, sig1_p_errm_d = np.loadtxt(
    'examples/aux/generic_wimp_oh2-planck-sigmav-pwave-dirac.dat', unpack=True, skiprows=7, usecols=range(0, 5))

i_s, m_s_d, sig1_s_d, sig1_s_errp_d, sig1_s_errm_d = np.loadtxt(
    'examples/aux/generic_wimp_oh2-planck-sigmav-swave-dirac.dat', unpack=True, skiprows=7, usecols=range(0, 5))


# MAJORANA
i_c, m_c_p, sig1_c_p, sig1_c_errp_p, sig1_c_errm_p = np.loadtxt(
    'examples/aux/thermal_sigmav_Majorana_SM_pwave.dat', unpack=True, skiprows=2, usecols=range(0, 5))

i_c, m_c_s, sig1_c_s, sig1_c_errp_s, sig1_c_errm_s = np.loadtxt(
    'examples/aux/thermal_sigmav_Majorana_SM_swave.dat', unpack=True, skiprows=2, usecols=range(0, 5))

# DIRAC
i_c, m_c_p_d, sig1_c_p_d, sig1_c_errp_p_d, sig1_c_errm_p_d = np.loadtxt(
    'examples/aux/thermal_sigmav_Dirac_SM_pwave.dat', unpack=True, skiprows=2, usecols=range(0, 5))

i_c, m_c_s_d, sig1_c_s_d, sig1_c_errp_s_d, sig1_c_errm_s_d = np.loadtxt(
    'examples/aux/thermal_sigmav_Dirac_SM_swave.dat', unpack=True, skiprows=2, usecols=range(0, 5))


# S WAVE MAJORANA

plt.figure(1)

plt.semilogx(m_s, sig1_s, label='s-wave', color='red', linewidth=.7)
plt.semilogx(m_s, sig1_s+sig1_s_errp, color='red', linewidth=.5)
plt.semilogx(m_s, sig1_s+sig1_s_errm, color='red', linewidth=.5)
plt.fill_between(m_s, sig1_s+sig1_s_errm, sig1_s +
                 sig1_s_errp, color='red', alpha=0.2)


plt.semilogx(m_c_s, sig1_c_s, label='s-wave correct',
             color='blue', linewidth=.7)
plt.semilogx(m_c_s, sig1_c_errp_s, color='blue', linewidth=.5)
plt.semilogx(m_c_s, sig1_c_errm_s, color='blue', linewidth=.5)
plt.fill_between(m_c_s, sig1_c_errm_s, sig1_c_errp_s, color='blue', alpha=0.2)

#plt.ylim(1E-25, 3.5E-25)
plt.title("Majorana s-wave")
plt.xlabel('WIMP mass')
plt.ylabel('($\sigma v$) [cm$^3$ s$^{-1}$]')
plt.legend()
plt.savefig('fig/majorana-swave.pdf', format='pdf')

plt.show()

# PERCENT PLOT

plt.figure(5)

lgmass = np.linspace(np.log(0.1e-2), np.log(0.3e6), 140)
mass = np.exp(lgmass)
sig1s = np.interp(mass, m_s, sig1_s)
sig1sc = np.interp(mass, m_c_s, sig1_c_s)

plt.semilogx(mass, np.absolute(sig1s-sig1sc)/sig1sc, color='red', linewidth=.7)


plt.title("Majorana s-wave difference")
plt.xlabel('WIMP mass')
plt.ylabel(
    r'$\frac{|\langle\sigma v\rangle_{correct}-\langle\sigma v\rangle|}{\langle\sigma v\rangle_{correct}}$')
plt.savefig('fig/majorana-swave-diff.pdf', format='pdf')
plt.show()


# P WAVE MAJORANA

plt.figure(2)
plt.semilogx(m_p, sig1_p, label='p-wave', color='red', linewidth=.7)
plt.semilogx(m_p, sig1_p+sig1_p_errp, color='red', linewidth=.5)
plt.semilogx(m_p, sig1_p+sig1_p_errm, color='red', linewidth=.5)
plt.fill_between(m_p, sig1_p+sig1_p_errm, sig1_p +
                 sig1_p_errp, color='red', alpha=0.2)

plt.semilogx(m_c_p, sig1_c_p, label='p-wave correct',
             color='blue', linewidth=.7)
plt.semilogx(m_c_p, sig1_c_errp_p, color='blue', linewidth=.5)
plt.semilogx(m_c_p, sig1_c_errm_p, color='blue', linewidth=.5)
plt.fill_between(m_c_p, sig1_c_errm_p, sig1_c_errp_p, color='blue', alpha=0.2)


# plt.ylim(1E-25, 3.5E-25)
plt.title("Majorana p-wave")
plt.xlabel('WIMP mass')
plt.ylabel('$b$')
plt.legend()
plt.savefig('fig/majorana-pwave.pdf', format='pdf')
plt.show()

# PERCENT PLOT

plt.figure(6)

lgmass = np.linspace(np.log(0.1e-2), np.log(0.3e6), 140)
mass = np.exp(lgmass)
sig1p = np.interp(mass, m_p, sig1_p)
sig1pc = np.interp(mass, m_c_p, sig1_c_p)

plt.semilogx(mass, np.abs(sig1p-sig1pc)/sig1pc, color='red', linewidth=.7)


plt.title("Majorana p-wave difference")
plt.xlabel('WIMP mass')
plt.ylabel(
    r'$\frac{|\langle\sigma v\rangle_{correct}-\langle\sigma v\rangle|}{\langle\sigma v\rangle_{correct}}$')
plt.savefig('fig/majorana-pwave-diff.pdf', format='pdf')

plt.show()


# S WAVE DIRAC

plt.figure(3)

plt.semilogx(m_s_d, sig1_s_d, label='s-wave', color='red', linewidth=.7)
plt.semilogx(m_s_d, sig1_s_d+sig1_s_errp_d, color='red', linewidth=.5)
plt.semilogx(m_s_d, sig1_s_d+sig1_s_errm_d, color='red', linewidth=.5)
plt.fill_between(m_s_d, sig1_s_d+sig1_s_errm_d, sig1_s_d +
                 sig1_s_errp_d, color='red', alpha=0.2)


plt.semilogx(m_c_s_d, sig1_c_s_d, label='s-wave correct',
             color='blue', linewidth=.7)
plt.semilogx(m_c_s_d, sig1_c_errp_s_d, color='blue', linewidth=.5)
plt.semilogx(m_c_s_d, sig1_c_errm_s_d, color='blue', linewidth=.5)
plt.fill_between(m_c_s_d, sig1_c_errm_s_d,
                 sig1_c_errp_s_d, color='blue', alpha=0.2)

#plt.ylim(1E-25, 3.5E-25)
plt.title("Dirac s-wave")
plt.xlabel('WIMP mass')
plt.ylabel('($\sigma v$) [cm$^3$ s$^{-1}$]')
plt.legend()
plt.savefig('fig/dirac-swave.pdf', format='pdf')

plt.show()

# PERCENT PLOT

plt.figure(7)

lgmass = np.linspace(np.log(0.1e-2), np.log(0.3e6), 140)
mass = np.exp(lgmass)
sig1sd = np.interp(mass, m_s_d, sig1_s_d)
sig1scd = np.interp(mass, m_c_s_d, sig1_c_s_d)

plt.semilogx(mass, np.abs(sig1sd-sig1scd)/sig1scd, color='red', linewidth=.7)


plt.title("Dirac s-wave difference")
plt.xlabel('WIMP mass')
plt.ylabel(
    r'$\frac{|\langle\sigma v\rangle_{correct}-\langle\sigma v\rangle|}{\langle\sigma v\rangle_{correct}}$')
plt.savefig('fig/dirac-swave-diff.pdf', format='pdf')
plt.show()


# P WAVE DIRAC

plt.figure(4)
plt.semilogx(m_p_d, sig1_p_d, label='p-wave', color='red', linewidth=.7)
plt.semilogx(m_p_d, sig1_p_d+sig1_p_errp_d, color='red', linewidth=.5)
plt.semilogx(m_p_d, sig1_p_d+sig1_p_errm_d, color='red', linewidth=.5)
plt.fill_between(m_p_d, sig1_p_d+sig1_p_errm_d, sig1_p_d +
                 sig1_p_errp_d, color='red', alpha=0.2)

plt.semilogx(m_c_p_d, sig1_c_p_d, label='p-wave correct',
             color='blue', linewidth=.7)
plt.semilogx(m_c_p_d, sig1_c_errp_p_d, color='blue', linewidth=.5)
plt.semilogx(m_c_p_d, sig1_c_errm_p_d, color='blue', linewidth=.5)
plt.fill_between(m_c_p_d, sig1_c_errm_p_d,
                 sig1_c_errp_p_d, color='blue', alpha=0.2)


# plt.ylim(1E-25, 3.5E-25)
plt.title("Dirac p-wave")
plt.xlabel('WIMP mass')
plt.ylabel('$b$')
plt.legend()
plt.savefig('fig/dirac-pwave.pdf', format='pdf')
plt.show()

# PERCENT PLOT

plt.figure(8)

lgmass = np.linspace(np.log(0.1e-2), np.log(0.3e6), 140)
mass = np.exp(lgmass)
sig1pd = np.interp(mass, m_p_d, sig1_p_d)
sig1pcd = np.interp(mass, m_c_p_d, sig1_c_p_d)

plt.semilogx(mass, np.abs(sig1pd-sig1pcd)/sig1pcd, color='red', linewidth=.7)


plt.title("Dirac p-wave difference")
plt.xlabel('WIMP mass')
plt.ylabel(
    r'$\frac{|\langle\sigma v\rangle_{correct}-\langle\sigma v\rangle|}{\langle\sigma v\rangle_{correct}}$')
plt.savefig('fig/dirac-pwave-diff.pdf', format='pdf')
plt.show()
