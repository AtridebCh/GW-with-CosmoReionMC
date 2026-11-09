import numpy as np
import scipy as sp
from numpy import genfromtxt
from scipy.integrate import solve_ivp, simpson
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt


# ======================================================================
# Physical constants — CGS throughout
# ======================================================================
C_LIGHT      = 2.99792458e+10   # cm/s
H_PLANCK     = 6.6260755e-27    # erg·s
K_BOLTZ      = 1.380658e-16     # erg/K
M_ELECT      = 9.1093897e-28    # g
M_H_ATOM     = 1.673725e-24     # g
G_NEWT       = 6.67259e-08      # cm³/(g·s²)
A_RAD        = 7.565914e-15     # erg/(cm³·K⁴)
SIGMA_T      = 6.65e-25         # cm²
MPC_TO_CM    = 3.0857e+24       # cm/Mpc
SOLAR_MASS_G = 1.989e+33        # g
YR_TO_SEC    = 3.1557e+07       # yr/s
DELTA_NU = 8.1*1e14

# Astrophysical
Y_HE         = 0.25
LAMBDA_ALPHA = 1216e-8          # cm, Lyman-alpha
LAMBDA_LL    = 912e-8           # cm, Lyman-limit


class Generate21cmSignal:
    """
    Computes the 21-cm brightness temperature signal.
    All internal calculations in CGS.

    Parameters
    ----------
    H0           : float — Hubble parameter (km/s/Mpc)
    ombh2        : float — baryon density parameter
    omch2        : float — cold dark matter density parameter
    redshift     : array — redshift grid from reionization model
    logEps_X     : array — log10 X-ray emissivity (erg/s/cMpc³) on redshift grid
    lognDotAlpha : array — log10 Ly-alpha photon emissivity (/s/cMpc³) on redshift grid
    Q_HII        : array — ionized fraction on redshift grid
    Z_21cm       : array — redshift grid for the 21-cm signal output
    f_X          : float — X-ray heating efficiency (default 0.2, Furlanetto 06)
    f_alpha      : float — Ly-alpha coupling efficiency (default 1.0)
    """

    def __init__(self, H0, ombh2, omch2, redshift,
                 logEps_X, lognDotAlpha, Q_HII, Z_21cm,
                 f_X, f_alpha):

        self.h       = H0 / 100.0
        self.omega_b = ombh2 / self.h ** 2
        self.omega_m = (ombh2 + omch2) / self.h ** 2
        self.f_X     = 10**f_X
        self.f_alpha = 10**f_alpha
        self.Z_21cm  = Z_21cm

        # Interpolators from reionization model output
        self.log10epsilonX_interp  = sp.interpolate.interp1d(redshift, logEps_X, fill_value = 'extrapolate')
        self.Q_HII_interp          = sp.interpolate.interp1d(redshift, Q_HII, fill_value = 'extrapolate')
        self.log10nDotAlpha_interp = sp.interpolate.interp1d(redshift, lognDotAlpha, fill_value = 'extrapolate')

        # Comoving hydrogen number density (cm^-3)
        # n_H = rho_b * (1 - Y_He) / m_H  at z=0
        self.n_H = 2.0e-7 * (ombh2 / 0.022) * (1.0 - Y_HE)

        # Ionized fraction on the 21-cm redshift grid
        self.Q_HII_21cm = self.Q_HII_interp(self.Z_21cm)
        

    # ------------------------------------------------------------------
    # Hubble parameter (s^-1)
    # ------------------------------------------------------------------

    def Hubble(self, z):
        # H0: km/s/Mpc → cm/s/cm
        H0_cgs = self.h * 100.0 * 1.0e5 / MPC_TO_CM
        return H0_cgs * np.sqrt(self.omega_m * (1.0 + z) ** 3
                                + (1.0 - self.omega_m))

    # ------------------------------------------------------------------
    # CMB temperature (K)
    # ------------------------------------------------------------------

    def CMB_temp(self, z):
        return 2.73 * (1.0 + z)
    
    # ------------------------------------------------------------------
    #collisional coefficient
    # ------------------------------------------------------------------
    def collision_coeff(self, z):
        X_c=(0.068/(2.85*10**-15*self.T_bg))*10**-7*(1+z)**3*(3.1*10**-11*self.T_k**0.357*np.exp(-31/self.T_bg)+self.Q_HII_interp(z)*10**(-9.607+0.5*np.log10(self.T_k)*np.exp(-(np.log10(self.T_k**4.5)/1800.0))))
        return X_c

    # ------------------------------------------------------------------
    # Ly-alpha background
    # ------------------------------------------------------------------

    def _dtbydz(self, z):
        """dt/dz (s)."""
        return 1.0 / ((1.0 + z) * self.Hubble(z))

    def _integrand_lyalpha(self, z):
        """Ly-alpha emissivity × dt/dz  (/cm³)."""
        return 10.0 ** self.log10nDotAlpha_interp(z) * self._dtbydz(z)

    def J_alpha(self, Z_array):
        """
        Ly-alpha specific intensity at each redshift (photons/cm²/s/Hz/sr).
        Integrates over photons emitted between z and z_max (Lyman-limit redshift).
        """
        J = np.zeros(len(Z_array))
        for i, z in enumerate(Z_array):
            z_max  = (LAMBDA_ALPHA / LAMBDA_LL) * (1.0 + z) - 1.0
            z_grid = np.linspace(z, z_max, int(abs(z_max - z) / 0.01) + 1)
            integrand = self._integrand_lyalpha(z_grid)
            J[i] = self.f_alpha * (simpson(integrand, z_grid)
                    * C_LIGHT * (1.0 + z) ** 3
                    / (DELTA_NU*4.0 * np.pi * MPC_TO_CM ** 3))
        return J

    def x_alpha(self, z):
        """Wouthuysen-Field Ly-alpha coupling coefficient (dimensionless)."""
        return 1.81e11 * self.J_alpha(z) * (1.0 + z) ** -1

    # ------------------------------------------------------------------
    # Gas kinetic temperature ODE
    # ------------------------------------------------------------------

    def _dTk_dz(self, z, T_k):
        """
        dT_k/dz including adiabatic cooling and X-ray heating. CGS.
        Furlanetto 2006.

        epsilon_X  : erg/s/cm³
        n_H        : cm^-3
        K_BOLTZ    : erg/K
        Hubble     : s^-1
        """
        adiabatic = 2.0 * T_k / (1.0 + z)

        epsilon_X_cgs = (10.0 ** self.log10epsilonX_interp(z))/(MPC_TO_CM)**3           # erg/s/cm³

        xray_heat = (2.0 / 3.0 * self.f_X
                     * epsilon_X_cgs*self._dtbydz(z)
                     / (K_BOLTZ *self.n_H ))
        
        dT_kdz   =  adiabatic - xray_heat 
        return dT_kdz

    # ------------------------------------------------------------------
    # Main signal generator
    # ------------------------------------------------------------------

    def signal_generator(self):
        """
        Integrate T_k and compute the 21-cm brightness temperature.

        Returns
        -------
        T_b : array — brightness temperature in mK on Z_21cm grid
        """
        z_init = self.Z_21cm[0]
        z_end  = self.Z_21cm[-1]

        # Initial condition: adiabatic cooling from z=50
        T_k0 = [(1.0 + z_init) ** 2 / (1.0 + 50.0) ** 2 * 50.0]

        sol = solve_ivp(
            self._dTk_dz,
            (z_init, z_end),
            T_k0,
            method='Radau',
            t_eval=self.Z_21cm,
        )
        self.T_k = sol.y[0]

        self.T_bg = self.CMB_temp(self.Z_21cm)
        x_a  = self.x_alpha(self.Z_21cm)
        x_c  = self.collision_coeff(self.Z_21cm)

        # Spin temperature
        T_s = (1.0 + x_a + x_c) / (self.T_k ** -1 * (x_a + x_c) + self.T_bg ** -1)
        

        # Brightness temperature (mK)
        T_b = (10.1 * (1.0 - self.Q_HII_21cm)
               * (1.0 - self.T_bg / T_s)
               * (1.0 + self.Z_21cm) ** 0.5)

        return T_b


'''
redshift_edges, bright_temp = genfromtxt(
    '/home/atrideb/Downloads/figure2_plotdata.csv',
    delimiter=',', skip_header=27, usecols=(1, 4), unpack=True
)

Z, x_HII, tau, dNLLdz, Gamma_PI_12, pop2, pop3, x_HI5 = running_code_with_21(
    71.22007, 0.02259486, 0.1259488, 0.819696266003826,
    0.9480701, 0.004602334, 0.0001, 5.451408
)

signal = Generating21cmSignal(
    H0=71.22007, ombh2=0.02259486, omch2=0.1259488,
    redshift=Z, logEps_X=pop2, lognDotAlpha=pop3,
    Q_HII=x_HII, Z_21cm=redshift_edges,
)
T_bright = signal.signal_generator()

plt.plot(redshift_edges, T_bright,         label='fiducial')
plt.plot(redshift_edges, bright_temp*1000, label='EDGES')
plt.xlabel('z')
plt.ylabel(r'$T_b$ (mK)')
plt.legend()
plt.show()
'''
