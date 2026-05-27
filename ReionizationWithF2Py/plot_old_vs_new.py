import numpy as np
import matplotlib.pyplot as plt

# ---------------------------------------------------
# LOAD DATA
# ---------------------------------------------------

z, dfcolldt_ion_old, dfcolldt_neut_old = np.loadtxt(
    'output_old_method.txt',
    unpack=True
)

z, dfcolldt_ion_new, dfcolldt_neut_new = np.loadtxt(
    'output_new_method.txt',
    unpack=True
)

# ---------------------------------------------------
# PLOT
# ---------------------------------------------------

fig, axes = plt.subplots(
    1, 2,
    figsize=(10, 4),
    sharex=True
)

# ---------------------------------------------------
# IONIZED
# ---------------------------------------------------

axes[0].plot(
    z,
    dfcolldt_ion_old/1e-19,
    label='Old'
)

axes[0].plot(
    z,
    dfcolldt_ion_new/1e-19,
    label='New'
)

axes[0].set_xlabel('z')
axes[0].set_ylabel(r'$df_{\rm coll}/dt$')
axes[0].set_title('Ionized')
axes[0].legend()

# ---------------------------------------------------
# NEUTRAL
# ---------------------------------------------------

axes[1].plot(
    z,
    dfcolldt_neut_old/1e-19,
    label='Old'
)

axes[1].plot(
    z,
    dfcolldt_neut_new/1e-19,
    label='New'
)

axes[1].set_xlabel('z')
axes[1].set_ylabel(r'$df_{\rm coll}/dt$')
axes[1].set_title('Neutral')
axes[1].legend()

# ---------------------------------------------------

plt.tight_layout()
plt.savefig('old_vs_new_method.pdf')
plt.show()

