import emcee

backend = emcee.backends.HDFBackend("test_run_CMB_Reion_GW.h5", read_only=True)

print("iterations =", backend.iteration)

chain = backend.get_chain()
print("chain shape =", chain.shape)

logp = backend.get_log_prob()
print("log_prob shape =", logp.shape)

last = backend.get_last_sample()
print(type(last))
print(last.coords.shape)
print(last.log_prob.shape)
