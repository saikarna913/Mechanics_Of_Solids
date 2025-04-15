from Euler_Bernoulli import EulerBernoulliBeam

# Common beam properties
length = 5.0  # m
E = 200e9  # Pa (steel)
I = 8.33e-6  # m^4
beam_height = 0.2  # m
udl_start = 0.0  # m
udl_end = 5.0  # m
udl_intensity = -2000  # N/m (downward)
x_analysis = 2.5
y_analysis = 0.05
z_analysis = 0.0

# Support types to analyze
support_types = ['simply_supported', 'cantilever', 'fixed_fixed']

for support in support_types:
    print(f"\n=== {support.replace('_', ' ').title()} Beam with UDL ===")
    
    beam = EulerBernoulliBeam(length, E, I, support_type=support)
    beam.add_distributed_load(udl_start, udl_end, udl_intensity)

    # Analyze
    beam.calculate_reactions()
    print("Reactions:", beam.reactions)

    # Plots
    beam.plot_sfd_bmd()
    beam.plot_deflection()

    # Stress and Mohr's Circle
    sigma_1, sigma_2, tau_max = beam.plot_mohr_circle(x_analysis, y_analysis, z_analysis, beam_height)
    print(f"Principal stresses at x={x_analysis}m: σ1={sigma_1:.2e} Pa, σ2={sigma_2:.2e} Pa")
    print(f"Maximum shear stress: {tau_max:.2e} Pa")
