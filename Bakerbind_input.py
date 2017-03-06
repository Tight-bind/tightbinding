# General parameters
job_name = "Bulk Silicon"
sim_box_x = 5.431
sim_box_y = 5.431
sim_box_z = 5.431
dist_cut_off = 2.7
smear_temp = 0.001

# Atomic parameters
coordinate_file = "bulkSi.coord"
isfractionalcoord = True

eta_dict = {"s_s_sigma": -1.938,
            "s_p_sigma": 1.745,
            "p_s_sigma": 1.745,
            "p_p_sigma": 3.050,
            "p_p_pi": -1.075}

on_site_dict = {"Cp": -8.97,
                "Cs": -17.52,
                "Hs": -13.61,
                "Sip": -5.75,
                "Sis": -12.2}

# Solution parameters
sol_method = "SP2"
sparsity = "sparse"

# MP mesh dimensions
mpx = 2
mpy = 2
mpz = 2
