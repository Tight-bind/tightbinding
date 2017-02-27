
cdef slater_koster_table(orb_i, orb_j,
                        const double x_dir_cos,const double y_dir_cos,
                        const double z_dir_cos, const double dist_ij)


cdef on_site_energy_table(species, orbital)

cdef bloch_phase_factor(const double [:] k_point, const double [:] r_ij)
