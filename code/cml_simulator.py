"""
This code simulates CML progression based on a hierarchical model
for hematopoiesis.
"""


def construct_hierarchy(num_lev_var,
                        sc_number,
                        level_cell_increase,
                        total_daily_output,
                        mitotic_pool_loc,
                        gamma_pool_1):
    """
    Creating the arrays saving the different rates and the number of
    cell in homeostasis
    The rates will include differentiation and division
    """
    wild_type = np.zeros(num_lev_var)
    wild_diff = np.zeros(num_lev_var)
    delta_array = np.zeros(num_lev_var)
    wild_division = np.zeros(num_lev_var)

    for lev in range(num_lev_var):
        num_cells_mit_1 = sc_number * level_cell_increase ** (mitotic_pool_loc - 1)
        if lev < mitotic_pool_loc:
            wild_type[lev] = sc_number * level_cell_increase ** lev
        else:
            wild_type[lev] = num_cells_mit_1 * 2.0 ** (lev - (mitotic_pool_loc - 1))
    delta_array[-1] = total_daily_output

    # Compute the value of delta and compute the rates based on gamma
    for lev in range(num_lev_var - 2, -1, -1):
        if lev < mitotic_pool_loc - 1:
            delta_array[lev] = delta_array[lev + 1] / gamma_pool_1
        else:
            delta_array[lev] = delta_array[lev + 1] / 2.0

    wild_division[0] = 0.5 * delta_array[0] / wild_type[0]
    wild_diff[0] = wild_division[0]

    for i in range(1, num_lev_var):
        wild_division[i] = (delta_array[i]) * (0.5 - delta_array[i - 1] / delta_array[i]) / wild_type[i]
    for i in range(1, num_lev_var):
        wild_diff[i] = 0.5 * delta_array[i] / wild_type[i]

    return wild_type, wild_diff, wild_division


def mobil_dist(mito_1_loc, expo_param, sum_r):
    """
    Create the mobilization vector
    """
    mob_d = np.zeros_like(sum_r)
    for ind, tot_r in enumerate(sum_r):
        if ind < mito_1_loc:
            mob_d[ind] = tot_r
        else:
            mob_d[ind] = tot_r * expo_param
    return np.array([mob_d]).T


def get_gamma(num_lev_var,
              sc_number,
              level_cell_increase,
              total_daily_output,
              mitotic_pool_loc,
              sc_div_y
              ):
    """
    Find the appropriate gamma value 1 to run the simulation,
    Then we have the proper values for the homeostatic values
    """
    gamma_initial = 2.0
    high_div_flag = True
    while high_div_flag:
        _, wild_diff, _ = construct_hierarchy(num_lev_var,
                                              sc_number,
                                              level_cell_increase,
                                              total_daily_output,
                                              mitotic_pool_loc,
                                              gamma_initial)
        if 2 * wild_diff[0] * 365 > sc_div_y:
            gamma_initial += 0.001
        else:
            high_div_flag = False
    return construct_hierarchy(num_lev_var,
                               sc_number,
                               level_cell_increase,
                               total_daily_output,
                               mitotic_pool_loc,
                               gamma_initial)


def pois_fun(ele):
    """
    Simple application of a poisson distribution
    for low values
    """
    entries = np.where((ele > 0) & (ele < 10))
    for item in zip(entries[0], entries[1]):
        ele[item] = np.random.poisson(ele[item])
    return ele.round()


def mob_dyn(cell_r, mobi_var, mobi_base):
    """
    Function to get the dynamics of the mobilization
    """
    value = mobi_base * (np.exp((cell_r - 1) * mobi_var))
    if cell_r < 1.5:
        return 0
    return value


def bld_death_fun(bld_turn, diff_p, mito_1_loc):
    """
    Function to construct the death
    rate values for the bloodstream
    """
    bld_d = np.zeros_like(diff_p)
    for ind, d_p in enumerate(diff_p):
        if ind < mito_1_loc:
            bld_d[ind] = d_p / bld_turn
        else:
            bld_d[ind] = d_p
    bld_d[-1] = diff_p[-1] * 20
    return bld_d


def ind_simulation(mob_param,
                   blood_turn,
                   mut_rate,
                   sc_div_y,
                   fit_bcr,
                   mobi_intensity,
                   mobi_baseline,
                   point_id,
                   inst_id):
    """
    Function to perform the individual simulation and
    output the summary statistics
    """
    output_names = ['Latent time',
                    'Time to end the simulation',
                    "Crit group time",
                    'Bloodstream blast percentage at critical mutation',
                    'Bone marrow blast percentage at critical mutation',
                    'Cells in blood at critical mutation',
                    'Bone marrow blast percentage at the end',
                    "Bloodstream blast percentage at the end",
                    'Cells in blood at the end',
                    'Simulation code',
                    'Advanced phase duration',
                    'Chronic phase duration',
                    'star advanced phase duration',
                    'star chronic phase duration',
                    'star latent time',
                    'star blast phase duration',
                    'Mut prof lat',
                    'Mut prof chro',
                    'Mut prof adv',
                    'Mut prof crit',
                    "group per adv phase"
                    ]
    num_sc = 1e4
    mit_pool_1_loc = 15
    cell_growth_level = 2.769
    daily_output = 1e11
    num_levs = 20
    w_type, w_diff, w_division = get_gamma(num_levs,
                                           num_sc,
                                           cell_growth_level,
                                           daily_output,
                                           mit_pool_1_loc,
                                           sc_div_y
                                           )
    sum_rates = w_diff + w_division
    mob_dis = mobil_dist(mit_pool_1_loc, mob_param, sum_rates)
    max_fit = (w_diff[2] - w_division[2]) / (2 * sum_rates[2])
    fit_bcr_f = max_fit * fit_bcr
    num_crit_mut = 4
    fit = (max_fit - fit_bcr_f) / num_crit_mut
    fit += fit * 0.1
    num_mut = 15
    out_det_freq = 1e3
    inst_id_out_freq = 1
    bld_file_out = results_detail + 'blood_' + str(point_id) + '_' + str(inst_id) + '.txt'
    bm_file_out = results_detail + 'bone_' + str(point_id) + '_' + str(inst_id) + '.txt'
    out_flag = 999
    while out_flag == 999:
        bld_d_v = bld_death_fun(blood_turn, w_diff, mit_pool_1_loc)
        # Creating the arrays to store the cells and rates
        bone_cells = np.zeros((num_levs, num_mut))
        bld_cells = np.zeros((num_levs, num_mut))
        bld_death = np.zeros((num_levs, num_mut))
        bone_diff = np.zeros((num_levs, num_mut))
        bone_div = np.zeros((num_levs, num_mut))
        bone_cells[:, 0] = w_type
        bone_diff[:, 0] = w_diff
        bone_div[:, 0] = w_division
        # Creating the number of mature cells in the bloodstream
        bld_cells[-1, 0] = 5e10
        bld_death[:, 0] = bld_d_v
        simulation_time = 0
        num_steps = 0
        bm_coll_out = []
        bld_coll_out = []
        for clone in range(1, num_mut):
            if clone == 1:
                bone_div[:, 1] = bone_div[:, 0] + sum_rates * fit_bcr_f
                bone_diff[:, 1] = bone_diff[:, 0] - sum_rates * fit_bcr_f
                bld_death[:, 1] = bld_d_v
            else:
                bone_div[:, clone] = bone_div[:, 1] + sum_rates * fit * (clone - 1)
                bone_diff[:, clone] = bone_diff[:, 1] - sum_rates * fit * (clone - 1)
                bld_death[:, clone] = bld_d_v
        homeo_cells = np.sum(bone_cells)
        bone_cells[0, 1] = 1
        flag_current = 0
        phase_flag = 0
        star_flag = 0
        state_dic = {}
        for out_par in output_names:
            state_dic[out_par] = None
        while simulation_time < 365 * 20 and flag_current == 0:
            num_cells_bm = np.sum(bone_cells)
            num_cells_bld = np.sum(bld_cells)
            blast_cells_bm = np.sum(bone_cells[:15, :])
            blast_cells_bld = np.sum(bld_cells[:15, :])
            bld_per = (blast_cells_bld / num_cells_bld) * 100
            bone_per = (blast_cells_bm / num_cells_bm) * 100
            cell_ratio = num_cells_bm / homeo_cells
            max_cell_rate = np.max([np.max(bld_death),
                                    np.max(bone_div),
                                    np.max(bone_diff),
                                    np.max(mob_dis * mob_dyn(cell_ratio, mobi_intensity, mobi_baseline))
                                    ]
                                   )
            delta_val = max_cell_rate / 10
            if inst_id % inst_id_out_freq == 0 and num_steps % out_det_freq == 0:
                bm_out = np.concatenate([[np.full(num_mut, simulation_time)],
                                         [np.linspace(0, num_mut - 1, num=num_mut)],
                                         [np.full(num_mut, mob_dyn(cell_ratio, mobi_intensity, mobi_baseline))],
                                         [np.full(num_mut, bone_per)],
                                         bone_cells]).T.tolist()
                bm_coll_out.extend(bm_out)
                bld_out = np.concatenate([[np.full(num_mut, simulation_time)],
                                          [np.linspace(0, num_mut - 1, num=num_mut)],
                                          [np.full(num_mut, mob_dyn(cell_ratio, mobi_intensity, mobi_baseline))],
                                          [np.full(num_mut, bld_per)],
                                          bld_cells]).T.tolist()
                bld_coll_out.extend(bld_out)

            # Compute the number of events by division
            num_b_div = pois_fun(bone_cells * bone_div * delta_val)
            # Compute mutants by division from num of events except from first and last clones
            b_div_mut = pois_fun(num_b_div[:, 1:-1] * mut_rate)
            num_b_div[:, 1:-1] -= b_div_mut

            # Compute the number of events by differentiation
            num_b_diff = pois_fun(bone_cells * bone_diff * delta_val)
            # Compute the number of mutants by differentiation, except in the last and first mutants clone.
            b_dif_mut = pois_fun(num_b_diff[:, 1:-1] * mut_rate)
            num_b_diff[:, 1:-1] -= b_dif_mut

            # Compute the number of mobilization and death events
            num_b_mov = pois_fun(
                (bone_cells * mob_dis) * delta_val * mob_dyn(cell_ratio, mobi_intensity, mobi_baseline))
            num_bld_death = pois_fun(bld_death * bld_cells * delta_val)

            # Modifying the number of cells in the arrays
            bone_cells += num_b_div
            bone_cells[1:, :] += 2 * num_b_diff[:-1, :]
            bld_cells[-1, :] += 2 * num_b_diff[-1, :]
            bone_cells -= num_b_diff

            # Modifying the number of cells by mutants
            # add mutants in the next clone by division
            bone_cells[:, 2:] += b_div_mut
            # Add mutants in the next clone by differentiation
            bone_cells[1:, 2:] += 2 * b_dif_mut[:-1, :]

            # Move clones in the last level to the blood compartment
            bld_cells[-1, 2:] += 2 * b_dif_mut[-1, :]
            bone_cells[:, 2:] -= b_dif_mut

            # Apply the leaking
            bone_cells -= num_b_mov
            bld_cells += num_b_mov

            # Kill cells in the blood
            bld_cells -= num_bld_death

            simulation_time += delta_val
            num_steps += 1
            bone_cells[bone_cells < 0] = 0
            bld_cells[bld_cells < 0] = 0

            if bld_per > 2 and phase_flag == 0:
                state_dic['Latent time'] = simulation_time / 365.0
                mut_array = bone_cells.nonzero()[0]
                mut_array = np.unique(mut_array, return_counts=True)[1]
                state_dic['Mut prof lat'] = mut_array - 2
                phase_flag = 1
            if bld_per > 10 and phase_flag == 1:
                state_dic['Chronic phase duration'] = simulation_time / 365.0 - state_dic["Latent time"]
                mut_array = bone_cells.nonzero()[0]
                mut_array = np.unique(mut_array, return_counts=True)[1]
                tot_cells = np.sum(bone_cells)
                group_sum = np.sum(bone_cells, axis=0)
                group_sum /= tot_cells / 100
                state_dic["group per adv phase"] = group_sum[num_crit_mut + 2]
                state_dic['Mut prof chro'] = mut_array - 2
                phase_flag = 2
            if bld_per > 20 and phase_flag == 2:
                state_dic['Advanced phase duration'] = simulation_time / 365.0 - state_dic["Chronic phase duration"] - \
                                                       state_dic["Latent time"]
                mut_array = bone_cells.nonzero()[0]
                mut_array = np.unique(mut_array, return_counts=True)[1]
                state_dic['Mut prof adv'] = mut_array - 2
                phase_flag = 3
            if bone_cells[mit_pool_1_loc, num_crit_mut + 2] > 1e4 and star_flag == 0:
                state_dic["Crit group time"] = simulation_time / 365.0
                mut_array = bone_cells.nonzero()[0]
                mut_array = np.unique(mut_array, return_counts=True)[1]
                state_dic['Mut prof crit'] = mut_array - 2
                star_flag = 1
                state_dic["Bloodstream blast percentage at critical mutation"] = bld_per
                state_dic['Bone marrow blast percentage at critical mutation'] = bone_per
                state_dic['Cells in blood at critical mutation'] = num_cells_bld
                if bld_per <= 2:
                    state_dic["star latent time"] = simulation_time / 365.0
                else:
                    if bld_per <= 10:
                        state_dic['star chronic phase duration'] = simulation_time / 365.0 - state_dic['Latent time']
                        state_dic["star latent time"] = state_dic["Latent time"]
                    else:
                        if bld_per <= 20:
                            state_dic['star advanced phase duration'] = simulation_time / 365.0 - state_dic[
                                'Chronic phase duration'] - state_dic["Latent time"]
                            state_dic['star chronic phase duration'] = state_dic['Chronic phase duration']
                            state_dic["star latent time"] = state_dic["Latent time"]
                        else:
                            state_dic['star blast phase duration'] = simulation_time / 365.0 - state_dic[
                                'Advanced phase duration'] - state_dic['Chronic phase duration'] - state_dic[
                                                                         "Latent time"]
                            state_dic['star chronic phase duration'] = state_dic['Chronic phase duration']
                            state_dic['star advanced phase duration'] = state_dic['Advanced phase duration']
                            state_dic["star latent time"] = state_dic["Latent time"]

            if bone_cells[0, 1] == 0:
                # print('BCR-ABL stem cells extincted')
                flag_current = 999
            if bld_per > 30:
                # print('Blood percentage too high')
                flag_current = 1
        #            if bone_per > 30:
        # print('Bone marrow percentage too high')
        #                flag_current = 3
        #            if bone_cells[mit_pool_1_loc,6] > 1e4 and bld_per>20:
        # print('Neoplastic cells found')
        #                flag_current = 1
        state_dic['Time to end the simulation'] = simulation_time / 365.0
        state_dic['Bone marrow blast percentage at the end'] = bone_per
        state_dic["Bloodstream blast percentage at the end"] = bld_per
        state_dic['Cells in blood at the end'] = num_cells_bld
        state_dic['Simulation code'] = flag_current
        if flag_current != 999 and inst_id % inst_id_out_freq == 0:
            with open(bm_file_out, 'w+', newline='') as csvfile:
                bm_writer = csv.writer(csvfile)
                bm_writer.writerows(bm_coll_out)
            with open(bld_file_out, 'w+', newline='') as csvfile_d:
                bm_writer = csv.writer(csvfile_d)
                bm_writer.writerows(bld_coll_out)
        out_flag = flag_current
        # print(f"out_flag {out_flag}")
    return state_dic


def initiate_sim(vec_input):
    """
    function to initiate a set of 100 instances in each one
    of the state points
    """
    num_levs_var = 20
    list_names = ['Mobilization parameter',
                  'Blood turnover',
                  'Mutation rate',
                  'Stem cell division',
                  'Fitness BCR-ABL',
                  'Mobilization intensity',
                  "Mobilization baseline",
                  'Point space',
                  'Instance']

    output_names = ['Latent time',
                    'Time to end the simulation',
                    "Crit group time",
                    'Bloodstream blast percentage at critical mutation',
                    'Bone marrow blast percentage at critical mutation',
                    'Cells in blood at critical mutation',
                    'Bone marrow blast percentage at the end',
                    "Bloodstream blast percentage at the end",
                    'Cells in blood at the end',
                    'Simulation code',
                    'Advanced phase duration',
                    'Chronic phase duration',
                    'star advanced phase duration',
                    'star chronic phase duration',
                    'star latent time',
                    'star blast phase duration',
                    'Mut prof lat',
                    'Mut prof chro',
                    'Mut prof adv',
                    'Mut prof crit',
                    "group per adv phase"
                    ]
    crucial_output = ['Latent time',
                      'Time to end the simulation',
                      "Crit group time",
                      'Bloodstream blast percentage at critical mutation',
                      'Bone marrow blast percentage at critical mutation',
                      'Cells in blood at critical mutation',
                      'Bone marrow blast percentage at the end',
                      "Bloodstream blast percentage at the end",
                      'Cells in blood at the end',
                      'Simulation code',
                      'Advanced phase duration',
                      'Chronic phase duration',
                      'star advanced phase duration',
                      'star chronic phase duration',
                      'star latent time',
                      'star blast phase duration']

    expli_out_names = []
    for para in output_names:
        if "prof" in para:
            for ele in range(num_levs_var):
                expli_out_names.append(para + "_" + str(int(ele)))
        else:
            expli_out_names.append(para)

    full_list_names = list_names + expli_out_names
    # print(len(expli_out_names))
    copy_point = copy.deepcopy(vec_input)
    point_df = pd.DataFrame(columns=full_list_names)
    error_flag = 0
    num_inst = 0
    while num_inst < 1 and error_flag == 0:
        # print(f"Num inst {num_inst}")
        new_vec = copy_point + [num_inst]
        out_dic = ind_simulation(*new_vec)
        out_list = []
        for para in output_names:
            if "prof" in para:
                if out_dic[para] is not None:
                    for ele in range(num_levs_var):
                        # Minus one is due to the fact of not counting the BCR-ABL mutant
                        out_list.append(out_dic[para][ele])
                else:
                    out_list.extend(np.full(num_levs_var, None))
            else:
                out_list.append(out_dic[para])
                # print(para,out_dic[para])
        #        print(full_list_names)
        #        print(out_list)
        fin_list = new_vec + out_list
        # print(len(fin_list))
        # print(len(full_list_names))
        point_df.loc[num_inst] = fin_list
        # print(point_df.loc[num_inst][crucial_output])
        # print(point_df.loc[num_inst][["Chronic phase duration","Advanced phase duration"]])
        # print(f"chro median {point_df['Chronic phase duration'].median()}")
        # print(f"adv median {point_df['Advanced phase duration'].median()}")
        num_err = point_df[point_df['Simulation code'] != 1.0].shape[0]
        num_rows = point_df.shape[0]
        if num_rows > 50:
            med_chro = point_df["Chronic phase duration"].median()
            med_adv = point_df["Advanced phase duration"].median()
            if med_chro < 3.0 or med_chro > 5.0:
                print("bad median chro")
                error_flag = 1
                point_df.to_csv(results_summary + 'summ_dat_' + str(vec_input[-1]) + '.txt',
                                sep=',',
                                header=False,
                                index=False)
            if med_adv < 1.0 or med_adv > 2.0:
                print("bad adv median")
                error_flag = 1
                point_df.to_csv(results_summary + 'summ_dat_' + str(vec_input[-1]) + '.txt',
                                sep=',',
                                header=False,
                                index=False)
        if num_err > 10:
            print('Error point')
            error_flag = 1
            point_df.to_csv(results_summary + 'summ_dat_' + str(vec_input[-1]) + '.txt',
                            sep=',',
                            header=False,
                            index=False)
            break
        num_inst += 1
    if error_flag == 0:
        # print('Successful point')
        # point_df.to_csv('/home/mario/Desktop/summ_data_'+str(vec_input[-1])+'.txt',
        point_df.to_csv(results_summary + 'summ_dat_' + str(vec_input[-1]) + '.txt',
                        sep=',',
                        header=False,
                        index=False)
    return None


def output_args_vec(job_var):
    """
    Determine the combination of parameters to start
    """
    mut_rate_v = [0]
    blood_turn_v = [50]
    expo_mobi_v = [50]
    sc_div_year_v = np.linspace(0.5, 7.0, 20)
    fit_b_v = np.linspace(0.2, 0.8, 20)
    mobi_base_v = [0.03]
    mobi_inte_v = [0.1]

    grid_list = []
    counter = 0
    for mut_rate in mut_rate_v:
        for fit_bcr in fit_b_v:
            for blood_turn in blood_turn_v:
                for mob_param in expo_mobi_v:
                    for sc_d_y in sc_div_year_v:
                        for mobi_i in mobi_inte_v:
                            for mobi_base in mobi_base_v:
                                counter += 1
                                grid_list.append([mob_param,
                                                  blood_turn,
                                                  mut_rate,
                                                  sc_d_y,
                                                  fit_bcr,
                                                  mobi_i,
                                                  mobi_base,
                                                  counter])
    # print(f'len of array{len(grid_list)}')
    # print(grid_list[job_var])
    return grid_list[job_var]


if __name__ == "__main__":
    import copy
    import csv
    import time
    import numpy as np
    import pandas as pd
    # import argparse
    import os
    import csv

    print(f"Simulation starts")

    code_starts = time.time()
    # pr = cProfile.Profile()
    # pr.enable()

    # parser = argparse.ArgumentParser()
    # parser.add_argument('job_id', action='store', type=int)
    # args = parser.parse_args()

    # This array holds the model parameters of the individual simulation with the following order
    # beta, alpha, mu, stem cell division rate per year, BCR-ABL mutation strength, K, kappa, group_sim_id
    test_arr = [50, 50, 0, 6, 0.4, 0.1, 0.03, 3]
    # initiate_sim(test_arr)
    if not os.path.exists("../results"):
        os.mkdir("../results")
    results_summary = '../results/'
    results_detail = '../results/'
    # initiate_sim(load_parameter_fun(args.job_id, "/scratch/mape1416/cml/current_data/driver_main_points.csv"))
    initiate_sim(test_arr)
    # pr.disable()
    # s = io.StringIO()
    # sortby = "cumulative"
    # ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
    # ps.print_stats()
    code_ends = time.time()
    print(f" Time for the simulation {code_ends - code_starts} s")
