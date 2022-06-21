import os
import pandas as pd
import plotly.graph_objects as go
from matplotlib import colors
import numpy as np
from plotly.subplots import make_subplots
import matplotlib
import matplotlib.pyplot as plt

crit_names = ['Mobilization parameter',
              'Blood turnover',
              'Mutation rate',
              'Stem cell division',
              'Fitness BCR-ABL',
              'Mobilization intensity',
              "Mobilization baseline",
              'Latent time',
              'Simulation code',
              'Advanced phase duration',
              'Chronic phase duration',
              'Instance']


def phase_filter(dataset):
    median = dataset.groupby("Point space").median()
    median = median[(median["Chronic phase duration"]>3.0) & (median["Chronic phase duration"] < 5.0)]
    median = median[(median["Advanced phase duration"]>1.0) & (median["Advanced phase duration"] < 2.0)]
    return median


def import_sum_dataset(name):
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
                    "crit output per"
                    ]

    expli_out_names = []
    for para in output_names:
        if "prof" in para:
            for ele in range(20):
                expli_out_names.append(para + "_" + str(int(ele)))
        else:
            expli_out_names.append(para)
    full_list_names = list_names + expli_out_names
    # expli_out_names

    raw_data = pd.read_csv(name,
                           header=None,
                           names=full_list_names
                           )
    return raw_data


def plot_medians(f_g_data):
    fig6 = go.Figure()
    x_lat = f_g_data[(f_g_data['Latent time'] > 0)]['Latent time']
    print(f"Percentage of simulations with a latent phase larger than {len(x_lat[x_lat > 10]) * 100 / len(x_lat)}")
    x_chro = f_g_data[(f_g_data['Chronic phase duration'] > 0)]['Chronic phase duration']
    x_adv = f_g_data['Advanced phase duration']
    print(
        f"Latent phase quartiles first {x_lat.quantile(0.25)} "
        f"median {x_lat.quantile(0.5)} third {x_lat.quantile(0.75)} "
        f"std_dev {np.std(x_lat)} mean {np.mean(x_lat)}")
    print("Chronic quartiles first {} median {} third {} std {}".format(x_chro.quantile(0.25),
                                                                        x_chro.quantile(0.5),
                                                                        x_chro.quantile(0.75),
                                                                        np.std(x_chro)))
    print("Advanced quartiles first {} median {} third {} std {}".format(
        f_g_data['Advanced phase duration'].quantile(0.25),
        f_g_data['Advanced phase duration'].quantile(0.5),
        f_g_data['Advanced phase duration'].quantile(0.75),
        np.std(f_g_data['Advanced phase duration'])
    ))
    fig6.add_trace(go.Box(y=x_lat,
                          name='Latent',
                          marker_color="#91cf60"
                          # boxpoints="all"
                          )
                   )
    fig6.add_trace(go.Box(y=x_chro,
                          name='Chronic',
                          marker_color="#fecc5c"
                          # boxpoints="all"
                          )
                   )
    fig6.add_trace(go.Box(y=x_adv,
                          name='Advanced',
                          marker_color="#fc8d59"
                          # boxpoints="all"
                          )
                   )
    fig6.add_annotation(
        x=0.42,
        y=9.9,
        text="median: 9.7",
        showarrow=False
    )
    fig6.add_annotation(
        x=0.36,
        y=8.9,
        text="q1: 8.6",
        showarrow=False
    )
    fig6.add_annotation(
        x=0.36,
        y=10.9,
        text="q3: 10.9",
        showarrow=False
    )
    fig6.add_annotation(
        x=1.36,
        y=3.9,
        text="q3: 3.5",
        showarrow=False
    )
    fig6.add_annotation(
        x=1.42,
        y=3.18,
        text="median: 3.2",
        showarrow=False
    )
    fig6.add_annotation(
        x=1.36,
        y=2.5,
        text="q1: 2.9",
        showarrow=False
    )
    fig6.add_annotation(
        x=2.36,
        y=2.33,
        text="q3: 1.9",
        showarrow=False
    )
    fig6.add_annotation(
        x=2.42,
        y=1.71,
        text="median: 1.7",
        showarrow=False
    )
    fig6.add_annotation(
        x=2.36,
        y=1.12,
        text="q1: 1.6",
        showarrow=False
    )
    # fig2.update_layout(barmode='overlay')
    # fig6.update_layout(title_text="Chronic phase duration histogram",
    # title_font_size=30,
    # legend_title_font_size=15,
    # overlay=None
    # )
    fig6.update_yaxes(title_text='Time (years)',
                      title_font_size=20,
                      showgrid=True,
                      tickvals=np.arange(0, 21, 2))
    fig6.update_layout(margin=dict(l=1, b=10, t=25, r=10),
                       width=800,
                       height=500,
                       showlegend=False)
    return fig6


def create_time_plot(group_var, inst_var):
    """
    In this function the first line is a path to the folder that contains the individual
    simulations data
    :param group_var:
    :param inst_var:
    :return:
    """
    ind_data_paths = '../data'
    list_files = os.listdir(ind_data_paths)
    bone_list = [ind_data_paths + file for file in list_files if "bone" in file]
    bld_list = [ind_data_paths + file for file in list_files if "blood" in file]
    specific_inst_bone = ind_data_paths + "/bone_" + str(group_var) + "_" + str(inst_var) + ".txt"
    specific_inst_blood = ind_data_paths + "/blood_" + str(group_var) + "_" + str(inst_var) + ".txt"
    bm_data = pd.read_csv(specific_inst_bone, header=None)
    blood_data = pd.read_csv(specific_inst_blood, header=None)
    bm_data.drop_duplicates(inplace=True)
    blood_data.drop_duplicates(inplace=True)

    time_vec = bm_data[0].unique()
    bm_blast_vec = []
    blood_blast_vec = []
    blood_mature = []
    blood_tot = []
    bm_tot = []
    num_sc_vec = []
    num_mut_14 = []
    # num_mut_15 = []
    # mobilization = []
    latent_flag = False
    latent_time = 0
    cp_flag = False
    cp_time = 0
    adv_flag = False
    adv_time = 0
    for step in time_vec:
        bm_filter_data = bm_data[bm_data[0] == step]
        blood_filter_data = blood_data[blood_data[0] == step]
        # mob_val = bm_data[bm_data[0] == step][3].iloc[0]
        num_mut_sc = len(bm_filter_data[4].to_numpy().nonzero()[0]) - 2
        mut_mit_1 = len(bm_filter_data[18].to_numpy().nonzero()[0]) - 2
        bm_mit2_cells = bm_filter_data.iloc[:, 15:].sum().sum()
        bld_mit2_cells = blood_filter_data.iloc[:, 15:].sum().sum()
        blood_true_mature = blood_filter_data.iloc[:, -1].sum().sum()
        bm_total_cells = bm_filter_data.sum().sum()
        blood_total_cells = blood_filter_data.sum().sum()
        bm_blast_perce = bm_filter_data.iloc[0, 3]
        blood_blast_perce = blood_filter_data.iloc[0, 3]
        # print(bm_blast_perce,blood_blast_perce)
        if not latent_flag and blood_blast_perce > 2:
            latent_time = step
            latent_flag = True
        if not cp_flag and blood_blast_perce > 10:
            cp_time = step
            cp_flag = True
        if not adv_flag and blood_blast_perce > 20:
            adv_time = step
            adv_flag = True
        bm_blast_vec.append(bm_blast_perce)
        blood_tot.append(blood_total_cells)
        blood_blast_vec.append(blood_blast_perce)
        num_sc_vec.append(num_mut_sc)
        num_mut_14.append(mut_mit_1)
        bm_tot.append(bm_total_cells)
        blood_mature.append(blood_true_mature)
        # mobilization.append(mob_val)

    # print(latent_time,cp_time,adv_time)
    num_mut_14 = [0 if ele < 0 else ele for ele in num_mut_14]

    time_vec = time_vec / 365

    # print(latent_time / 365, (cp_time - latent_time) / 365, (adv_time - cp_time) / 365)

    fig = make_subplots(
        rows=3, cols=1,
        shared_xaxes=True,
        vertical_spacing=0.07)
    m_size = 1
    line_s = 2
    mode = "lines"
    fig.add_trace(go.Scatter(x=time_vec,
                             y=num_sc_vec,
                             mode="markers",
                             marker=dict(color="black",
                                         symbol="line-ns",
                                         line_width=m_size
                                         ),
                             name='Stem cell level'),
                  row=1, col=1)
    fig.add_trace(go.Scatter(x=time_vec,
                             y=num_mut_14,
                             mode="markers",
                             marker=dict(color="black",
                                         symbol="line-ew",
                                         line_width=m_size
                                         ),
                             name='Last level of mitotic pool I',
                             ), row=1, col=1)
    fig.add_trace(go.Scatter(x=time_vec,
                             y=bm_blast_vec,
                             mode=mode,
                             line=dict(color="#0000B2",
                                       width=line_s,
                                       dash="solid"
                                       ),
                             name='Bone marrow',
                             ), row=2, col=1)
    fig.add_trace(go.Scatter(x=time_vec,
                             y=blood_blast_vec,
                             mode=mode,
                             line=dict(color="#B20000",
                                       width=line_s,
                                       dash="solid"
                                       ),
                             name='Bloodstream',
                             ), row=2, col=1)
    fig.add_trace(go.Scatter(x=time_vec,
                             y=bm_tot,
                             mode=mode,
                             line=dict(color="#0000B2",
                                       width=line_s,
                                       dash="solid"
                                       ),
                             showlegend=False,
                             ), row=3, col=1)
    fig.add_trace(go.Scatter(x=time_vec,
                             y=blood_tot,
                             mode=mode,
                             line=dict(color="#B20000",
                                       width=line_s,
                                       dash="solid"
                                       ),
                             # marker=dict(color="#B20000",
                             #             size=m_size),
                             showlegend=False,
                             ), row=3, col=1)
    fig.add_trace(go.Scatter(x=time_vec,
                             y=blood_mature,
                             mode=mode,
                             line=dict(color="#B20000",
                                       width=line_s,
                                       dash="dot"
                                       ),
                             name='Mature cells in the bloodstream',
                             ), row=3, col=1)
    fig.add_vrect(x0=0,
                  x1=latent_time / 365.0,
                  line_dash="dash",
                  annotation_text="Latent phase",
                  annotation_position="top left",
                  fillcolor="#4dac26",
                  opacity=0.3,
                  line_width=0)
    fig.add_vrect(x0=latent_time / 365.0,
                  x1=cp_time / 365.0,
                  line_dash="dash",
                  annotation_text="CP",
                  annotation_position="left",
                  fillcolor="#b8e186",
                  opacity=0.3,
                  line_width=0)
    fig.add_vrect(x0=cp_time / 365.0,
                  x1=adv_time / 365.0,
                  line_dash="dash",
                  annotation_text="AP",
                  annotation_position="left",
                  fillcolor="#f1b6da",
                  opacity=0.3,
                  line_width=0)
    fig.add_vrect(x0=adv_time / 365.0,
                  x1=step / 365.0,
                  line_dash="dash",
                  annotation_text="BP",
                  annotation_position="left",
                  fillcolor="#d01c8b",
                  opacity=0.3,
                  line_width=0)
    fig.update_layout(
        # title_text='Computer simulation of CML progression',
        margin=dict(l=0, r=0, t=25, b=0),
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=-0.25,
            xanchor="left",
            x=0,
            font_size=13,
            itemsizing='constant'
        )
    )
    fig.update_xaxes(title='Time (years)',
                     dtick=1,
                     row=3,
                     col=1,
                     range=[-0.01, 18.4]
                     )
    fig.update_yaxes(title='Additional mutations',
                     dtick=2,
                     title_font_size=15,
                     row=1,
                     col=1)
    fig.update_yaxes(title='Blast percentage',
                     title_font_size=15,
                     dtick=5,
                     row=2,
                     col=1)
    fig.update_yaxes(title='Cells',
                     title_font_size=15,
                     exponentformat="power",
                     showexponent="first",
                     row=3,
                     col=1,
                     type='log')
    fig.write_image("../results/time_graph.eps",
                    width=800,
                    height=500,
                    engine="kaleido"
                    )
    # fig.show(width=1000,
    #          height=500)


def get_histogram_level(data, phase):
    final_array = np.full((14, 20), np.nan, dtype=float)
    lat_col = [ele for ele in data.columns.tolist() if phase + '_' in ele]
    lat_data = data.loc[:, lat_col]
    len_norm = lat_data.shape[0]
    for level in range(20):
        lev_dat = lat_data.iloc[:, level]
        for mutation in range(14):
            val_count = np.count_nonzero(lev_dat == mutation)
            if val_count > 0.0:
                final_array[mutation][level] = (val_count/len_norm)*100
    lat_data.dropna(inplace=True)
    lat_data = lat_data.apply(lambda x: x.value_counts(), axis=0)
    lat_data = lat_data.apply(lambda x: x / x.sum(), axis=0)
    lat_data *= 100
    return final_array


def create_mutation_profile(data_driver_sim):
    data_driver_sim = data_driver_sim.round({"Fitness BCR-ABL": 3})
    first_row = data_driver_sim[data_driver_sim["Fitness BCR-ABL"] == 0.2]
    second_row = data_driver_sim[data_driver_sim["Fitness BCR-ABL"] == 0.421]
    third_row = data_driver_sim[data_driver_sim["Fitness BCR-ABL"] == 0.611]
    data_list = [first_row, second_row, third_row]
    plot_list = []
    for dataset in data_list:
        aux_list = []
        for phase in ["lat", "chro", "adv"]:
            aux_list.append(get_histogram_level(dataset, phase))
        plot_list.append(aux_list)
    x_ticks = np.arange(0, 14)
    y_ticks = np.arange(0, 20)
    tick_label_size = 10

    fig, axs = plt.subplots(3, 3, figsize=(12, 12), sharey=True)
    fig.tight_layout()
    for vec in axs:
        for axis in vec:
            axis.set_xticks(x_ticks)
            axis.set_xticklabels(axis.get_xticks(), rotation=-60)
            axis.set_yticks(y_ticks)
            for tick in axis.get_yticklabels()[-5:]:
                tick.set_color("red")
            axis.tick_params(axis='both', which='major', labelsize=tick_label_size)
            axis.grid(zorder=-5.0)
            for tick in axis.get_ygridlines()[-5:]:
                tick.set_color("red")
    already_flag = 0
    for axis, data in zip(axs, plot_list):
        for ele in range(3):
            if already_flag == 0:
                im_col = axis[ele].imshow(data[ele].T, origin="lower", zorder=2.0, cmap="Blues", vmin=0, vmax=100)
                already_flag = 1
            else:
                axis[ele].imshow(data[ele].T, origin="lower", zorder=2.0, cmap="Blues", vmin=0, vmax=100)

    for item in range(3):
        axs[2][item].set_xlabel("Number of driver mutations", fontsize=10)
        axs[item][0].set_ylabel("Level in the hierarchy")
    axs[0][0].set_title("Latent")
    axs[0][1].set_title("Chronic")
    axs[0][2].set_title("Advanced")
    fig.subplots_adjust(right=0.8, wspace=0.05)
    cbar_ax = fig.add_axes([0.81, 0.35, 0.03, 0.29])
    cbar = fig.colorbar(im_col, cax=cbar_ax)
    cbar.ax.set_title('Percentage', fontsize=tick_label_size, loc="left")
    cbar_ax.tick_params(axis="both", which="major", labelsize=tick_label_size)
    plt.savefig("../results/mutational_dinamics.svg")
    # plt.show()


def create_param_comparison(driv_data, bcr_data):
    fig, ax = plt.subplots(figsize=(16,8))
    graph = driv_data[["Stem cell division", "Fitness BCR-ABL"]].drop_duplicates().to_numpy()

    mut_rate = []
    for sc, bcr in graph:
        mut_rate.append(min(driv_data[(driv_data["Stem cell division"] == sc) & (driv_data["Fitness BCR-ABL"] == bcr)]["Mutation rate"]))

    driver = ax.scatter(x=graph[:,1],
                        y=graph[:,0],
                        c=mut_rate,
                        norm=matplotlib.colors.LogNorm(),
                        # cmap=plt.get_cmap("RdBu"),
                        label="\n Drivers\n",
                        edgecolor="black",
                        linewidth=1,
                        s=100
                        )

    nutral = ax.scatter(x=bcr_data["Fitness BCR-ABL"],
                        y=bcr_data["Stem cell division"],
                        norm=matplotlib.colors.LogNorm(),
                        label="\n Only BCR-ABL\n",
                        marker = 'D',
                        facecolors = "none",
                        edgecolor="black",
                        linewidth=1,
                        s=200
                        )

    plt.legend(prop={'size': 12})
    leg = ax.get_legend()
    leg.legendHandles[0].set_edgecolor('black')
    leg.legendHandles[0].set_facecolor('none')
    leg.legendHandles[0]._sizes = [200]
    cbar = plt.colorbar(driver)
    cbar.ax.set_title('Mutation rate', size=15)
    ax.set_ylabel('Stem cell total divisions per year',fontsize=20)
    ax.set_xlabel('BCR-ABL strenght',fontsize=20)
    fig.tight_layout()
    ax.tick_params(axis='both', which='major', labelsize=15)
    cbar.ax.tick_params(labelsize=15)
    cbar.ax.minorticks_on()
    fig.savefig("../results/new_relationship.svg", format="svg", transparent=True)


if __name__ == "__main__":
    if not os.path.exists("../results"):
        os.mkdir("../results")
    create_time_plot(3180, 80)
    driv_data = import_sum_dataset("../data/filtered_full_final_driver.csv")
    only_bcr_data = import_sum_dataset("../data/filtered_full_final_no_mut.csv")

    driv_phase_filter = phase_filter(driv_data)
    bcr_phase_filter = phase_filter(only_bcr_data)

    create_param_comparison(driv_phase_filter, bcr_phase_filter)

    image = plot_medians(driv_data)
    image.write_image("../results/median_durations.eps", engine="kaleido")
    create_mutation_profile(driv_data)
