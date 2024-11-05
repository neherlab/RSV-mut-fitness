"""Computing mutation rates"""

import pandas as pd
import numpy as np
from modules import glm


# Nucleotide alphabet and mutations
letters = ["A", "C", "G", "T"]
mut_types = ["AC", "AG", "AT", "CA", "CG", "CT", "GA", "GC", "GT", "TA", "TC", "TG"]


#
class Rates:
    # Initialize all rates to zeros
    def __init__(self, lb1=13467, lb2=21562):
        unpaired = [0, 1]
        # light_switch = [l_min, l_max]
        light_switch = [True, False]

        unique_muts = []

        for m in mut_types:
            for x_l in letters:
                for x_r in letters:
                    for p in unpaired:
                        if m in ["CT", "AT", "CG", "GC"]:
                            for l in light_switch:
                                unique_muts.append([m, x_l + m[0] + x_r, p, l, 0.0])
                        else:
                            unique_muts.append([m, x_l + m[0] + x_r, p, False, 0.0])

        self.rates = pd.DataFrame(
            unique_muts,
            columns=[
                "mut_type",
                "motif",
                "unpaired",
                "nt_site_before_boundary",
                "rate",
            ],
        )
        self.lb1 = lb1
        self.lb2 = lb2

    # Populate rates according to inferred GLM
    def populate_rates(self, count_df):
        # Infer GLM
        general_linear_model = glm.GeneralLinearModel(
            ["global_context", "rna_structure", "local_context"]
        )
        general_linear_model.train(count_df)
        s_max = count_df.nt_site.max()

        rates = self.rates
        rates["nt_site"] = rates["nt_site_before_boundary"].apply(
            lambda x: 0 if x else s_max
        )

        # Populate rates according to GLM
        for m in mut_types:
            rates.loc[rates["mut_type"] == m, "rate"] = (
                np.exp(
                    general_linear_model.create_data_matrix(
                        rates[rates["mut_type"] == m], m
                    )
                    @ general_linear_model.W[m]
                )
                - 0.5
            )

        rates.drop(columns=["nt_site"], inplace=True)

        # Rescale counts by total number of mutations and number of synonymous sites
        tot_mut = count_df.actual_count.sum()
        n_ss = count_df.shape[0]

        rates.rate *= n_ss / tot_mut

        # Add column `condition` with tuple summary of mutation conditions
        cond_cols = ["mut_type", "motif", "unpaired", "nt_site_before_boundary"]
        rates["condition"] = rates[cond_cols].apply(tuple, axis=1)

    def predicted_counts(self, count_df_syn):
        T = self.evol_time(count_df_syn)
        gen_comp = self.genome_composition(count_df_syn)

        assert self.rates.shape[0] == len(gen_comp)

        new_cols = {"predicted_count": (T * self.rates.rate), "cond_count": gen_comp}
        self.rates = self.rates.assign(**new_cols)

    # Computing condition specific predicted counts for all rows of clade counts table
    def predicted_counts_by_clade(self, df):
        # Look-up dictionary condition -> predicted counts
        count_dict = self.rates.set_index("condition")["predicted_count"].to_dict()

        # List of row-wise conditions
        cond_list = np.column_stack(
            [df["mut_type"], df["motif"], df["unpaired"], df["nt_site_before_boundary"]]
        )

        # Vector of predicted counts
        pred_count = np.array(list(map(lambda x: count_dict[tuple(x)], cond_list)))

        return pred_count

    def genome_composition(self, count_df_syn):
        gen_comp = np.array(
            self.rates.apply(
                lambda x: sum(
                    (x.mut_type == count_df_syn.mut_type)
                    & (x.motif == count_df_syn.motif)
                    & (x.unpaired == count_df_syn.unpaired)
                    & (
                        x.nt_site_before_boundary
                        == count_df_syn.nt_site_before_boundary
                    )
                ),
                axis=1,
            )
        )

        return gen_comp

    def evol_time(self, count_df_syn):
        gen_comp = self.genome_composition(count_df_syn)

        N_mut = sum(count_df_syn.actual_count)

        T = N_mut / (gen_comp @ np.array(self.rates.rate))

        return T

    def residual_variance(self, count_df, tau):
        self.rates["residual"] = np.zeros(self.rates.shape[0], float)

        idx = self.rates[self.rates.cond_count != 0].index

        emp_counts = (
            count_df.groupby(
                ["mut_type", "motif", "unpaired", "nt_site_before_boundary"]
            )
            .apply(lambda x: x.actual_count.to_numpy())
            .to_list()
        )

        assert len(idx) == len(emp_counts)

        self.rates.loc[idx, "residual"] = list(
            map(
                lambda x, y: np.mean((np.log(x + 0.5) - np.log(y + 0.5)) ** 2),
                emp_counts,
                self.rates.loc[idx].predicted_count.to_numpy(),
            )
        )

        av_n_cond = np.mean(self.rates.cond_count)

        self.rates["residual"] = (
            self.rates.groupby("mut_type")
            .apply(
                lambda x: np.exp(-x.cond_count / av_n_cond) * tau[x.name]
                + (1 - np.exp(-x.cond_count / av_n_cond)) * x.residual
            )
            .to_list()
        )

    # Populate array with condition specific residuals for each row of clade-wise table of counts
    def residual_by_clade(self, df):
        # Look-up dictionary condition -> residuals
        res_dict = self.rates.set_index("condition")["residual"].to_dict()

        # List of row-wise conditions
        cond_list = np.column_stack(
            [df["mut_type"], df["motif"], df["unpaired"], df["nt_site_before_boundary"]]
        )

        # Vector of residuals counts
        res = np.array(list(map(lambda x: res_dict[tuple(x)], cond_list)))

        return res


def add_predicted_count(train_df, count_df, clades):
    rate = Rates()

    rate.populate_rates(train_df)

    for c in clades:
        count_clade = count_df.loc[count_df.clade == c, :].copy()
        count_clade_syn = count_clade.loc[count_clade.mut_class == "synonymous"].copy()

        rate.predicted_counts(count_clade_syn)
        pred_count_clade = rate.predicted_counts_by_clade(count_clade)
        count_clade_syn["predicted_count"] = rate.predicted_counts_by_clade(
            count_clade_syn
        )

        tau = count_clade_syn.groupby("mut_type").apply(
            lambda x: np.mean(
                (np.log(x.actual_count + 0.5) - np.log(x.predicted_count + 0.5)) ** 2
            )
        )

        rate.residual_variance(count_clade_syn, tau)
        residual_clade = rate.residual_by_clade(count_clade)

        count_df.loc[count_clade.index, "predicted_count"] = pred_count_clade
        count_df.loc[count_clade.index, "tau_squared"] = residual_clade


# Add predicted counts to `count_df` for all clades
def add_predicted_count_all_clades(train_pre_o, train_o, count_df):
    clades = count_df.clade.unique()
    clades_pre_o = clades[0:11]
    clades_o = clades[11:]

    add_predicted_count(train_pre_o, count_df, clades_pre_o)
    add_predicted_count(train_o, count_df, clades_o)
