"""
Credits:
- Freek de Kreek, 2022

Description:
Determine mutation rate with the use of a linear mixed model, a pairwise comparison with metadata is used as input,
the pairwise comparison of the same patient and different collection date is used as input for the model, the cloud of
diversity is determined with pairwise comparisons of the same patient with the same collection data as input.

Statistics: Shapiro-Wilk Test, D'Agostino-Pearson K2 Test, White’s Lagrange Multiplier Test for Heteroscedasticity

Figures:
- scatterplot of genetic distance in days
- cloud of diversity
- normality plots of residuals linear mixed model:
    - KDE plot
    - Q-Q plot
    - RVF plot
    - boxplot by group
- random effects forest plot

Changes:
@date-of-change, description of added feature/edit

"""

# Import packages here
import argparse
import timeit
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
import statsmodels.formula.api as smf
import statsmodels.api as sm
from statsmodels.stats.diagnostic import het_white
import locale
import statsmodels.api as sm

locale.setlocale(locale.LC_TIME, 'nl_NL')

# Initialize parser, for all arguments given in command-line
parser = argparse.ArgumentParser()
# Main arguments
parser.add_argument("-i", "--Input",
                    help="Provide input file name of the distance matrix in '.xlsx' format, see help for other input "
                         "options", required=False)
parser.add_argument("-o", "--Output", help="Provide path to the output folder", required=False)

# Optional arguments
parser.add_argument("-v", "--Verbose", help="When this argument is given, it will turn on debug, this will print A LOT",
                    action="store_true")

# Variables
args = parser.parse_args()
main_input = args.Input
main_output = args.Output
verbose = args.Verbose
start = timeit.default_timer()
all_excpt_st131 = False
duplicates = False
no_outliers = True
debug = False
transform_data = False
glmm_mode = False

# Input and output
genetic_distance = 'wgMLST_stable_seqsphere'
input_study = 'riethorst'
main_input = f'./{input_study}/epi_{genetic_distance}_metadata.xlsx'
main_output = f'./{input_study}/results/{genetic_distance}/'

if transform_data:
    main_output = f'{main_output}transformed_'
if glmm_mode:
    main_output = f'{main_output}glmm_'
print('RESULTS:\n', file=open(f'{main_output}output.txt', 'w'))
######################################## Functions #####################################################################


def get_mutation_rate(pws_comp):
    """
    genetic difference / datetime difference in months * 6
    :return:
    """
    if no_outliers:
        pws_comp = pws_comp[pws_comp[genetic_distance] < 100].copy()  # tijdelijk voor outliers BRMO_surveilance, weghalen ivm gekoppeld aan genetic distance

    # calculate the difference in days between the two dates
    pws_comp['Collection_date_query'] = pd.to_datetime(pws_comp['Collection_date_query'], format='%d-%b-%Y')
    pws_comp['Collection_date_reference'] = pd.to_datetime(pws_comp['Collection_date_reference'], format='%d-%b-%Y')
    pws_comp['day_diff'] = abs((pws_comp['Collection_date_query'] - pws_comp['Collection_date_reference']).dt.days)

    # make selection of same patient and same patient at the same date
    df_same_pt = pws_comp[(pws_comp['patient_id_query'] == pws_comp['patient_id_reference']) & (pws_comp['Collection_date_query'] != pws_comp['Collection_date_reference'])].copy()
    df_same_pt_date = pws_comp[(pws_comp['patient_id_query'] == pws_comp['patient_id_reference']) & (pws_comp['Collection_date_query'] == pws_comp['Collection_date_reference'])].copy()

    df_same_pt_grouped = df_same_pt.groupby('patient_id_query').agg('count')
    print(df_same_pt[['day_diff', genetic_distance]].describe(), file=open(f'{main_output}output.txt', 'a'))

    if duplicates:
        df_same_pt = df_same_pt.drop_duplicates(subset=['patient_id_query']).copy()
        df_same_pt_date = df_same_pt_date.drop_duplicates(subset=['patient_id_query']).copy()

    """
    Linear mixed model
    Fit the linear mixed model
    """
    if glmm_mode:
        # Fit a poisson distributed generalized linear mixed model
        X = sm.add_constant(df_same_pt[['day_diff']])
        glmm_mutation_rate = sm.GLM(df_same_pt[genetic_distance], X,
                                      groups=df_same_pt['patient_id_query'],
                                      family=sm.families.NegativeBinomial(link=sm.families.links.log())).fit()

        # Back-transform the estimates and standard errors
        coef = np.exp(glmm_mutation_rate.params)
        se = np.exp(glmm_mutation_rate.bse)
        ci = np.exp(glmm_mutation_rate.conf_int())

        # Create a new summary table with the back-transformed estimates and standard errors
        ci_cols = [f'ci_{i + 1}' for i in range(ci.shape[1])]
        summary_table = pd.DataFrame({'coef': coef, 'se': se, **dict(zip(ci_cols, ci.T))})

        # Print the summary table
        print(summary_table, file=open(f'{main_output}output.txt', 'a'))
        print(glmm_mutation_rate.summary(), file=open(f'{main_output}output.txt', 'a'))
    else:
        if transform_data:
            df_same_pt[genetic_distance] = df_same_pt[genetic_distance].apply(lambda x: x**0.5)

        # Fit a linear mixed model
        lmm_mutation_rate = smf.mixedlm(f'{genetic_distance} ~ day_diff', data=df_same_pt, groups=df_same_pt['patient_id_query']).fit()

        # Print the summary of the model
        print(lmm_mutation_rate.summary(), file=open(f'{main_output}output.txt', 'a'))

    """ 
    Start plotting figures
    Create figure and axis objects with subplots()
    """
    #
    # Plot the random effects
    #

    # extract the random effects
    if glmm_mode:
        pass
    else:
        random_effects = lmm_mutation_rate.random_effects

        # create a dataframe with the random effects and their standard errors
        df_random_effects = pd.DataFrame({
            'Random Effect': random_effects.keys(),
            'Estimate': [re['Group'] for re in random_effects.values()],
            'SE': [re['Group'].std() for re in random_effects.values()]
        })

        fig, ax = plt.subplots()
        ax.errorbar(df_random_effects['Estimate'], np.arange(len(df_random_effects)), xerr=df_random_effects['SE'], fmt='o', capsize=3)
        ax.set_yticks(np.arange(len(df_random_effects)))
        ax.set_yticklabels(df_random_effects['Random Effect'])
        ax.axvline(0, linestyle='dashed', color='black')
        ax.set_xlabel(genetic_distance)
        ax.set_ylabel('Random Effect')
        ax.set_title('Random Effects Forest Plot')
        plt.tight_layout()
        plt.savefig(f'{main_output}random_effects_plot.png')
        plt.clf()

    #
    # Plot the interaction between day_diff and patient_id_query
    #
    
    fig, axis = plt.subplots(2, 1, figsize=[19.2, 14.4])
    ax1 = axis[0]
    ax2 = axis[1]

    # determine cloud of diversity
    diversity_list = [row[genetic_distance] for i, row in df_same_pt_date.iterrows()]
    diversity_list.sort()
    n_pws_comp_list = [i + 1 for i, x in enumerate(diversity_list)]

    div_per95_x = round(len(diversity_list) * 0.95)
    try:
        div_per95_y = diversity_list[div_per95_x - 1]  # -1 due to list index starting at 0

        #
        # plot figure cloud of diversity
        #

        ax1.scatter(n_pws_comp_list, diversity_list, label='allele difference')
        ax1.hlines(div_per95_y, color='grey', linestyle='dashed', label='95th percentile', xmin=0, xmax=div_per95_x)
        ax1.vlines(div_per95_x, color='grey', linestyle='dashed', ymin=0, ymax=div_per95_y)
        ax1.text(2, div_per95_y, f'95th percentile: {int(div_per95_y)} alleles')
        ax1.set_xlabel('n patients')
        ax1.set_ylabel('n allele diff')
        # ax1.set_xlim([0, max(n_pws_comp_list)])
        # ax1.set_ylim([0, max(diversity_list)])
        ax1.set_title('cloud of diversity')
    except IndexError:
        print('no cloud of diversity plot, because there are no same patient day 0 samples', file=open(f'{main_output}output.txt', 'a'))

    #
    # scatterplot of time in days and distances in alleles
    #
    
    if debug:
        print('distance in time outliers:', file=open(f'{main_output}output.txt', 'a'))
        print(df_same_pt[df_same_pt[genetic_distance] > 200][['query', 'reference', genetic_distance]], file=open(f'{main_output}output.txt', 'a'))
        print('cloud of diversity outliers:', file=open(f'{main_output}output.txt', 'a'))
        print(df_same_pt_date[df_same_pt_date[genetic_distance] > 30][['query', 'reference', genetic_distance]], file=open(f'{main_output}output.txt', 'a'))

    ax2.scatter(df_same_pt['day_diff'], df_same_pt[genetic_distance])
    # ax2.set_xlim([0, df_same_pt['day_diff'].max()])
    # ax2.set_ylim([0, df_same_pt[genetic_distance].max()])
    ax2.set_xlabel('time in days')
    ax2.set_ylabel('allele diff')
    ax2.set_title('days and distance')

    corr = df_same_pt[genetic_distance].corr(df_same_pt['day_diff'])
    plt.text(df_same_pt['day_diff'].min(), df_same_pt[genetic_distance].max(), f'Correlation: {corr:.2f}')
    # save figure
    fig.tight_layout()
    plt.savefig(f'{main_output}mutation_rate_plot.png')
    plt.clf()

    """ 
    Repeating the linear mixed model
    """
    repetitions = False
    if repetitions:
        # Define number of subsamples and iterations
        n_subsamples = 10
        n_iterations = 100

        # Define variables for storing results
        mutation_rates_lmm = []
        std_devs_lmm = []
        cutoffs_lmm = []

        # Loop through subsamples
        for i in range(n_subsamples):
            # Create random subsample of data
            subsample = df_same_pt.sample(frac=0.5)

            # Loop through iterations
            for j in range(n_iterations):
                # Create random seed for reproducibility
                np.random.seed(j)

                # Fit linear mixed model
                model = sm.MixedLM.from_formula('wgMLST_stable_seqsphere ~ day_diff',
                                                data=subsample,
                                                groups=subsample['patient_id_query']).fit()

                # Extract mutation rate
                mutation_rate_lmm = -model.params['day_diff']

                # Append results to lists
                mutation_rates_lmm.append(mutation_rate_lmm)
                std_devs_lmm.append(model.bse['day_diff'])
                cutoffs_lmm.append(mutation_rate_lmm * 180)  # 2 * mutation rate

        print(f"Mutation rate: {np.mean(mutation_rates_lmm):.3f} +/- {np.mean(std_devs_lmm):.3f} (std dev), cutoff "
              f"{np.mean(cutoffs_lmm)} all values are the means of {n_subsamples = } and {n_iterations = }", file=open(f'{main_output}output.txt', 'a'))

    if glmm_mode:
        end_result = glmm_mutation_rate._results
    else:
        end_result = lmm_mutation_rate
    return end_result


def test_normality(model):
    # create axis for normality plots
    fig, axis = plt.subplots(4, 1, figsize=[19.2, 14.4])
    ax1 = axis[0]
    ax2 = axis[1]
    ax3 = axis[2]
    ax4 = axis[3]

    if glmm_mode:
        model_residues = model.resid_deviance
    else:
        model_residues = model.resid

    # sns.distplot(model_residues, hist=False, kde_kws={"fill": True, "lw": 1}, fit=stats.norm, ax=ax1)  # deprecated
    sns.kdeplot(model_residues, ax=ax1, fill=True)
    mu, std = stats.norm.fit(model_residues)
    xmin, xmax = ax1.get_xlim()
    x = np.linspace(xmin, xmax, 100)
    p = stats.norm.pdf(x, mu, std)
    ax1.plot(x, p, 'k', linewidth=2)

    ax1.set_title("KDE Plot of Model Residuals (Blue) and Normal Distribution (Black)")
    ax1.set_xlabel("Residuals")

    ## Q-Q PLot
    sm.qqplot(model_residues, dist=stats.norm, line='s', ax=ax2)
    ax2.set_title("Q-Q Plot")

    ## Shapir-Wilk test of normality
    labels_SW = ["Shapiro-Wilk Test:", "p-value"]
    norm_res_SW = stats.shapiro(model_residues)

    for key, val in dict(zip(labels_SW, norm_res_SW)).items():
        print(key, val, file=open(f'{main_output}output.txt', 'a'))

    ## Anderson-Darling test of normality
    labels_AD = ["Anderson-Darling Test:", "p-value"]
    norm_res_AD = stats.anderson(model_residues)

    print(f"Anderson-Darling Test: {norm_res_AD.statistic: .4f}", file=open(f'{main_output}output.txt', 'a'))

    for i in range(len(norm_res_AD.critical_values)):
        sig, crit = norm_res_AD.significance_level[i], norm_res_AD.critical_values[i]

        if norm_res_AD.statistic < norm_res_AD.critical_values[i]:
            print(
                f"At {sig}% significance,{norm_res_AD.statistic: .4f} <{norm_res_AD.critical_values[i]: .4f} data looks normal (fail to reject H0)", file=open(f'{main_output}output.txt', 'a'))
        else:
            print(
                f"At {sig}% significance,{norm_res_AD.statistic: .4f} >{norm_res_AD.critical_values[i]: .4f} data does not look normal (reject H0)", file=open(f'{main_output}output.txt', 'a'))


    ## HOMOSKEDASTICITY OF VARIANCE
    # residuals versus fitted values (RVF) plot
    sns.scatterplot(y=model_residues, x=model.fittedvalues, ax=ax3)
    ax3.axhline(y=0, linestyle='--')

    ax3.set_title("RVF Plot")
    ax3.set_xlabel("Fitted Values")
    ax3.set_ylabel("Residuals")

    ## boxplot of residuals
    sns.boxplot(x=model.model.groups, y=model_residues, ax=ax4)
    ax4.set_title("Distribution of Residuals for Weight by Litter")
    ax4.set_ylabel("Residuals")
    ax4.set_xlabel("Litter")
    plt.tight_layout()
    plt.savefig(f'{main_output}normality_plots.png')

    ## White’s Lagrange Multiplier Test for Heteroscedasticity
    het_white_res = het_white(model_residues, model.model.exog)

    labels = ["LM Statistic", "LM-Test p-value", "F-Statistic", "F-Test p-value"]

    for key, val in dict(zip(labels, het_white_res)).items():
        print(key, val, file=open(f'{main_output}output.txt', 'a'))
    return


######################################## Main ##########################################################################
# Call functions here

mutation_model = get_mutation_rate(pd.read_excel(main_input))
test_normality(mutation_model)


# End of script
stop = timeit.default_timer()  # ends timer
print("Finished succesfully, time taken: ", stop - start, file=open(f'{main_output}output.txt', 'a'))  # prints script run duration
