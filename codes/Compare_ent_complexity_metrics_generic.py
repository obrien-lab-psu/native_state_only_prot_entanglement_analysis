import sys, os, re, time, logging
from sklearn.utils import shuffle
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LassoCV, Lasso
from sklearn.datasets import make_regression
from sklearn.model_selection import StratifiedKFold, KFold, cross_validate, GridSearchCV
from sklearn.metrics import mean_squared_error
from sklearn.metrics import confusion_matrix
from sklearn.metrics import accuracy_score, balanced_accuracy_score, average_precision_score,f1_score,recall_score,precision_score,roc_auc_score
from sklearn.model_selection import train_test_split
from sklearn.linear_model import LogisticRegression, LogisticRegressionCV
from sklearn import metrics, preprocessing
from sklearn.preprocessing import LabelEncoder
import argparse
import numpy as np
import pandas as pd
from glob import glob
import pickle
from itertools import product, combinations
from scipy.stats import permutation_test, ttest_ind, false_discovery_control, mode, bootstrap, mannwhitneyu
import matplotlib.pyplot as plt
#pd.set_option('display.max_rows', 500)

class DataAnalysis:
    """
    This class is meant to calculate the distribution of various entanglement complexity metrics between Group1 and Group2 gene list

    """

    def __init__(self, outpath, uent_files):
        """
        Initializing the DataAnalysis object and ensure atleast the top level outpath is present and makes it if not. 
        """

        # Make outpath if it doesnt exists
        self.outpath = outpath
        self.DataAnalysisOutpath = os.path.join(self.outpath, 'DataAnalysisOutput/')
        if not os.path.exists(self.DataAnalysisOutpath):
            os.makedirs(self.DataAnalysisOutpath)
            print(f'Made directory: {self.DataAnalysisOutpath}')

        self.DataAnalysisDistPlotsOutpath = os.path.join(self.outpath, 'DataAnalysisOutput/Dist_plots/')
        if not os.path.exists(self.DataAnalysisDistPlotsOutpath):
            os.makedirs(self.DataAnalysisDistPlotsOutpath)
            print(f'Made directory: {self.DataAnalysisDistPlotsOutpath}')

        self.three_to_one_letter = {
        'ALA': 'A', 'ARG': 'R', 'ASN': 'N', 'ASP': 'D', 'CYS': 'C',
        'GLU': 'E', 'GLN': 'Q', 'GLY': 'G', 'HIS': 'H', 'ILE': 'I',
        'LEU': 'L', 'LYS': 'K', 'MET': 'M', 'MSE': 'M', 'PHE': 'F', 
        'PRO': 'P', 'SER': 'S', 'THR': 'T', 'TRP': 'W', 'TYR': 'Y', 
        'VAL': 'V'}
        
        self.uent_files = glob(os.path.join(uent_files, '*'))

        self.keys = ['Gn', 'N_term_thread', 'Gc', 'C_term_thread',
        'loopsize', 'num_zipper_nc', 'perc_bb_loop',
        'num_loop_contacting_res', 'num_cross_nearest_neighbors',
        'ent_coverage', 'min_N_prot_depth_left', 'min_N_thread_depth_left',
        'min_N_thread_slippage_left', 'min_C_prot_depth_right',
        'min_C_thread_depth_right', 'min_C_thread_slippage_right', 'ACO', 'RCO']


    #################################################################################################################
    def load_files(self, mask):
        """
        Method to load unique entanglement files and keep only certain columns for analysis
        'Gn', 
        'N_term_thread', 
        'Gc', 
        'C_term_thread',
        'loopsize', 
        'num_zipper_nc', 
        'perc_bb_loop',
        'num_loop_contacting_res', 
        'num_cross_nearest_neighbors',
        'ent_coverage', 
        'min_N_prot_depth_left', 
        'min_N_thread_depth_left',
        'min_N_thread_slippage_left', 
        'min_C_prot_depth_right',
        'min_C_thread_depth_right', 
        'min_C_thread_slippage_right', 
        'prot_size'

        """
        dfs = []
        for gene in mask:
            print(gene)
            gene_uent_file = [f for f in self.uent_files if gene in f][0]

            gene_uent = pd.read_csv(gene_uent_file, sep='|')
            #print(gene, gene_uent_file)
            #print(gene_uent)
            num_uent = len(gene_uent)

            if num_uent != 0:
                gene_uent['gene'] = gene
                dfs += [gene_uent]

        uent_df = pd.concat(dfs) 
        return uent_df


    def DistStats(self, df, keys, dist_tag, n_resamples=1000, alpha=0.05):
        """
        Calculate various parameters of the data distributions
        (1) mean, median, mode
        (2) confidence intervals for mean, median, and mode
        * change the self.keys parameters to suit your needs
        """
        df_copy = df.copy()
        df_copy = df_copy[keys]
        print(f'df_copy:\n{df_copy}')

        results = {'metric':[], 'mean':[], 'mean_lb':[], 'mean_ub':[], 'median':[], 'median_lb':[], 'median_ub':[]}

        def calculate_bootstrap_ci(data, statistic_func):
            if statistic_func == 'mean':
                res = bootstrap((data,), np.mean, confidence_level=1-alpha, n_resamples=n_resamples)
                return res.confidence_interval.low, res.confidence_interval.high

            elif statistic_func == 'median':
                medians = []
                for b in range(n_resamples):
                    boot = np.random.choice(data, replace=True)
                    medians += [np.median(boot)]
                lb = np.percentile(medians, 2.5)
                ub = np.percentile(medians, 97.5)
                return lb, ub


        for column in df_copy.columns:
            print(f'column: {column}')
            col_data = df_copy[column].dropna().values

            ## drop 0s if column == C_term_thread or N_term_thread
            if column in ['N_term_thread', 'C_term_thread']:
                col_data = col_data[np.where(col_data != 0)]
            print(col_data)

            mean_val = np.mean(col_data)
            median_val = np.median(col_data)

            mean_ci = calculate_bootstrap_ci(col_data, 'mean')
            print('Mean', mean_val, mean_ci)
            median_ci = calculate_bootstrap_ci(col_data, 'median')
            print('Median', median_val, median_ci)

            results['metric'] += [column]
            results['mean'] += [mean_val]
            results['mean_lb'] += [mean_ci[0]]
            results['mean_ub'] += [mean_ci[1]]
            results['median'] += [median_val]
            results['median_lb'] += [median_ci[0]]
            results['median_ub'] += [median_ci[1]]

            # plot histogram
            plot_filename = f'{self.DataAnalysisDistPlotsOutpath}{dist_tag}_{column}.png'
            plt.hist(col_data, bins=100, color='blue', edgecolor='black', density=True)  # 100 bins
            plt.xlabel(column)
            plt.ylabel('PDF')
            plt.title(dist_tag)
            plt.savefig(plot_filename)
            print(f'SAVED: {plot_filename}')
            plt.close()

        return pd.DataFrame(results)


    def Permutation(self, df1_full, df2_full, keys, n_resamples=10000):
        """
        For the columns in the two dataframes calculate the pairwise pvalue by permutation. 
        """

        df1 = df1_full.copy()
        df2 = df2_full.copy()
        df1 = df1[keys]
        df2 = df2[keys]

        if not all(df1.columns == df2.columns):
            raise ValueError("Both DataFrames must have identical column names")

        results = []

        for column in df1.columns:
            data1 = df1[column].dropna().values
            data2 = df2[column].dropna().values

            # Define the statistic function for the permutation test
            def statistic(x, y, axis):
                return np.mean(x, axis=axis) - np.mean(y, axis=axis)

            # Perform the permutation test
            #res = permutation_test((data1, data2), statistic, vectorized=True, n_resamples=n_resamples)
            stat = statistic(data1, data2, 0)
            if stat > 0:
                res = mannwhitneyu(data1, data2, alternative='greater')
            elif stat < 0:
                res = mannwhitneyu(data1, data2, alternative='less')
            else:
                res = mannwhitneyu(data1, data2)

            #results[column] = res.pvalue
            results += [res.pvalue]

            print(f'column: {res.pvalue}')

        return results

        
    def LassoRegression(self, df, X_keys, Y_key):

        ### Get input data 
        X = df[X_keys]
        X = X.fillna(0)
        # Initialize the LabelEncoder and make sure y categories are integer incoded
        y = df[Y_key]
        label_encoder = LabelEncoder()
        y_encoded = label_encoder.fit_transform(y)

        scaler = preprocessing.StandardScaler()
        std_scaled_df = scaler.fit_transform(X)
        std_scaled_df = pd.DataFrame(std_scaled_df, columns=X_keys)
        X = std_scaled_df
        X = X.values
        y = y.values
        #print(f'X:\n{X}')
        #print(f'y:\n{y}')

        logistic_regression =  LogisticRegression(penalty='l1', solver='liblinear')
        Cs = np.linspace(0.00001, 10, 1000)

        fit_data = {'C':[], 'fold':[], 'balanced_accuracy':[], 'accuracy':[]}
        for C in Cs:

            print(f'{"#"*100}\nTESTING C: {C}')

            for col in X_keys:
                if col not in fit_data:
                    fit_data[col] = []

            ## make folds and fit model
            skf = StratifiedKFold(n_splits=5, shuffle=True)
            #skf = StratifiedKFold(n_splits=5)

            for i, (train_index, test_index) in enumerate(skf.split(X, y)):
                #print(f"Fold {i}:")
                #print(f"  Train: index={train_index} {len(train_index)}")
                #print(f"  Test:  index={test_index} {len(test_index)}")

                X_train = X[train_index]
                y_train = y[train_index]
                X_test = X[test_index]
                y_test = y[test_index]

                # Get features for optimal regularization ceof
                logistic_regression =  LogisticRegression(penalty='l1', solver='liblinear', C=C)
                logistic_regression.fit(X_train, y_train)

                coefs = logistic_regression.coef_[0].tolist()
                #print(f'coefs: {coefs}')

                # Predict on the testing data
                y_pred = logistic_regression.predict(X_test)
                # Calculate balanced accuracy
                balanced_accuracy = balanced_accuracy_score(y_test, y_pred)
                accuracy = accuracy_score(y_test, y_pred)

                fit_data['C'] += [C]
                fit_data['fold'] += [i]
                fit_data['balanced_accuracy'] += [balanced_accuracy]
                fit_data['accuracy'] += [accuracy]
                for col_i, col in enumerate(X_keys):
                    fit_data[col] += [coefs[col_i]]
                    
        fit_data = pd.DataFrame(fit_data)
        print(f'fit_data:\n{fit_data}')
        return fit_data

    def Plot_Lasso(self, df, keys, outfile):
        
        #        C  fold  balanced_accuracy  accuracy  
        # Create subplots
        fig, axes = plt.subplots(1, 3, figsize=(12, 4))
        X = []
        num_nonzero = []
        num_robust_nonzero = []
        BA = []
        for C, C_df in df.groupby('C'):
            #print(C_df)

            num_nonzero_coef = []
            num_robust_nonzero_coef = []
            for key in keys:
                coefs = C_df[key].values
                all_nonzero = np.all(coefs != 0)
                same_sign = np.all(coefs >= 0) or np.all(coefs <= 0)
                #print(C, key, coefs, all_nonzero, same_sign)
                num_nonzero_coef += [all_nonzero]
                num_robust_nonzero_coef += [same_sign]

            X += [C]
            num_nonzero += [np.sum(num_nonzero_coef)]
            num_robust_nonzero += [np.sum(num_robust_nonzero_coef)]
            BA += [np.mean(C_df['balanced_accuracy'].values)]

        axes[0].plot(X, num_nonzero)
        axes[0].set_ylabel('# non-zero ceof.')
        axes[0].set_xlabel('inverse regularization strength')
        axes[0].set_ylim(0, 20)
        axes[1].plot(X, num_robust_nonzero)
        axes[1].set_ylabel('# non-zero & robust ceof.')
        axes[1].set_xlabel('inverse regularization strength')
        axes[1].set_ylim(0, 20)
        axes[2].plot(X, BA)
        axes[2].set_ylabel('Balanced Accuracy')
        axes[2].set_xlabel('inverse regularization strength')
        axes[2].set_ylim(0, 1)

        plt.tight_layout()
        plt.savefig(outfile)
        print(f'SAVED: {outfile}')

#################################################################################################################
def main():
    """
    Main function to control workflow. 
    """

    # Parse the user supplied arguments
    parser = argparse.ArgumentParser(description="Process user specified arguments")
    parser.add_argument("-g1", "--Group1_gene_list", type=str, required=True, help=f"path to Group1 gene list used for mask")
    parser.add_argument("-g2", "--Group2_gene_list", type=str, required=True, help=f"path to Group2 gene list used for mask")
    parser.add_argument("-l", "--log_file", type=str, required=True, help="Path to logging file")
    parser.add_argument("-e", "--uent_files", type=str, required=True, help="path to unique entanglement files")
    parser.add_argument("-o", "--outpath", type=str, required=True, help="path to output directory. will be made if doesnt exist")
    parser.add_argument("-p", "--num_permute", type=int, required=True, help="Number of permutations")
    args = parser.parse_args()

    uent_files = args.uent_files
    log_file = args.log_file
    outpath = args.outpath
    Group1_gene_list = args.Group1_gene_list
    Group2_gene_list = args.Group2_gene_list
    num_permute = args.num_permute

    # Setup logging configuration
    logging.basicConfig(filename=log_file, level=logging.INFO, format='%(asctime)s %(message)s') 
    logging.info(f'{"#"*100}\nNEW RUN')

    # Make outpath if it doesnt exists
    if not os.path.exists(outpath):
        os.makedirs(outpath)
        print(f'Made directory: {outpath}')

    # Initalize the DataAnalysis class object
    Analyzer = DataAnalysis(outpath, uent_files)
    print(f'Analyzer: {Analyzer}')

    ## Get the Group1  gene data and stats
    Group1_gene_mask = np.loadtxt(Group1_gene_list, dtype=str)
    print(f'Group1_gene_mask: {Group1_gene_mask} {len(Group1_gene_mask)}')

    Group1_ent_data = Analyzer.load_files(Group1_gene_mask)
    Group1_ent_data['Gn'] = Group1_ent_data['Gn'].abs()
    Group1_ent_data['Gn'] = Group1_ent_data['Gn'].apply(lambda x: np.nan if x < 0.6 else x)
    Group1_ent_data['Gc'] = Group1_ent_data['Gc'].abs()
    Group1_ent_data['Gc'] = Group1_ent_data['Gc'].apply(lambda x: np.nan if x < 0.6 else x)
    Group1_ent_data['Group'] = 1
    print(f'Group1_ent_data: {Group1_ent_data}')
    Group1_ent_data_outfile = f'{Analyzer.outpath}Group1_uent_data.csv'
    Group1_ent_data.to_csv(Group1_ent_data_outfile, sep='|', index=False)
    print(f'SAVED: {Group1_ent_data_outfile}')

    Group1_stats = Analyzer.DistStats(Group1_ent_data, Analyzer.keys, 'Group1', n_resamples=num_permute)
    print(f'Group1_stats:\n{Group1_stats}')
    Group1_stats_data_outfile = f'{Analyzer.outpath}Group1_stats_uent_data.csv'
    Group1_stats.to_csv(Group1_stats_data_outfile, sep='|', index=False)
    print(f'SAVED: {Group1_stats_data_outfile}')


    ## Get the Group 2 gene data and stats
    Group2_gene_mask = np.loadtxt(Group2_gene_list, dtype=str)
    print(f'Group2_gene_mask: {Group2_gene_mask} {len(Group2_gene_mask)}')

    Group2_ent_data = Analyzer.load_files(Group2_gene_mask)
    Group2_ent_data['Gn'] = Group2_ent_data['Gn'].abs()
    Group2_ent_data['Gn'] = Group2_ent_data['Gn'].apply(lambda x: np.nan if x < 0.6 else x)
    Group2_ent_data['Gc'] = Group2_ent_data['Gc'].abs()
    Group2_ent_data['Gc'] = Group2_ent_data['Gc'].apply(lambda x: np.nan if x < 0.6 else x)
    Group2_ent_data['Group'] = 2
    print(f'Group2_ent_data: {Group2_ent_data}')
    Group2_ent_data_outfile = f'{Analyzer.outpath}Group2_uent_data.csv'
    Group2_ent_data.to_csv(Group2_ent_data_outfile, sep='|', index=False)
    print(f'SAVED: {Group2_ent_data_outfile}')

    Group2_stats = Analyzer.DistStats(Group2_ent_data, Analyzer.keys, 'Group2', n_resamples=num_permute)
    print(f'Group2_stats:\n{Group2_stats}')
    Group2_stats_data_outfile = f'{Analyzer.outpath}Group2_stats_uent_data.csv'
    Group2_stats.to_csv(Group2_stats_data_outfile, sep='|', index=False)
    print(f'SAVED: {Group2_stats_data_outfile}')


    ## Compare the entanglement complexity between Group1 and Group2 data
    Group1VGroup2_pvalues = Analyzer.Permutation(Group1_ent_data, Group2_ent_data, Analyzer.keys, n_resamples=num_permute)
    print(f'Group1VGroup2_pvalues:\n{Group1VGroup2_pvalues}')


    # Merge the DataFrames along axis=1
    metrics = Group1_stats['metric']
    Group1_stats.drop(columns=['metric'])
    Group2_stats.drop(columns=['metric'])

    Group1_stats = Group1_stats.add_prefix('Group1_')
    Group2_stats = Group2_stats.add_prefix('Group2_')

    merged_df = pd.concat([metrics, Group1_stats, Group2_stats], axis=1)
    merged_df['pvalues'] = Group1VGroup2_pvalues
    print(f'merged_df:\n{merged_df}')
    merged_data_outfile = f'{Analyzer.outpath}merged_stats_uent_data.csv'
    merged_df.to_csv(merged_data_outfile, sep='|', index=False)
    print(f'SAVED: {merged_data_outfile}')

    # make the combined raw ent dataframes
    combined_df = pd.concat([Group1_ent_data, Group2_ent_data])
    print(f'combined_df:\n{combined_df}')
    combined_data_outfile = f'{Analyzer.outpath}combined_uent_data.csv'
    combined_df.to_csv(combined_data_outfile, sep='|', index=False)
    print(f'SAVED: {combined_data_outfile}')

    Lasso_results = Analyzer.LassoRegression(combined_df, Analyzer.keys, 'Group')
    print(f'Lasso_results:\n{Lasso_results}')
    Lasso_results_outfile = f'{Analyzer.outpath}Lasso_results.csv'
    Lasso_results.to_csv(Lasso_results_outfile, sep='|', index=False)
    print(f'SAVED: {Lasso_results_outfile}')

    Lasso_results_plot_outfile = f'{Analyzer.outpath}Lasso_results.png'
    Analyzer.Plot_Lasso(Lasso_results, Analyzer.keys, Lasso_results_plot_outfile)

if __name__ == "__main__":
    main()

print('NORMAL TERMINATION')
