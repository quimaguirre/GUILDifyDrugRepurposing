import os, sys
from collections import defaultdict
import rpy2.robjects as robjects
from statsmodels.sandbox.stats.multicomp import multipletests
import scipy.stats as stats

def main():

    type_top = 'p' # p: percent / e: enrich
    top1 = 1
    top2 = 1
    calculation = 'E1' # E1 / E2 / FCG / FCF
    #selected_method = 'GObp' # GObp / GOmf / Reactome
    selected_methods = ['GObp', 'GOmf', 'Reactome']
    criteria = 'fdr_bh' # fdr_bh / bonferroni
    taxID = 9606 # human:9606, mouse:10090, rat: 10116, worm:6239, fly:7227, plant:3702, yeast:4932
    associations_dir = '/Users/quim/Dropbox/UPF/PhD/Projects/GUILDify/case_studies_JMB/data/{}'.format(taxID)
    working_dir = '/Users/quim/Dropbox/UPF/PhD/Projects/GUILDify/case_studies_JMB/NSCLC'
    go_obo_file = '/Users/quim/Documents/Projects/goatools/go-basic.obo'
    disease1 = 'cell_lung_cancer'
    disease2 = 'erlotinib'
    #working_dir = '/Users/quim/Dropbox/UPF/PhD/Projects/GUILDify/case_studies_JMB/NSCLC'
    #disease1 = 'cell_lung_cancer'
    #diseases2 = ['afatinib', 'ceritinib', 'crizotinib', 'erlotinib', 'gefitinib', 'palbociclib', 'lapatinib']

    for selected_method in selected_methods:

        dicbp, dicbp2, term2name = load_functional_terms(selected_method, associations_dir, taxID)
        #print(dicbp2)
        N = len(dicbp.keys())
        NB = len(dicbp2.keys())

        # Parse top nodes in disease 1
        disease_dir1 = os.path.join(working_dir, disease1)
        scores_file_n1 = os.path.join(disease_dir1, 'guild_scores.txt')
        ranking_n1, node_to_info_n1 = read_scores_file(scores_file_n1)
        if type_top == 'p':
            length_top_nodes_n1 = int(round(top1 * float(len(ranking_n1)) / float(100)))
            top_nodes_n1 = set(ranking_n1[0: length_top_nodes_n1])
        elif type_top == 'e':
            top_nodes_n1 = set(ranking_n1[0: top1])
        print(len(top_nodes_n1))

        # Parse top nodes in disease 2
        disease_dir2 = os.path.join(working_dir, disease2)
        scores_file_n2 = os.path.join(disease_dir2, 'guild_scores.txt')
        ranking_n2, node_to_info_n2 = read_scores_file(scores_file_n2)
        if type_top == 'p':
            length_top_nodes_n2 = int(round(top2 * float(len(ranking_n2)) / float(100)))
            top_nodes_n2 = set(ranking_n2[0: length_top_nodes_n2])
        elif type_top == 'e':
            top_nodes_n2 = set(ranking_n2[0: top2])
        print(len(top_nodes_n2))

        if calculation == 'E1':
            print(disease1)
            top_geneIDs_n1 = set()
            for user_entity in top_nodes_n1:
                if user_entity in node_to_info_n1:
                    geneID = node_to_info_n1[user_entity][5]
                    if geneID != '-':
                        top_geneIDs_n1.add(geneID)
            query = top_geneIDs_n1
            output_file = os.path.join(disease_dir1, 'enrichment_{}.{}.{}.{}.txt'.format(disease1, selected_method, criteria, type_top+str(top1)))
        elif calculation == 'E2':
            print(disease2)
            top_geneIDs_n2 = set()
            for user_entity in top_nodes_n2:
                if user_entity in node_to_info_n2:
                    geneID = node_to_info_n2[user_entity][5]
                    if geneID != '-':
                        top_geneIDs_n2.add(geneID)
            query = top_geneIDs_n2
            output_file = os.path.join(disease_dir2, 'enrichment_{}.{}.{}.{}.txt'.format(disease2, selected_method, criteria, type_top+str(top2)))
        elif calculation == 'FCG':
            print(disease1, disease2)
            intersection_top = top_nodes_n1 & top_nodes_n2
            intersection_top_geneIDs = set()
            for user_entity in intersection_top:
                if user_entity in node_to_info_n1:
                    geneID = node_to_info_n1[user_entity][5]
                    if geneID != '-':
                        intersection_top_geneIDs.add(geneID)
            output_file = os.path.join(working_dir, 'CG_{}_{}.{}.{}.txt'.format(disease1, disease2, type_top+str(top1), type_top+str(top2)))
            print(intersection_top_geneIDs)
            with open(output_file, 'w') as out_fd:
                for user_entity in intersection_top:
                    out_fd.write('{}\n'.format('\t'.join(node_to_info_n1[user_entity])))
            query = intersection_top_geneIDs
            output_file = os.path.join(working_dir, 'FCG_{}_{}.{}.{}.{}.{}.txt'.format(disease1, disease2, selected_method, criteria, type_top+str(top1), type_top+str(top2)))
        elif calculation == 'FCF':
            print(disease1, disease2)
            from scipy.stats import combine_pvalues
            functions_file_n1 = os.path.join(disease_dir1, 'enrichment_{}.{}.{}.{}.txt'.format(disease1, selected_method, criteria, type_top+str(top1)))
            functions_file_n2 = os.path.join(disease_dir2, 'enrichment_{}.{}.{}.{}.txt'.format(disease2, selected_method, criteria, type_top+str(top2)))
            term_to_values_n1 = parse_functions_file(functions_file_n1)
            term_to_values_n2 = parse_functions_file(functions_file_n2)
            common = set(term_to_values_n1.keys()) & set(term_to_values_n2.keys())
            print(common)
            common_values = []
            for term in common:
                term_name = term_to_values_n1[term][0]
                x1 = term_to_values_n1[term][1]
                x2 = term_to_values_n2[term][1]
                m = term_to_values_n1[term][2]
                pvalue1 = float(term_to_values_n1[term][3])
                pvalue2 = float(term_to_values_n2[term][3])
                chi, comb_pval = combine_pvalues([pvalue1, pvalue2], method='fisher')
                common_values.append([term, term_name, x1, x2, m, pvalue1, pvalue2, comb_pval])
            common_values = sorted(common_values, key=lambda x: x[7])
            output_file = os.path.join(working_dir, 'FCF_{}_{}.{}.{}.{}.{}.txt'.format(disease1, disease2, selected_method, criteria, type_top+str(top1), type_top+str(top2)))
            with open(output_file, 'w') as out_fd:
                for values in common_values:
                    term, term_name, x1, x2, m, pvalue1, pvalue2, comb_pval = values
                    out_fd.write('{}\t{}\t{}\t{}\t{}\t{:.1E}\t{:.1E}\t{:.1E}\n'.format(term, term_name, x1, x2, m, pvalue1, pvalue2, comb_pval))
        else:
            print('Wrong calculation.')
            sys.exit(10)

        if calculation != 'FCF':
            test_passed_f1, terms_l1, term_to_values = functional_enrichment(dicbp2, N, list(query), criteria)
            print(test_passed_f1)
            print(terms_l1)
            with open(output_file, 'w') as out_fd:
                for term, values in sorted(term_to_values.items(),key=lambda x: x[1][1],reverse=False):
                    pval, pval_corrected, x, m = values
                    if term in term2name:
                        term_name = term2name[term]
                    else:
                        term_name = '-'
                    out_fd.write('{}\t{}\t{}\t{}\t{}\t{}\n'.format(term, term_name, x, m, pval, pval_corrected))

    return


#######################
#######################
# SECONDARY FUNCTIONS #
#######################
#######################

def load_functional_terms(selected_method, associations_dir, taxID):
    #file_method = {'GOmf': 'GOmf_DSGPPI.csv', 'GObp': 'GObp_DSGPPI.csv', 'Reactome': 'Reactome_DSGPPI.csv'}
    file_method = {'GOmf': 'GOmf2gene.{}.txt'.format(taxID), 'GObp': 'GObp2gene.{}.txt'.format(taxID), 'Reactome': 'Reactome2gene.{}.txt'.format(taxID)}
    data_file = os.path.join(associations_dir, '{}'.format(file_method[selected_method]))
    dicbp = defaultdict(list)
    dicbp2 = defaultdict(list)
    term2name = {}
    fi3 = open(data_file)
    for line in fi3.readlines():
        gene = line.rstrip().split("\t")
        if len(gene) == 2:
            gobp = gene[0]
            geneid = gene[1]
        elif len(gene) == 3:
            gobp = gene[0]
            name = gene[1]
            geneid = gene[2]
            term2name[gobp] = name

        if gobp not in dicbp[geneid]:
            dicbp[geneid].append(gobp)

        if geneid not in dicbp2[gobp]:
            dicbp2[gobp].append(geneid)

    fi3.close()
    print(selected_method, 'data loaded')
    return dicbp, dicbp2, term2name

def functional_enrichment(dicbp2, N, geneset_toenrich, criteria):
    pvals = {}
    genes_enriched = {}
    go_enriched = False
    k = len(geneset_toenrich)
    terms_l = []
    test_passed_f = False
    term_to_values = {}
    for term in dicbp2.keys():
        m = len(dicbp2[term])
        xl = [y for y in dicbp2[term] if y in geneset_toenrich]
        x = len(xl)

        if x != 0:
            go_enriched = True
            xlist = []
            for i in range(x, m + 1):
                xlist.append(i)

            # calculation of the hypervalue
            dhyper = robjects.r['dhyper']
            xboh = robjects.IntVector(xlist)
            dhypervalue = dhyper(xboh, m, (N - m), k, log=False)
            # threshold of enrichment
            pvals[term] = sum(dhypervalue)
            hypervalue_scipy = sum(stats.hypergeom.pmf(xlist, N, m, k))
            genes_enriched[term] = [x, m] # Quim: addition to get genes enriched
    if go_enriched:
        pvals_values = list(pvals.values())
        terms = list(pvals.keys())
        pvals_corrected = multipletests(pvals_values, alpha=0.05, method=criteria, is_sorted=False,
                                        returnsorted=False)

        for i in range(0, len(terms)):
            if list(pvals_corrected[1])[i] < 0.05:
                # Quim: addition to get the p-values and genes enriched
                pval = pvals_values[i]
                pval_corrected = list(pvals_corrected[1])[i]
                term = terms[i]
                x, m = genes_enriched[term]
                term_to_values[term] = [pval, pval_corrected, x, m]
                ####################################
                test_passed_f = True
                terms_l.append(terms[i])
    return test_passed_f, terms_l, term_to_values

def read_scores_file(scores_file):
    ranking_nodes = []
    ranking_nodes_to_info = {}
    with open(scores_file, 'r') as scores_fd:
        first_line = scores_fd.readline()
        for line in scores_fd:
            fields = line.strip().split('\t')
            biana_id, uniprot, gene, seed, description, geneID, eq_entries, score, degree = fields
            ranking_nodes.append(biana_id)
            ranking_nodes_to_info[biana_id] = [biana_id, uniprot, gene, seed, description, geneID, eq_entries, score, degree]
    return ranking_nodes, ranking_nodes_to_info

def parse_functions_file(functions_file, cutoff = 0.05):
    term_to_values = {}
    with open(functions_file, 'r') as functions_fd:
        for line in functions_fd:
            term, term_name, x, m, pval, pval_corrected = line.strip().split('\t')
            if float(pval_corrected) < cutoff:
                term_to_values[term] = [term_name, x, m, pval_corrected]
    return term_to_values

if  __name__ == "__main__":
    main()
