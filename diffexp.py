
import mne
import pandas as pd
import json

#import seaborn as sns
#import matplotlib.pyplot as plt 
import pandas as pd
import numpy as np
#reading the json file

import numpy as np

def min_acc(array, key):
    min_val = array[0][key]
    minimized_arr = []
    for item in array:
        if item[key] > min_val:
            item[key] = min_val
            minimized_arr.append(item)
        else:
            minimized_arr.append(item)
            min_val = item[key]
    return minimized_arr

def _ecdf(x):
    """No frills empirical cdf used in fdrcorrection."""
    nobs = len(x)
    return np.arange(1, nobs + 1) / float(nobs)

def modded_fdr_correction(items, pval_key, alpha=0.05, method='indep'):
    items.sort(key=lambda x: x[pval_key])

    #items_size = len(items)

    if method in ['i', 'indep', 'p', 'poscorr']:
        factors = list(_ecdf(items))
        #factors.reverse()
        ecdffactors = factors
    elif method in ['n', 'negcorr']:
        cm = np.sum(1. / np.arange(1, len(items) + 1))
        ecdffactors = _ecdf(items) / cm
    else:
        raise ValueError("Method should be 'indep' and 'negcorr'")

    with open("ecdf_modded.json","w") as file:
      file.write(json.dumps(list(ecdffactors)))

    rejectmax = 0
    for idx in range(0, len(items)):
        item = items[idx]
        ecdffactor = ecdffactors[idx]
        #print(item[pval_key], ecdffactor, alpha)
        reject = item[pval_key] < (ecdffactor * alpha)
        item['fdr_reject'] = bool(reject)
        if reject:
            rejectmax += 1
    
    for idx in range(0, rejectmax):
        item = items[idx]
        item['fdr_reject'] = bool(True)

    with open("pvals_sorted_fdr_modded.json","w") as file:
      file.write(json.dumps(list(items)))

    raw_min_pval = 1.0
    for idx in range(0, len(items)):
        item = items[idx]
        ecdffactor = ecdffactors[idx]

        item['fdr_pval'] = item[pval_key] / ecdffactor
        #raw_min_pval = min(raw_min_pval, (item[pval_key] / ecdffactor))

    #with open("pvals_divisor_modded.json","w") as file:
    #  file.write(json.dumps(list(items)))

    min_items = min_acc(items[::-1], "fdr_pval")[::-1]

    with open("pvals_divisor_min.json","w") as file:
      file.write(json.dumps(list(min_items)))

    for idx in range(0, len(min_items)):
        item = min_items[idx]

        if item['fdr_pval'] > 1.0:
            item['fdr_pval'] = 1.0

    return min_items
    #for idx in range(0, len(items)):
    #    item = items[idx]
    #    if item['fdr_pval'] > 1.0:
    #        item['fdr_pval'] = 1.0

    '''
    pvals_corrected = np.minimum.accumulate(pvals_corrected_raw[::-1])[::-1]
    pvals_corrected[pvals_corrected > 1.0] = 1.0
    pvals_corrected = pvals_corrected[sortrevind].reshape(shape_init)
    reject = reject[sortrevind].reshape(shape_init)
    '''

def fdr_correction(pvals, alpha=0.05, method='indep'):
    """P-value correction with False Discovery Rate (FDR).
    Correction for multiple comparison using FDR :footcite:`GenoveseEtAl2002`.
    This covers Benjamini/Hochberg for independent or positively correlated and
    Benjamini/Yekutieli for general or negatively correlated tests.
    Parameters
    ----------
    pvals : array_like
        Set of p-values of the individual tests.
    alpha : float
        Error rate.
    method : 'indep' | 'negcorr'
        If 'indep' it implements Benjamini/Hochberg for independent or if
        'negcorr' it corresponds to Benjamini/Yekutieli.
    Returns
    -------
    reject : array, bool
        True if a hypothesis is rejected, False if not.
    pval_corrected : array
        P-values adjusted for multiple hypothesis testing to limit FDR.
    References
    ----------
    .. footbibliography::
    """
    pvals = np.asarray(pvals)
    shape_init = pvals.shape
    pvals = pvals.ravel()

    pvals_sortind = np.argsort(pvals)
    pvals_sorted = pvals[pvals_sortind]
    sortrevind = pvals_sortind.argsort()

    if method in ['i', 'indep', 'p', 'poscorr']:
        ecdffactor = _ecdf(pvals_sorted)
    elif method in ['n', 'negcorr']:
        cm = np.sum(1. / np.arange(1, len(pvals_sorted) + 1))
        ecdffactor = _ecdf(pvals_sorted) / cm
    else:
        raise ValueError("Method should be 'indep' and 'negcorr'")

    with open("ecdf.json","w") as file:
      file.write(json.dumps(list(ecdffactor)))

    reject = pvals_sorted < (ecdffactor * alpha)
    if reject.any():
        rejectmax = max(np.nonzero(reject)[0])
    else:
        rejectmax = 0
    reject[:rejectmax] = True

    with open("pvals_sorted_fdr.json","w") as file:
      file.write(json.dumps(list(pvals_sorted)))

    pvals_corrected_raw = pvals_sorted / ecdffactor
    
    with open("pvals_corrected_raw_fdr.json","w") as file:
      file.write(json.dumps(list(pvals_corrected_raw)))
    
    pvals_corrected = np.minimum.accumulate(pvals_corrected_raw[::-1])[::-1]

    with open("pvals_corrected_min.json","w") as file:
      file.write(json.dumps(list(pvals_corrected)))
    
    pvals_corrected[pvals_corrected > 1.0] = 1.0
    with open("pvals_corrected_cap.json","w") as file:
      file.write(json.dumps(list(pvals_corrected)))
    
    pvals_corrected = pvals_corrected[sortrevind].reshape(shape_init)

    with open("pvals_corrected_shaped.json","w") as file:
      file.write(json.dumps(list(pvals_corrected)))

    reject = reject[sortrevind].reshape(shape_init)
    return reject, pvals_corrected


def bonferroni_correction(pval, alpha=0.05):
    """P-value correction with Bonferroni method.
    Parameters
    ----------
    pval : array_like
        Set of p-values of the individual tests.
    alpha : float
        Error rate.
    Returns
    -------
    reject : array, bool
        True if a hypothesis is rejected, False if not.
    pval_corrected : array
        P-values adjusted for multiple hypothesis testing to limit FDR.
    """
    pval = np.asarray(pval)
    pval_corrected = pval * float(pval.size)
    # p-values must not be larger than 1.
    pval_corrected = pval_corrected.clip(max=1.)
    reject = pval_corrected < alpha
    return reject, pval_corrected

def modded_bonferroni_correction(items, pval_key, alpha=0.05):
    """P-value correction with Bonferroni method.
    Parameters
    ----------
    pval : array_like
        Set of p-values of the individual tests.
    alpha : float
        Error rate.
    Returns
    -------
    reject : array, bool
        True if a hypothesis is rejected, False if not.
    pval_corrected : array
        P-values adjusted for multiple hypothesis testing to limit FDR.
    """
    arr_len = len(items)
    #pval = np.asarray(pval)
    for item in items:
        corrected_pval = item[pval_key] * float(arr_len)
        if corrected_pval > 1:
            corrected_pval = 1
        item['bonferroni_pval'] = corrected_pval

        item['bonferroni_reject'] = bool(corrected_pval < alpha)

def read_soft_file(filepath):
    with open(filepath) as f:
        lines = f.readlines()

        phenotypes = {
            'responder': [],
            'nonresponder': []
        }

        counter = 0
        for l in lines:
            line = l.strip()
            if line.startswith("!dataset_table_begin"):
                counter +=1
                break
            elif line.startswith('#GSM'):
                subject_metadata = line.split(";")
                phen = subject_metadata[4].strip()
                #if 'GSM125238' in line:
                #    print('seen')
                if phen != "not assessable":
                    responder_id = subject_metadata[0].split(" ")[0][1:]

                    phenotypes[phen].append(responder_id)
            counter += 1

        column_headers = lines[counter].strip().split("\t")
        #if 'GSM125238' in column_headers:
        #    print('seen in headers')
        counter +=1

        probes = []

        for i in range(counter, len(lines)):
            named_probe_data = {}
            probe_data = lines[i].strip().split("\t")
            for j in range(0, len(probe_data)):
                named_probe_data[column_headers[j]] = probe_data[j]
            probes.append(named_probe_data)

        #print(len(phenotypes['responder']), len(phenotypes['nonresponder']))
        return probes, phenotypes

def index(array, item):
    #print(array)
    for idx, val in np.ndenumerate(array):
        if val == item:
            return idx

def process_data(probes, phenotypes):

    from scipy import stats

    analyzed_probes = []

    counter = 0

    for probe in probes:
        #if counter == 3:
        #    return analyzed_probes

        responder_scores = []
        nonresponder_scores = []

        not_seen = []
        score_data = []

        for rid in phenotypes['responder']:
            if rid in probe:
                probe_val = probe[rid]
                if probe_val != 'null':
                    responder_scores.append(float(probe_val))
                    score_data.append({'group': 'responder', 'value': float(probe_val)})
            else:
                not_seen.append(rid)

        for nrid in phenotypes['nonresponder']:
            if nrid in probe:
                probe_val = probe[nrid]
                if probe_val != 'null':
                    nonresponder_scores.append(float(probe_val))
                    score_data.append({'group': 'nonresponder', 'value': float(probe_val)})
            else:
                not_seen.append(nrid)

        #probe_name = probe['ID_REF'].replace("/", "_")
        #pd.DataFrame(score_data).to_csv(f"data/{probe_name}.csv")
        statistic, p_value = stats.ttest_ind(responder_scores, nonresponder_scores, equal_var=True)

        #bon =  mne.stats.bonferroni_correction(p_value, alpha=0.01)
        #print(bon)
        #bon_reject, bon_pval_corrected = bon
        
        #fdr = mne.stats.fdr_correction(p_value, alpha=0.01, method='indep')
        #print(fdr)
        #fdr_reject, fdr_pval_corrected = fdr

        #if probe['ID_REF'] == '41113_at':
        #    print(f'4113_at pvalue: {p_value}')
        #    print(len(responder_scores))
        #    print(len(nonresponder_scores))
        #    print(not_seen)

        #print(probe)
        #print(f'corrected: ' + str(type(fdr_pval_corrected)))
        #print(fdr_pval_corrected)
        if 'IDENTIFIER' not in probe:
            pass
            #print(probe)
        else:
            analyzed_probes.append({
                'counter': counter,
                'probeid': probe['ID_REF'],
                'gene': probe['IDENTIFIER'],
                'pval': p_value
                #'bon_pval': bon_pval_corrected,
                #'bon_reject': bool(bon_reject),
                #'fdr_pval': index(fdr_pval_corrected, 0),
                #'fdr_reject': bool(fdr_reject)
            })
        counter += 1
    return analyzed_probes



probes, phenotypes = read_soft_file("GDS3116.soft")

extracted_data = {
    'probes': probes,
    'phenotypes': phenotypes
}

with open("extracted_data.json","w") as file:
  file.write(json.dumps(extracted_data))

with open("extracted_data.json", "r") as file:
    extracted_data = json.loads(file.read())

p_values = []
processed_probes = process_data(extracted_data['probes'], extracted_data['phenotypes'])
for processed_probe in processed_probes:
    p_values.append(processed_probe['pval'])

with open("p_values.json","w") as file:
  file.write(json.dumps(list(p_values)))

alpha = 0.01
bon =  mne.stats.bonferroni_correction(p_values, alpha=alpha)
#print(bon)
bon_reject, bon_pval_corrected = (list(bon[0]), list(bon[1]))
modded_bonferroni_correction(processed_probes, 'pval', alpha)

with open("bon_pval_correct.json","w") as file:
  file.write(json.dumps(list(bon_pval_corrected)))

print(processed_probes)
with open("bon_modded.json","w") as file:
  file.write(json.dumps(processed_probes))

#for idx in range(0, len(p_values)):
#    if bon_reject[idx]:
#        probeid = processed_probes[idx]['probeid']
#        gene = processed_probes[idx]['gene']
#        print(f'{probeid} {gene} {list(bon_pval_corrected)[idx]}')

#for bon in list(bon_reject):
#    if bon:
#        print(bon)

print()

fdr = fdr_correction(p_values, alpha=alpha, method='indep')

#fdr = mne.stats.fdr_correction(p_values, alpha=alpha, method='indep')
fdr_reject, fdr_pval_corrected = (list(fdr[0]), list(fdr[1]))

#for idx in range(0, len(p_values)):
#    if fdr_reject[idx]:
#        probeid = processed_probes[idx]['probeid']
#        gene = processed_probes[idx]['gene']
#        print(f'{probeid} {gene} {list(fdr_pval_corrected)[idx]}')


with open("fdr_pval_corrected.json","w") as file:
  file.write(json.dumps(list(fdr_pval_corrected)))

#modded_fdr_correction(items, pval_key, alpha=0.05, method='indep')
md = modded_fdr_correction(processed_probes, 'pval', alpha, 'indep')
#print(processed_probes)
with open("fdr_modded.json","w") as file:
  file.write(json.dumps(md))


