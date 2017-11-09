import json
from os.path import join
import sys

import numpy as np


def main(dirname):
    result = {}
    # do coordinates
    with open(join(dirname, 'coordinates.json')) as handle:
        result['coordinates'] = list(int(x) for x in json.load(handle)['coordinates'])

    with open(join(dirname, 'dates.json')) as handle:
        result['dates'] = list({'date': k, 'name': v} for k, v in json.load(handle).items())

    with open(join(dirname, 'trees.json')) as handle:
        result['tree'] = json.load(handle)['tree']

    with open(join(dirname, 'sequences.json')) as handle:
        sequences = json.load(handle)

    result['mrca'] = sequences['MRCA']
    result['reference'] = sequences['Reference']

    result['ancestors'] = list({'name': k, 'sequence': v} for k, v in sequences['Ancestors'].items())
    
    with open(join(dirname, 'copynumbers.json')) as handle:
        copynumbers = json.load(handle)

    with open(join(dirname, 'manifold.json')) as handle:
        embedding = json.load(handle)

    result['sequences'] = list({
        'name': k,
        'date': date,
        'sequence': v,
        'copynumber': copynumbers[k],
        'x': embedding[k][0],
        'y': embedding[k][1]
    } for date, seqs in sequences['Observed'].items()
                               for k, v in seqs.items())

    protein_metrics = {
        'single': [],
        'paired': [],
    }

    # combine all protein metrics
    with open(join(dirname, 'js_divergence.json')) as handle:
        js_divergence = json.load(handle)

    combined_js = list(np.vstack(list(np.array(v) for v in js_divergence.values())).max(axis=0))
    js_divergence['Combined'] = combined_js
    protein_metrics['single'].append({
        'name': "JS divergence",
        'paired': False,
        'data': list({
            'date': k,
            'values': v
        } for k, v in js_divergence.items())
    })

    with open(join(dirname, 'rates.json')) as handle:
        rates = json.load(handle)

    protein_metrics['single'].append({
        'name': "Entropy",
        'paired': False,
        'data': list({
            'date': k,
            'values': list(row[4] for row in v)
        } for k, v in rates.items())
    })

    dnds = {
        'name': 'dNdS',
        'paired': True,
        'data': []
    }
    dnds['data'].append({
        'name': "Mean dS",
        'data': list({
            'date': k,
            'values': list(row[1] for row in v)
        } for k, v in rates.items())
    })
    dnds['data'].append({
        'name': "Mean dN",
        'data': list({
            'date': k,
            'values': list(row[0] for row in v)
        } for k, v in rates.items())
    })
    protein_metrics['paired'].append(dnds)
    result['protein_metrics'] = protein_metrics

    # combine all region metrics
    with open(join(dirname, 'rates_pheno.json')) as handle:
        rates_pheno = json.load(handle)

    regions = []
    for region in rates_pheno:
        regions.append({
            'name': region['Segment'],
            'date': region['Date'],
            'evo_metrics': [
                {
                    'name': 'S diversity',
                    'value': float(region['s_diversity'])
                },
                {
                    'name': 'nS diversity',
                    'value': float(region['ns_diversity'])
                },
                {
                    'name': 'Total diversity',
                    'value': float(region['total_diversity'])
                },
                {
                    'name': 'dS diversity',
                    'value': float(region['ds_diversity'])
                },
                {
                    'name': 'dN diversity',
                    'value': float(region['dn_diversity'])
                },
                {
                    'name': 'S divergence',
                    'value': float(region['s_divergence'])
                },
                {
                    'name': 'nS divergence',
                    'value': float(region['ns_divergence'])
                },
                {
                    'name': 'Total divergence',
                    'value': float(region['total_divergence'])
                },
                {
                    'name': 'dS divergence',
                    'value': float(region['ds_divergence'])
                },
                {
                    'name': 'dN divergence',
                    'value': float(region['dn_divergence'])
                }
            ],
            'pheno_metrics': [
                {
                    'name': 'Length',
                    'value': float(region['Length'])
                },
                {
                    'name': 'PNGS',
                    'value': float(region['PNGS'])
                },
                {
                    'name': 'Isoelectric Point',
                    'value': float(region['IsoelectricPoint'])
                }
            ]
        })
    result['region_metrics'] = regions
    with open(join(dirname, 'session.json'), 'w') as handle:
        json.dump(result, handle, sort_keys=True, indent=2)
        
if __name__ == '__main__':
    dirname = sys.argv[1]
    main(dirname)
