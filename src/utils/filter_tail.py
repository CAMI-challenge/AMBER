#!/usr/bin/env python

import numpy as np


def filter_tail(data, percentage):
    value = percentage / 100.0
    # sort bins by increasing predicted size
    data = sorted(data, key=lambda x: x['predicted_size'])
    sum_of_predicted_sizes = sum([int(float(bin['predicted_size'])) for bin in data])
    cumsum_of_predicted_sizes = 0
    for bin in data:
        predicted_size = int(float(bin['predicted_size']))
        cumsum_of_predicted_sizes += predicted_size
        if cumsum_of_predicted_sizes / sum_of_predicted_sizes <= value:
            bin['purity_bp'] = np.nan
            bin['purity_seq'] = np.nan
        else:
            break
    return data
