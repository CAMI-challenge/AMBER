#!/usr/bin/env python

import numpy as np


def filter_tail(data, percentage):
    value = float(percentage) / 100.0
    # sort bins by increasing predicted size
    data = sorted(data, key=lambda x: x['predicted_size'])
    sum_of_predicted_sizes = sum([int(float(bin['predicted_size'])) for bin in data])
    cumsum_of_predicted_sizes = 0
    for bin in data:
        predicted_size = int(float(bin['predicted_size']))
        cumsum_of_predicted_sizes += predicted_size
        if float(cumsum_of_predicted_sizes) / float(sum_of_predicted_sizes) <= value:
            bin['precision'] = np.nan
        else:
            break
    return data
