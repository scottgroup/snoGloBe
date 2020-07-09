#!/usr/bin/env python3

import pandas as pd
import numpy as np
import subprocess
from io import StringIO
import os


def merge_consecutive_windows(df_bed, step, nb_min, last_pos, agg_fct=['mean', 'min', 'max', 'count']):
    d_agg = {
        'mean': np.nanmean,
        'max': np.nanmax,
        'min': np.nanmin,
        'median': np.nanmedian,
        'count': 'count'
    }
    np_agg_fct = [d_agg[f] for f in agg_fct]
    wlength = 13
    # snoRNA and target windows are always in 5' -> 3' orientation
    # consecutive windows are in opposite numerical direction, if twindow increases, sno_window should decrease

    df_bed.loc[(df_bed.count1.isnull()), 'next_snowindow'] = df_bed.sno_window + step
    df_bed.loc[(df_bed.count1.isnull()), 'cons_twindow'] = df_bed.wstart - step
    df_merged = df_bed[['target_id', 'sno_id', 'sno_window', 'wstart', 'score', 'sno_window_start', 'sno_window_end',
                        'target_window_start', 'target_window_end',
                        'min_score1', 'max_score1', 'mean_score1', 'count1']].copy(deep=True)
    df_merged['border_snowindow'] = df_merged.sno_window
    df_merged['border_twindow'] = df_merged.wstart
    cols = ['target_id', 'sno_id', 'score', 'next_snowindow', 'cons_twindow']
    df_tomerge = df_bed[cols].copy(deep=True)
    left_cols = ['target_id', 'sno_id', 'sno_window', 'wstart']
    right_cols = ['target_id', 'sno_id', 'next_snowindow', 'cons_twindow']
    original_windows = df_bed[left_cols].values.tolist()
    cons_windows = df_bed[right_cols].values.tolist()
    i = 1
    while len(set(map(tuple, original_windows)) & set(map(tuple, cons_windows))) > 0:
        df_merged = df_merged.merge(df_tomerge,
                                    left_on=left_cols,
                                    right_on=right_cols,
                                    suffixes=['', '_%d' % i],
                                    how='left')

        df_merged.loc[~df_merged['score_%d' % i].isnull(), 'border_snowindow'] -= step
        df_merged.loc[~df_merged['score_%d' % i].isnull(), 'border_twindow'] += step
        df_merged.loc[~df_merged['score_%d' % i].isnull(), 'next_snowindow'] += step
        df_merged.loc[~df_merged['score_%d' % i].isnull(), 'cons_twindow'] -= step

        if i < nb_min:
            df_merged = df_merged[~(df_merged['score_%d' % i].isnull()) | (df_merged['border_twindow'] == last_pos)]
        original_windows = df_merged[left_cols].values.tolist()
        cons_windows = df_merged[right_cols].values.tolist()
        cols_i = [j + '_%d' % i if j.startswith('score') else j for j in cols]
        df_tomerge = df_merged[cols_i].copy(deep=True)
        df_tomerge = df_tomerge.rename(columns={'score_%d' % i: 'score'})
        cols_to_keep = ['target_id', 'sno_id', 'sno_window', 'wstart', 'border_snowindow', 'border_twindow',
                        'sno_window_start',
                        'sno_window_end', 'target_window_start', 'target_window_end',
                        'min_score1', 'max_score1', 'mean_score1', 'count1']
        cols_to_keep.extend([j for j in df_merged.columns if j.startswith('score')])
        df_merged = df_merged[cols_to_keep].copy(deep=True)
        df_merged = df_merged.sort_values(['score_%d' % i])
        df_merged = df_merged.drop_duplicates(subset=['target_id', 'sno_id', 'border_snowindow', 'border_twindow'],
                                              keep='first')
        i += 1
    agg_cols = [f + '_score' if f != 'count' else f for f in agg_fct]
    df_merged[agg_cols] = df_merged[[j for j in df_merged.columns if j.startswith('score')]].agg(np_agg_fct, axis=1)

    df_merged['border_snowindow'] = np.nanmin(df_merged[['border_snowindow', 'sno_window_start']], axis=1)
    df_merged['min_score'] = np.nanmin(df_merged[['min_score', 'min_score1']], axis=1)
    df_merged['max_score'] = np.nanmax(df_merged[['max_score', 'max_score1']], axis=1)
    df_merged.loc[(~df_merged.mean_score1.isnull())
                  & (~df_merged.mean_score.isnull()), 'mean_score'] = (df_merged['count1'] * df_merged['mean_score1']
                                                                       + df_merged['count'] * df_merged['mean_score']) \
                                                                      / (df_merged['count1'] + df_merged['count'])
    df_merged.loc[(~df_merged.mean_score1.isnull())
                  & (df_merged.mean_score.isnull()), 'mean_score'] = df_merged['mean_score1']
    df_merged.loc[~df_merged.count1.isnull(), 'count'] += df_merged.count1
    df_merged['sno_window_end'] -= wlength
    df_merged['sno_window'] = np.nanmax(df_merged[['sno_window', 'sno_window_end']], axis=1)
    df_merged['wstart'] = np.nanmin(df_merged[['wstart', 'target_window_start']], axis=1)
    final_cols = ['target_id', 'sno_id', 'border_snowindow', 'sno_window', 'wstart', 'border_twindow']
    final_cols.extend(agg_cols)
    df_merged = df_merged[final_cols]

    df_merged = df_merged.rename(
        columns={
            'border_snowindow': 'sno_window_start',
            'sno_window': 'sno_window_end',
            'wstart': 'target_window_start',
            'border_twindow': 'target_window_end'
        }
    )
    df_merged['target_window_end'] += wlength
    df_merged['sno_window_end'] += wlength
    df_merged[agg_cols] = df_merged[agg_cols].astype(float).round(3)
    df_merged['count'] = df_merged['count'].astype(int)
    df_prev = df_merged[df_merged.target_window_end - 13 == last_pos].copy(deep=True)
    df_merged = df_merged[df_merged.target_window_end - 13 != last_pos]
    df_merged = df_merged[df_merged['count'] >= nb_min]
    return df_merged, df_prev


def cons_windows(infile, nb_min, chunksize, step, nb_threads):
    # File should be sorted by target start to be compatible with chunks,
    # otherwise, some consecutive windows could not be merged
    wlength = 13
    # sort file by target_id, sno_id, target_start
    sort_cmd = ['sort --parallel=%d -t, -k1,1 -k3,3 -k4,4n %s' % (nb_threads, infile)]
    result = subprocess.Popen(sort_cmd, stdout=subprocess.PIPE, shell=True)
    sorted_input = StringIO(result.communicate()[0].decode('utf-8'))
    # merge windows
    temp_output = infile + '.sorted.merged.csv'
    df_prev = pd.DataFrame(columns=['target_id', 'wstart', 'end', 'sno_id', 'sno_window',
                                    'sno_window_start', 'sno_window_end', 'target_window_start', 'target_window_end',
                                    'min_score1', 'max_score1', 'mean_score1', 'count1'])
    for i, df in enumerate(pd.read_csv(sorted_input,
                                       names=['sno_id', 'sno_window', 'target_id', 'wstart', 'score'],
                                       chunksize=chunksize)):
        if df.empty:
            continue
        last_pos = df.wstart.values[-1]
        df = pd.concat([df, df_prev], sort=False, ignore_index=True)

        df_merged, df_prev = merge_consecutive_windows(df, step, nb_min, last_pos)
        df_prev = df_prev.rename(
            columns={
                'min_score': 'min_score1',
                'max_score': 'max_score1',
                'mean_score': 'mean_score1',
                'count': 'count1'
            })

        df_prev['sno_window'] = df_prev.sno_window_start
        df_prev['wstart'] = df_prev.target_window_end - wlength
        df_prev['wend'] = df_prev.target_window_end
        if i == 0:
            mode = 'w'
            header = True
        else:
            mode = 'a'
            header = False
        df_merged[['target_id', 'target_window_start', 'target_window_end', 'sno_id',
                   'sno_window_start', 'sno_window_end', 'count',
                   'mean_score', 'min_score', 'max_score']].to_csv(temp_output,
                                                                   mode=mode,
                                                                   header=header,
                                                                   index=False)
    if df_prev.empty is False:
        df_prev[df_prev['count1'] >= nb_min][['target_id', 'target_window_start', 'target_window_end',
                                              'sno_id', 'sno_window_start', 'sno_window_end', 'count1',
                                              'mean_score1', 'min_score1', 'max_score1']].to_csv(temp_output,
                                                                                                 mode='a',
                                                                                                 header=False,
                                                                                                 index=False)
    os.remove(infile)
    os.rename(temp_output, infile)
