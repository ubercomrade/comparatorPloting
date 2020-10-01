import argparse
import sys
import csv
import os
import matplotlib.pyplot as plt
from bisect import bisect_left, bisect_right, bisect
import math
import numpy as np
import itertools
from matplotlib import rcParams
from operator import itemgetter
from matplotlib_venn import venn3, venn3_circles, venn2, venn2_circles

def read_list(path):
    container = list()
    append = container.append
    with open(path, "r") as file:
        for line in file:
            append(float(line.strip()))
    file.close()
    return(container)


def read_score_fpr_table(path):
    container = {'score': [],
                'fpr': []}
    with open(path) as file:
        for line in file:
            score, fpr = [float(i) for i in line.strip().split()]
            container['score'].append(score)
            container['fpr'].append(fpr)
    file.close()
    container['score'] = container['score'][::-1]
    container['fpr'] = container['fpr'][::-1]
    return(container)


def from_score_to_fpr(scores, table):
    container = []
    append = container.append
    for s in scores:
        try:
            append(table['fpr'][:bisect_right(table['score'], s)][-1])
        except:
            append(0.01)
    return(container)


def read_sitega_log_fpr(path):
    container = []
    with open(path) as file:
        for line in file:
            fpr = float(line.strip().split()[1])
            container.append(fpr)
    return(container)


def creat_log_fpr(wd, tools):
    d = dict()
    for t in tools:
        scores = read_list(wd + '/scan-best/{}.scores.txt'.format(t.lower()))
        thresholds = read_score_fpr_table(wd + '/models/thresholds/{}_model_thresholds.txt'.format(t.lower()))
        fpr = from_score_to_fpr(scores, thresholds)
        fpr = [i if i != 0 else 1*10**(-99) for i in fpr]
        fpr_log = [-math.log10(j) for j in fpr]
        d[t] = fpr_log
    return(d)


def read_counts(path):
    with open(path) as csvfile:
        csvreader = csv.DictReader(csvfile, delimiter='\t')
        for row in csvreader:
            out = row
    for k in out.keys():
        out[k] = int(out[k])
    return(out)


def tool_to_color(tools):
    colors_tools = dict()
    colors = ['#FF0018', '#FBD704', '#0000F9', '#10D8B8']
    for index, t in enumerate(tools):
        colors_tools[t] = colors[index]
    return(colors_tools)


def read_bed(path):
    container = []
    with open(path) as csvfile:
        csvreader = csv.reader(csvfile, delimiter='\t')
        for row in csvreader:
            container.append({
                'chr': row[0], 'start': int(row[1]), 'end': int(row[2])})
    container = sorted(container, key = itemgetter('chr', 'start'))
    return(container)


def read_scan_file(path, scan_id):
    container = []
    with open(path) as csvfile:
        csvreader = csv.reader(csvfile, delimiter='\t')
        for row in csvreader:
            container.append({
                'chromosome': row[0], 'start': int(row[1]), 'end': int(row[2]),
                'name': int(row[3].split('_')[1]), 'score': float(row[4]),
                'strand': row[5], 'site': row[6], 'type': scan_id
            })
    return(container)

    
def is_intersect(interval, intervals):
    for i in intervals:
        if interval['start'] < i['end'] and interval['end'] > i['start']:
            return(1)
    return(0)
  
    
def get_indexes(peaks, sites):
    container = []
    append = container.append
    chrs = list(set([i['chr'] for i in peaks]))
    chrs.sort()
    index = 0
    for chr_ in chrs:
        sub_peaks = [i for i in peaks if i['chr'] == chr_]
        sub_sites = [i for i in sites if i['chr'] == chr_]
        if len(sub_sites) == 0:
            index += len(sub_peaks)
            continue      
        for p in sub_peaks:
            if is_intersect(p, sub_sites):
                append(index)
            index += 1               
    return(container)


def write_table(path, data):
    with open(path, 'w') as csvfile:
        fieldnames = list(data.keys())
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        writer.writerow(data)
    pass

        
def creat_petal_3(first_names, second_names, third_names):
    petal_labels = {'001': len(third_names - (second_names | first_names)),
     '010': len(second_names - (third_names | first_names)),
     '011': len((second_names & third_names) - (first_names)),
     '100': len(first_names - (second_names | third_names)),
     '101': len((first_names & third_names) - second_names),
     '110': len((first_names & second_names) - third_names),
     '111': len(first_names & second_names & third_names)}
    return(petal_labels)


def creat_petal_2(first_names, second_names):
    petal_labels = {'01': len(second_names - first_names),
     '10': len(first_names - second_names),
     '11': len(second_names & first_names)}
    return(petal_labels)


def heatmap_ploting(wd, tools, number_of_tools, path):
    d = creat_log_fpr(wd, tools)
    matrix = []
    coordinates = list(itertools.product(range(1, number_of_tools + 1), repeat=2))
    acc = 0

    for i in range(1, number_of_tools + 1):
        matrix.append(coordinates[acc:acc + number_of_tools])
        acc += number_of_tools

    coordinates = []
    for index, i in enumerate(matrix):
        coordinates += (i[:len(i) - index - 1])

    ftools = tools
    rtools = ftools[::-1]

    fig, axs = plt.subplots(squeeze=False, figsize=(14, 12.5))
    c = []
    for i, j in coordinates:
        tool1, tool2 = ftools[i - 1], rtools[j - 1]
        ax = plt.subplot2grid((len(tools) - 1, len(tools) - 1), (j - 1, i - 1))
        pcm = ax.hist2d(d[tool1], d[tool2], bins=(30, 30), density=True, vmin=0, vmax=1,
                        cmap=plt.cm.jet, range=np.array([(2.5, 6), (2.5, 6)]))
        ax.set_xticks([3,4,5,6])
        ax.set_yticks([3,4,5,6])
        ax.set_xlim(2.5, 6)
        ax.set_ylim(2.5, 6)
        if i - 1 == 0 and j - 1 == 0:
            ax.set_title(tool1.upper())
            ax.set_ylabel(tool2, size=30)
        elif i - 1 == 0:
            ax.set_ylabel(tool2, size=30)
        elif j - 1 == 0:
            ax.set_title(tool1)
        c.append(pcm)
    fig.subplots_adjust(right=0.82)
    #fig.colorbar(pcm[3], cax=cbar_ax, orientation='horizontal')
    fig.savefig(path, dpi=400, bbox_inches='tight', pad_inches = 0)
    return(0)


def pie_plots(wd, tools, number_of_tools, path):
    matrix = []
    coordinates = list(itertools.product(range(1, number_of_tools + 1), repeat=2))
    acc = 0

    for i in range(1, number_of_tools + 1):
        matrix.append(coordinates[acc:acc + number_of_tools])
        acc += number_of_tools

    coordinates = []
    for index, i in enumerate(matrix):
        coordinates += (i[:len(i) - index - 1])

    ftools = tools
    rtools = ftools[::-1]

    colors_back =['#4f4f4f', '#C6C6C6', '#FFFFFF']
    colors_tools = tool_to_color(tools)

    fig, axs = plt.subplots(squeeze=False, figsize=(5.5, 5.5), dpi=400)
    container = []
    for i, j in coordinates:
        tool1, tool2 = ftools[i - 1], rtools[j - 1]
        data = read_counts(wd + '/results/compare_{0}.{1}_counts.tsv'.format(tool1.lower(), tool2.lower()))
        v = list(data.values())
        vals = np.array([v[0], v[1], v[2], v[3], v[4]])
        l = list(data.keys())
        l = [a.split(':')[0] for a in l]
        labels = np.array([l[0], l[1], l[2], l[3], l[4]])
        ax = plt.subplot2grid((len(tools) - 1, len(tools) - 1), (j - 1, i - 1))
        pie = ax.pie(vals, autopct='%1.1f%%', labels=labels,
                    labeldistance=None, pctdistance=0.9,
                    colors= [colors_tools[tool1], colors_tools[tool2]]  + colors_back,
                    wedgeprops = {'linewidth': 1, 'edgecolor': 'black'})
        if i - 1 == 0 and j - 1 == 0:
            ax.set_title(tool1.upper())
            ax.set_ylabel(tool2.upper(), size=20)
        elif i - 1 == 0:
            ax.set_ylabel(tool2.upper(), size=20)
        elif j - 1 == 0:
            ax.set_title(tool1.upper())
        container.append(ax)

    # LEGEND #    
    handles, labels = list(), list()
    for i in container:
        out = i.get_legend_handles_labels()
        handles += out[0]
        labels += out[1]

    common_legend = dict()
    for h, l in zip(handles, labels):
        if l not in common_legend:
            common_legend[l] = h
    labels = [i.lower() for i in tools] + ['overlapped', 'not_overlapped', 'no_sites']
    handles = [common_legend[i] for i in labels]
    fig.legend(handles, labels, loc='lower right', fontsize=14)
    fig.savefig(path,
                   dpi=400, format='pdf', bbox_inches='tight', pad_inches = 0)
    return(0)


def venn_ploting(wd, tools, fpr, path):
    if os.path.isfile(wd + '/bed/peaks.bed'):
        peaks = read_bed(wd + '/bed/peaks.bed')
    else:
        peaks = read_bed(wd + '/bed/test_sample.bed')

    tools_peaks = dict()
    for t in tools:
        bed = read_bed(wd + '/scan/{0}_{1:.2e}.bed'.format(t.lower(), fpr) )
        tools_peaks[t] = set(get_indexes(peaks, bed))

    if len(tools) == 3:
        id_to_name_3 = {'001': '{0}'.format(tools[2]),
         '010': '{0}'.format(tools[1]),
         '011': '{0}&{1}'.format(tools[1], tools[2]),
         '100': '{0}'.format(tools[0]),
         '101': '{0}&{1}'.format(tools[0], tools[2]),
         '110': '{0}&{1}'.format(tools[0], tools[1]),
         '111': '{0}&{1}&{2}'.format(tools[0], tools[1], tools[2])}

        petal_labels = creat_petal_3(tools_peaks[tools[0]], tools_peaks[tools[1]], tools_peaks[tools[2]])
        data = dict()
        for k in petal_labels.keys():
            data[id_to_name_3[k]] = petal_labels[k]
        venn3(subsets=petal_labels, set_labels=tools, alpha=0.5, set_colors=['#FF0018', '#FBD704', '#0000F9'])
        venn3_circles(petal_labels, linewidth=.6)
    elif len(tools) == 2:
        id_to_name_2 = {'01': '{0}'.format(tools[1]),
         '10': '{0}'.format(tools[0]),
         '11': '{0}&{1}'.format(tools[0], tools[1])}

        petal_labels = creat_petal_2(tools_peaks[tools[0]], tools_peaks[tools[1]])
        data = dict()
        for k in petal_labels.keys():
            data[id_to_name_2[k]] = petal_labels[k]
        venn2(subsets=petal_labels, set_labels=tools, alpha=0.5, set_colors=['#FF0018', '#FBD704', '#0000F9'])
        venn2_circles(petal_labels, linewidth=.6)
    plt.savefig(path,
               format="pdf", dpi=400, bbox_inches='tight', pad_inches = 0)
    return(0)


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('directory', action='store', help='directory with data')
    parser.add_argument('models', action='store', metavar='N', nargs='+',
         help='list of models names used in analisys')
    parser.add_argument('-f', '--FPR', action='store', type=float, dest='fpr',
                        required=False, default=1.9*10**(-4), help='FPR, def=1.9*10^(-4)')
    if len(sys.argv) == 1:
        parser.print_help(sys.stderr)
        sys.exit(1)
    return(parser.parse_args())


def main():
    args = parse_args()
    wd = args.directory
    tools = args.models
    fpr = args.fpr
    number_of_tools = len(tools)
    plots_path = wd + '/plots'
    if not os.path.exists(plots_path):
        os.makedirs(plots_path)

    rcParams['font.family'] = 'arial'
    rcParams['font.monospace'] = ['Sans serif']
    rcParams['font.size'] = 18

    venn_ploting(wd, tools, fpr, plots_path + '/venn.pdf')
    heatmap_ploting(wd, tools, number_of_tools, plots_path + '/heatmap.pdf')
    pie_plots(wd, tools, number_of_tools, plots_path + '/pie.pdf')
    return(0)


if __name__ == "__main__":
    main()
    