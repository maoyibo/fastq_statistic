from operator import index
from matplotlib import pyplot as plt
import subprocess
from typing import DefaultDict, List, Optional
from io import StringIO
from pathlib import Path
import typer
import pandas as pd
import numpy as np
import matplotlib as mpl
from matplotlib.axes import Axes
mpl.use('agg')


app = typer.Typer()


def run_statistic(filepath: Path) -> pd.DataFrame:
    script_folder = Path(__file__).absolute().parent
    binary_path = script_folder/"fastq_statistic"
    stat = subprocess.check_output(f"{binary_path} {filepath}", shell=True)
    df = pd.read_csv(StringIO(stat.decode()), header=0, index_col=0)
    return df


def report(df1: pd.DataFrame, df2: pd.DataFrame, sampleid: str):
    df = pd.DataFrame(0, index=[], columns=[])
    df.loc[0, "#SampleID"] = sampleid
    df.loc[0, "ReadPair"] = df1["Count"].sum()
    df.loc[0, "Reads"] = df1["Count"].sum()+df2["Count"].sum()
    df.loc[0, "Bases"] = df1[['A', 'G', 'C', 'T', 'N']
                             ].sum().sum() + df2[['A', 'G', 'C', 'T', 'N']].sum().sum()
    df.loc[0, "GCBases"] = df1[['G', 'C']].sum().sum() + \
        df2[['G', 'C']].sum().sum()
    df.loc[0, "NBases"] = df1['N'].sum() + df2['N'].sum()
    df.loc[0, "Q20Bases"] = df1[[x for x in df1.columns[6:] if int(
        x)-33 >= 20]].sum().sum() + df2[[x for x in df1.columns[6:] if int(x)-33 >= 20]].sum().sum()
    df.loc[0, "Q30Bases"] = df1[[x for x in df1.columns[6:] if int(
        x)-33 >= 30]].sum().sum() + df2[[x for x in df1.columns[6:] if int(x)-33 >= 30]].sum().sum()
    df.loc[0, "GC"] = df.loc[0, "GCBases"]/df.loc[0, "Bases"]
    df.loc[0, "N"] = df.loc[0, "NBases"]/df.loc[0, "Bases"]
    df.loc[0, "Q20"] = df.loc[0, "Q20Bases"]/df.loc[0, "Bases"]
    df.loc[0, "Q30"] = df.loc[0, "Q30Bases"]/df.loc[0, "Bases"]
    read1_length = []
    for (cycle, count) in df1.loc[df1["Count"] != 0, ["Count"]].iterrows():
        read1_length.extend([cycle]*count["Count"])
    read2_length = []
    for (cycle, count) in df2.loc[df2["Count"] != 0, ["Count"]].iterrows():
        read2_length.extend([cycle]*count["Count"])
    df.loc[0, "MinLength"] = np.min(read1_length+read2_length)
    df.loc[0, "MeanLength"] = np.mean(read1_length+read2_length)
    df.loc[0, "MedianLength"] = np.median(read1_length+read2_length)
    df.loc[0, "MaxLength"] = np.max(read1_length+read2_length)
    # Read1
    df.loc[0, "Read1Bases"] = df1[['A', 'G', 'C', 'T', 'N']].sum().sum()
    df.loc[0, "Read1GCBases"] = df1[['G', 'C']].sum().sum()
    df.loc[0, "Read1NBases"] = df1['N'].sum()
    df.loc[0, "Read1Q20Bases"] = df1[[
        x for x in df1.columns[6:] if int(x)-33 >= 20]].sum().sum()
    df.loc[0, "Read1Q30Bases"] = df1[[
        x for x in df1.columns[6:] if int(x)-33 >= 30]].sum().sum()
    df.loc[0, "Read1GC"] = df.loc[0, "Read1GCBases"]/df.loc[0, "Read1Bases"]
    df.loc[0, "Read1N"] = df.loc[0, "Read1NBases"]/df.loc[0, "Read1Bases"]
    df.loc[0, "Read1Q20"] = df.loc[0, "Read1Q20Bases"]/df.loc[0, "Read1Bases"]
    df.loc[0, "Read1Q30"] = df.loc[0, "Read1Q30Bases"]/df.loc[0, "Read1Bases"]
    df.loc[0, "Read1MinLength"] = np.min(read1_length)
    df.loc[0, "Read1MeanLength"] = np.mean(read1_length)
    df.loc[0, "Read1MedianLength"] = np.median(read1_length)
    df.loc[0, "Read1MaxLength"] = np.max(read1_length)
    # Read2
    if df2["Count"].sum() != 0:
        df.loc[0, "Read2Bases"] = df2[['A', 'G', 'C', 'T', 'N']].sum().sum()
        df.loc[0, "Read2GCBases"] = df2[['G', 'C']].sum().sum()
        df.loc[0, "Read2NBases"] = df2['N'].sum()
        df.loc[0, "Read2Q20Bases"] = df2[[
            x for x in df2.columns[6:] if int(x)-33 >= 20]].sum().sum()
        df.loc[0, "Read2Q30Bases"] = df2[[
            x for x in df2.columns[6:] if int(x)-33 >= 30]].sum().sum()
        df.loc[0, "Read2GC"] = df.loc[0, "Read2GCBases"] / \
            df.loc[0, "Read2Bases"] if df.loc[0, "Read2Bases"] != 0 else 0
        df.loc[0, "Read2N"] = df.loc[0, "Read2NBases"] / \
            df.loc[0, "Read2Bases"] if df.loc[0, "Read2Bases"] != 0 else 0
        df.loc[0, "Read2Q20"] = df.loc[0, "Read2Q20Bases"] / \
            df.loc[0, "Read2Bases"] if df.loc[0, "Read2Bases"] != 0 else 0
        df.loc[0, "Read2Q30"] = df.loc[0, "Read2Q30Bases"] / \
            df.loc[0, "Read2Bases"] if df.loc[0, "Read2Bases"] != 0 else 0
        df.loc[0, "Read2MinLength"] = np.min(
            read2_length) if len(read2_length) != 0 else 0
        df.loc[0, "Read2MeanLength"] = np.mean(
            read2_length) if len(read2_length) != 0 else 0
        df.loc[0, "Read2MedianLength"] = np.median(
            read2_length) if len(read2_length) != 0 else 0
        df.loc[0, "Read2MaxLength"] = np.max(
            read2_length) if len(read2_length) != 0 else 0
    return df


def collect_plot_data(df: pd.DataFrame):
    x = []
    total_reads = df["Count"].sum()
    data = DefaultDict(list)
    q20_column = [x for x in df.columns[6:] if int(x)-33 >= 20]
    q30_column = [x for x in df.columns[6:] if int(x)-33 >= 30]
    for cycle, row in df.iterrows():
        x.append(cycle)
        reads_count = np.sum(row[['A', 'G', 'C', 'T', 'N']])
        data['Count'].append(row['Count']/total_reads)
        for n in ['A', 'G', 'C', 'T', 'N']:
            data[n].append(row[n]/reads_count)
        data["Q20"].append(np.sum(row[q20_column])/reads_count)
        data["Q30"].append(np.sum(row[q30_column])/reads_count)
    return x, data


def plot_df(ax: Axes, x: List[int], data: pd.DataFrame, label="Read1"):
    ax.bar(x, data["Count"], alpha=0.7,
           align='center', label=f"{label} Read Count")
    for i, l in zip(['A', 'G', 'C', 'T', ], ['-', '--', '-.', ':']):
        ax.plot(x, data[i], linestyle=l, label=f"{label} {i}")
    ax.scatter(x, [x*100 for x in data['N']], s=1.2, color='gray',
               alpha=0.7, linewidths=0.5, label=f"{label} N(%%)")
    ax.plot(x, data['Q20'], linestyle=(
        0, (3, 5, 1, 5, 1, 5)), label=f"{label} Q20")
    ax.plot(x, data['Q30'], linestyle=(
        0, (3, 1, 1, 1, 1, 1)), label=f"{label} Q30")
    ax.set_xlabel(f"{label} Cycles")
    ax.set_ylabel("Percentage")
    ax.set_yticks([x/100. for x in range(0, 101, 20)])
    ax.set_yticklabels([f'{x}%' for x in range(0, 101, 20)])


def plot(df1: pd.DataFrame, df2: pd.DataFrame, filepth: Path):
    if df2["Count"].sum() != 0:
        fig, (ax1, ax2) = plt.subplots(ncols=2)
        df1_x, df1_data = collect_plot_data(df1)
        df2_x, df2_data = collect_plot_data(df2)
        plot_df(ax1, df1_x, df1_data)
        plot_df(ax2, df2_x, df2_data, "Read2")
        ax1.legend(bbox_to_anchor=(0., 1.02, 1., .102),
                   loc='lower left', ncol=4, mode="expand", borderaxespad=0.)
        ax2.legend(bbox_to_anchor=(0., 1.02, 1., .102),
                   loc='lower left', ncol=4, mode="expand", borderaxespad=0.)
        fig.set_figwidth(fig.get_figwidth()*2.3)
    else:
        fig, ax1 = plt.subplots()
        df1_x, df1_data = collect_plot_data(df1)
        plot_df(ax1, df1_x, df1_data)
        ax1.legend(bbox_to_anchor=(0., 1.02, 1., .102),
                   loc='lower left', ncol=4, mode="expand", borderaxespad=0.)
        fig.set_figwidth(fig.get_figwidth()*1.2)
    fig.savefig(filepth, dpi=300, bbox_inches='tight')


@app.command()
def main(read1: Path = typer.Argument(..., help="Read1 fastq path or fastq path"),
         read2: Optional[Path] = typer.Argument(
             None, help="Read2 filepath or None"),
         sampleid: Optional[str] = typer.Option(
             None, help="SampleID, default is the first item of filename splited by underscore(_)"),
         result: Optional[Path] = typer.Option(None, help="Result csv file name, plot with use the same name.")):
    df1 = run_statistic(read1)
    if read2:
        df2 = run_statistic(read2)
    else:
        df2 = pd.DataFrame(0, index=df1.index, columns=df1.columns)
    # sampleid
    if not sampleid:
        sampleid = read1.stem.split('_')[0]
    df = report(df1, df2, sampleid)
    if result:
        df.to_csv(result, index=False)
        plot(df1, df2, result.with_suffix(".png"))
    else:
        result = Path(f'{sampleid}_statistic.csv')
        df.to_csv(result, index=False)
        plot(df1, df2, result.with_suffix(".png"))
