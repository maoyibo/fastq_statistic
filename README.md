Fastq Statistic
===============
Calculate statistics for Fastq Files.  

# Requirement
Python: 3.8 or upper  

# Install
Download the release whl file.  

```bash
user@linux:~$ python3 -m pip install fastq-statistic
```

# Usage
```bash
user@linux:~$ fastq-stat --help
Usage: fastq-stat [OPTIONS] READ1 [READ2]

Arguments:
  READ1    Read1 fastq path or fastq path  [required]
  [READ2]  Read2 filepath or None

Options:
  --sampleid TEXT                 SampleID, default is the first item of
                                  filename splited by underscore(_)

  --result PATH                   Result csv file name, plot with use the same
                                  name.

  --reserve-data / --no-reserve-data
                                  Reserve fastq statistic intermediate data.
                                  [default: False]

  --plot / --no-plot              Plot fastq statistic data.  [default: True]
  --install-completion [bash|zsh|fish|powershell|pwsh]
                                  Install completion for the specified shell.
  --show-completion [bash|zsh|fish|powershell|pwsh]
                                  Show completion for the specified shell, to
                                  copy it or customize the installation.

  --help                          Show this message and exit.

```