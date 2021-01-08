Fastq Statistic
===============
Calculate statistics for Fastq Files.  

# Requirement
Python: 3.6 or upper  

# Install
Download the release whl file.  

```bash
user@linux:~$ python3 -m pip install fastq-statistic
```

# Usage
```bash
user@linux:~$ fastq-statistic --help
Usage: fastq-statistic [OPTIONS] READ1 [READ2]

Arguments:
  READ1    Read1 fastq path or fastq path  [required]
  [READ2]  Read2 filepath or None

Options:
  --sampleid TEXT                 SampleID, default is the first item of
                                  filename splited by underscore(_)

  --result PATH                   Result csv file name, plot with use the same
                                  name.

  --help                          Show this message and exit.

```