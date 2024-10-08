# Slurm tips and tricks

Get an estimate of when your job might start

```
$ sbatch --test-only ...
sbatch: You are submitting jobs under your home directory.
sbatch: We recommend to submit jobs under your project directory or scratch directory.
sbatch: Job 55499212 to start at 02:59:39 12/02/24 using 1 processors on nodes spartan-bm005 in partition cascade
```

---

Check the [fairshare] of your project

```
$ sshare -l | sed -e '1p' -e '/punim2009/!d'
Account                    User  RawShares  NormShares    RawUsage   NormUsage  EffectvUsage  FairShare    LevelFS                    GrpTRESMins                    TRESRunMins
 punim2009                               1    0.000469     1625495    0.000021      0.000021             22.061798                                cpu=0,mem=0,energy=0,node=0,b+
  punim2009              mihall          1    0.500000     1625495    0.000021      1.000000   0.228625   0.500000                                cpu=0,mem=0,energy=0,node=0,b+
```

---

Detailed job [priority](https://dashboard.hpc.unimelb.edu.au/scheduler/#job-priority)

```
sprio -j 46085547
```

---

Position in the queue

```
squeue -p cascade -t pending
```

---

Cancel all (of your) jobs

```
scancel -u $USER
```

---

Cancel all jobs whose job name matches a pattern


```
squeue -u $USER -ho "%i %j" | grep pattern | cut -f1 -d' ' | xargs scancel
```

---

Downloading files from mediaflux shareable links on spartan

```
module load mediaflux-data-mover
token=<token from shareable link>
mediaflux-data-mover-cli -download $token <outdir>
```

alternatively

```
link="https://mediaflux.researchsoftware.unimelb.edu.au:443/mflux/share.mfjp?_token=<TOKEN>&browser=true&filename=<NAME>.zip"
curl "$link" -d browser=false -o <filename>
```

---

Check the status of all your jobs from the last two days

```
sacct -S now-2days
```

---

Get the full job name

```
sacct --format="JobID,JobName%150" -j <jobid>
```

---

Count the number of different job states in your job history

```
sacct -nb | tr -s ' ' | cut -d' ' -f2 | sort | uniq -c
```

`-n` means no header, `-b` means brief, `-s` in `tr` compresses runs of space into a single space

---

[fairshare]: https://slurm.schedmd.com/classic_fair_share.html
