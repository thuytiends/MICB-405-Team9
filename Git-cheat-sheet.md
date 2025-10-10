# Team 9 GitHub Cheat Sheet (Quick Start)

## Why GitHub?
  * It's basically a shared server for us to share files, scripts, notes, etc. :smiley: (e.g. I did an alignment and I want to share the code to use as a group)
  * You can use it directly from our university server by making yourself a clone of the Git repository (below), then you just have a normal directory to work on, and can push the files onto GitHub whenever needed
    - Notes: GitHub can't store large data (unless you pay), so I added a gitignore to exclude BAM, fasta, sequencing files. So  you won't be able to push any of them onto GitHub

## FIRST TIME ONLY: SET UP GIT FROM TERMINAL
  Make sure you clone the github repository in a persistent directory
git clone https://github.com/thuytiends/MICB-405-Team9.git
cd MICB-405-Team9

## Basic Workflow – Every Time You Work

```bash
# 1. Go to your project folder
cd MICB-405-Team9

# 2. Get the latest changes from GitHub
git pull origin main

# 3. Make your changes (edit scripts, add files, etc.)

# 4. Stage the files you changed
git add your_file_here

# 5. Commit with a short message
git commit -m "Describe what you did"

# 6. Push your changes to GitHub
git push origin main

## Common Command
| Task                      | Command                    |
| ------------------------- | -------------------------- |
| Check current status      | `git status`               |
| Add all changes           | `git add .`                |
| View commit history       | `git log --oneline`        |
| See changes before commit | `git diff`                 |
| Undo a file change        | `git checkout -- filename` |

## Git structure overview
├── data/ # Metadata and small data files
├── raw_data/ # Large raw data files (not in Git)
├── scripts/ # Unix bash scripts for alignment and processing
├── results/ # Output files (ignored by Git)
├── R/ # RStudio scripts and analysis
├── docs/ # Meeting notes and documentation
├── README.md # Project overview
├── .gitignore # Files/folders excluded from Git

