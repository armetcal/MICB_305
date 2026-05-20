# MICB 305 Scripts

This repository contains all the code covered in MICB 305 at the University of British Columbia.

# Creating your own GitHub repositories for MICB 305

Each team will be required to maintain a GitHub repository. Guidelines on how to do so will be provided through a Canvas module. A link to the GitHub repo needs to be shared with your TA so they can monitor progress in your repository throughout the term.  In general, the repository should contain the following: 

* A **readme.md** file. Readme.md files tell the reader what the purpose of the repository is. Your readme.md file should include the class (MICB 305) and term, group members, and title that summarizes your project.
* Meeting agendas and minutes, which should be updated on a weekly basis. Most groups add a link to a Google doc in their readme file or create a separate agenda/notes file or folder. Anything is fine as long as it is well organized/easy to find. 
* **Input datasets** for your R analyses (note: Github has a file size limit, so it is not strictly necessary to upload all of your Bash outputs)
* **Annotated and reproducible scripts.** Scripts should be able to be run on another person’s computer from start to finish and produce identical outputs. These should be organized in a logical manner. People should be able to read your paper and immediately find the corresponding code/results for each figure/table.
* **Output files (tables, figures)** in a format that can be visualized (.jpg, .png, .qzv, etc)`

The exact structure of your GitHub repo is up to you, but I recommend the following setup to start:

```text
├── README.md
├── Data/
│   └── input files for R
├── Scripts/
│   ├── QIIME2 code
│   └── R scripts (~1 per figure or 1 per analysis usually works well)
└── Results/
    ├── Figures
    └── Tables
```

Several other files may appear on your local computer when you create an R project and add it to GitHub:
- **.gitignore** (This is a plain text file that will automatically appear - add the relative file paths of any files/folders that you DON'T want to upload to GitHub)
- **.DS_Store** (This saves some settings specific to your computer - best to add this to .gitignore)
- **.Rproj** (Your R Project file, which acts as a manager for your R Project. You never need to edit this file, but it's helpful to add it to GitHub as it'll help all the code run smoothly)
