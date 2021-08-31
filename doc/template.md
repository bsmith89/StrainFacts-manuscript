# bsmith89/compbio-template2

A batteries-included Snakemake project template for computational biology.

# Creating or cloning a project

## Create a new project

1. Clone this repository
2. Check that that `snakemake` is in your path: `which snakemake`
3. `snakemake -j1 new_project`
4. Follow the "Configure a project" instructions below.

## Clone an existing project

1. Clone the repository
2. `snakemake -j1 initialize_project_config`
3. Follow the "Configure a project" instructions below.
4. Optional: Add the template as a remote:
    `git remote add template-origin <TEMPLATE-URL> && git fetch template-origin template`

## Configure a project
1. Edit the `env_local` script to load your Snakemake environment
    (e.g. `conda activate snakemake-env`)
2. Create a local profile (if necessary) and link it to `profiles/default`
    (e.g. `unlink profiles/default; ln -s profiles/bueno profiles/default`)
3. `source ./env`
4. Be sure to pre-link any outside filesystems before you start running things
    (e.g. `for dir in data raw ref; do ln -s /scratch/my-project-scratch/$dir $dir; done`)

## Merging in template updates
1. Make sure that the `template-origin` remote is loaded and the template
   branch is set to track it.
2. Update your template branch: `git checkout template && git pull`
3. Create a temporary branch for merging: `git checkout main && git checkout -b main-x-template`
4. Merge the template: `git merge template`, and then fix any conflicts.
5. Check what's changing: `git diff main main-x-template`
5. (If satisfied:) Fast-forward the main branch:
    `git checkout main && git merge --ff-only main-x-template`


# TODO:

-   [ ] Make the default conda environment truly generic.
-   [x] Write a base project initialization recipe in `Snakefile`.
-   [ ] Alternatively: switch to [cookiecutter](https://cookiecutter.readthedocs.io/en/latest/).
