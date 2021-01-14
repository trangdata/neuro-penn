# 4CE Phase 2.1 Neurological Condition Analysis

This repository performs missing data analysis on data at different sites and reports the results in the [`results/`](results/) folder.

## Which script should I run?

**Short answer: `elix-short-term.Rmd` and `elix-cpns.Rmd` ** 

The best way to run this analysis is to clone this repository on your local machine
(please ensure you're in a directory where you want the repository to be downloaded):

```git clone https://github.com/trang1618/neuro-penn.git```

Then, go inside the repository:

```cd neuro-penn```

and make a copy of `elix-short-term.Rmd`, name it with your site name, for example:

```cp elix-short-term.Rmd elix-short-term-penn.Rmd```

Then, open the R project

```open neuro-penn.Rproj```

and navigate to the newly created file (e.g. `elix-short-term-penn.Rmd`) to modify the code to run on the data at your specific site.
**Please be sure to change the second code chunk to read in your site-specific data, including replacing "penn" with your site name,** for example, `mysite = 'nwu'`.
Once everything runs, please hit the "Knit" button on top of the `.Rmd` file to create an `.html` file that will automatically be put into [`htmls/`](htmls/).

Knit the `elix-cpns-[your_site_id].Rmd` similarly.

Finally, please upload your results (in [`results/`](results/) and [`htmls/`](htmls/)) via a [pull request](https://github.com/trang1618/neuro-penn/pulls) or request @trang1618 to add you as a contributor.

If you run into any problem adapting this code to your data, let us (@meghutch and @trang1618) know via Slack or [submit an issue](https://github.com/trang1618/neuro-penn/issues/new).

## Troubleshooting
### Where is the `run_fish()` function?

It's under [`R/run_fish.R`](R/run_fish.R).

We also have other utility functions in the [`R/`](R/) directory.
If you did not clone the repo, you can simply download these files in the [`R/`](R/) directory and source them in the main notebook.
