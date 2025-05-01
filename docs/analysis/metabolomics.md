# Metabolomics Analysis

In order to explore the performance of learned chemical classes
against out-of-distribution data, we performed an enrichment analysis
over all samples in the EBI Metabolights database. On evaluating the
enrichment results against the study metadata using llm-matrix, 78%
were found to be meaningful. Note that this is a fairly crude metric
of validity, as we did not perform any processing of the Metabolights
MAF files, or make use of peak scores. We also did not attempt to
stratify targeted vs untargeted metabolomics, so over-representation
may reflect selection of targets.

Enriched MAF files and ORA analysis can be found here:

[![DOI](https://zenodo.org/badge/13996/cmungall/metabolights-enriched.svg)](https://doi.org/10.5281/zenodo.15233055)

* Notebook: [Metabolights-Overrepresentation-Analysis](../notebooks/Metabolights-Overrepresentation-Analysis.ipynb)
* Source repo: [metabolights-enriched](https://github.com/cmungall/metabolights-enriched)

We used [llm-matrix](https://github.com/monarch-initiative/llm-matrix) to qualitatively assess the enrichment results for each study.

Available here:

 * https://github.com/cmungall/metabolights-enriched/blob/main/analysis/results.csv
 * [google sheet](https://docs.google.com/spreadsheets/d/1Lns4XQpoBveC2bfmfROKDUjojq90hoRziXxMleZRpEY/edit?gid=1819805036#gid=1819805036)




