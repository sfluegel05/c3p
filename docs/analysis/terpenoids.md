# Terpenoid vision model analysis

Based on our analysis of the NPAtlas database, we sought to investigate the classification of terpenoids further. We took 107 terpenoid structures where the C3PO classification disagreed with what was in the ChEBI benchmarks. For each structure, we generated a 2D visualization using RDKit, and presented the image to a multimodal LLM (GPT-4o) and prompted the model to classify the structure into a terpenoid subtype, or another category. The model agreed with the ChEBI classification in only 40% of cases.

- Results as [google sheet](https://docs.google.com/spreadsheets/d/1lqHS2DSKax6TwbsBibN2Sw1bC6CJRwV3sYqk9FwUYGQ/edit?gid=1854221255#gid=1854221255)
- Results as [HTML table](terpenoids_table.html)

See also [issue #4651](https://github.com/ebi-chebi/ChEBI/issues/4651)


