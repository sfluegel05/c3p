# C3PO Benchmark

The benchmark dataset used to evaluate the performance of C3PO is
derived from [ChEBI](http://obofoundry.org/ontology/chebi), like most
ontologies, consists entirely of entries that are modeled as classes
in the OWL-DL representation.

Here we use an alternative system, the [CHEMROF
framework](https://w3id.org/chemrof/) which allows a candidate
partitioning of ChEBI into grouping classes (or simply “classes”, in
this manuscript), and chemical entities (the things to be classified,
also called structures here).

We use the following to partition ChEBI entries into two groupings, structures and classes:

 * If an entry has a SMILES string, it is included as a structure. These are included in ChEBI as property_value: http://purl.obolibrary.org/obo/chebi/smiles. We eliminate SMILES strings that have wildcards (“*”), as well as CHEBI entries with subclasses, focusing only on ground structures.
 * Otherwise, if the entry is a direct or indirect is-a parent (superclass) of at least one structure, it is included as a class. All structures grouped in this way are treated as positive examples, and all other structures count as negative examples.
 * Other entries are discarded.

## Availability on Hugging Face

* [datasets/MonarchInit/C3PO](https://huggingface.co/datasets/MonarchInit/C3PO)

## Summary Statistics

* See [notebooks/Summarize-Benchmark](notebooks/Summarize-Benchmark.ipynb)
