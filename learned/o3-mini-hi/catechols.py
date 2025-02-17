"""
Classifies: CHEBI:33566 catechols
"""
#!/usr/bin/env python3
"""
Classifies: Catechols (Any compound containing an o-diphenol component)

A catechol is defined as a molecule containing a benzene ring with two hydroxyl (-OH)
substituents in adjacent (ortho) positions. In this implementation we use a SMARTS pattern
that represents 1,2-dihydroxybenzene. Note that this pattern will match a benzene ring
with adjacent –OH groups even if there are extra substituents on the ring.
"""

from rdkit import Chem

def is_catechols(smiles: str):
    """
    Determines if a molecule is a catechol (contains an o-diphenol component)
    based on its SMILES string. Here a catechol is taken to be any benzene ring
    carrying two adjacent hydroxyl substituents.

    Args:
        smiles (str): SMILES string of the molecule.

    Returns:
        bool: True if the molecule contains a catechol moiety, False otherwise.
        str: Reason for the classification.
    """
    # Parse the SMILES string into an RDKit molecule.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # For substructure searches it is best if all hydrogens are explicit.
    mol = Chem.AddHs(mol)
    
    # Define a SMARTS pattern for a catechol substructure.
    # This pattern looks for a benzene ring (6-membered aromatic ring) that has
    # two adjacent -OH substituents. (The pattern "c1c(O)c(O)cc1" means that
    # on a ring numbered 1, atom 1 is aromatic (c) then atom 2 is aromatic and
    # carries an -OH, then atom 3 is aromatic and carries an -OH, followed by two more aromatic atoms,
    # closing the ring.)
    catechol_smarts = "c1c(O)c(O)cc1"
    catechol_query = Chem.MolFromSmarts(catechol_smarts)
    if catechol_query is None:
        return False, "Error in generating catechol SMARTS pattern"
    
    # Use the SMARTS substructure search to detect a catechol moiety.
    if mol.HasSubstructMatch(catechol_query):
        return True, "Contains o-diphenol (catechol) moiety on a benzene ring"
    else:
        return False, "No adjacent hydroxyl groups on a benzene ring found"

# Example testing (uncomment when running standalone):
# examples = [
#     ("[C@@H]1([C@@H]([C@H]([C@@H]([C@H](O1)CO)O)O)O)O[C@H]2/C(/[C@](C(=CO2)C(=O)OC)([H])CC(=O)OCCC=3C=CC(=C(C3)O)O)=C/C", "oleuropein"),
#     ("C1=CC(=C(C(=C1O)O)[N+]([O-])=O)C", "4-methyl-3-nitrocatechol"),
#     ("C1(=CC=C(C=C1O)[N+]([O-])=O)O", "4-nitrocatechol"),
#     ("C=1(C=CC(=C(C1)O)O)/C=C/C(OCC)=O", "ethyl trans-caffeate"),
#     ("S(OC1=C(O)C=C([C@@H](O)CN)C=C1)(O)(=O)=O", "norepinephrine sulfate"),
# ]
#
# for smi, name in examples:
#     result, reason = is_catechols(smi)
#     print(f"NAME: {name} -> {result}: {reason}")