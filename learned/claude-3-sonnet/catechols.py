"""
Classifies: CHEBI:33566 catechols
"""
"""
Classifies: CHEBI:27027 catechol
Any compound containing an o-diphenol component
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_catechol(smiles: str):
    """
    Determines if a molecule is a catechol based on its SMILES string.
    A catechol is any compound containing an o-diphenol component (two hydroxyl groups
    attached to adjacent carbon atoms, which may or may not be part of a ring system).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a catechol, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define o-diphenol substructure patterns
    o_diphenol_patterns = [
        Chem.MolFromSmarts("[c$(c(O)c1ccccc1),c$(c(O)c1ccccn1),c$(c(O)c1cccco1),c$(c(O)c1ccncc1)]1[c$(O)c(cccc1),c$(O)c(ccncc1),c$(O)c(cccco1),c$(O)c(ccccn1)]"),
        Chem.MolFromSmarts("c1c(O)c(O)ccc1"),
        Chem.MolFromSmarts("c1c(=O)c(O)cc(O)c1"),
        Chem.MolFromSmarts("c1c(O)c(=O)cc(O)c1"),
        # Add more patterns as needed
    ]

    # Check if molecule contains any of the o-diphenol substructures
    for pattern in o_diphenol_patterns:
        if mol.HasSubstructMatch(pattern):
            return True, "Contains an o-diphenol component"

    # Additional checks
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 100:
        return False, "Molecular weight too low for a catechol"

    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    if o_count < 2:
        return False, "Not enough oxygen atoms for a catechol"

    return False, "Does not contain an o-diphenol component"