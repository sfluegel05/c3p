"""
Classifies: CHEBI:22315 alkaloid
"""
"""
Classifies: CHEBI:22333 alkaloid
"""
from rdkit import Chem
from rdkit.Chem import Descriptors

def is_alkaloid(smiles: str):
    """
    Determines if a molecule is an alkaloid based on its SMILES string.
    Alkaloids are nitrogen-containing compounds, typically with heterocyclic structures,
    excluding amino acids, peptides, proteins, nucleotides, and simple amines.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alkaloid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Check for presence of nitrogen atoms
    nitrogens = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 7]
    if not nitrogens:
        return False, "No nitrogen atoms present"

    # Check for amino acids (carboxylic acid + amine group)
    carboxylic_acid = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    amine = Chem.MolFromSmarts("[NH2]")
    if mol.HasSubstructMatch(carboxylic_acid) and mol.HasSubstructMatch(amine):
        return False, "Amino acid detected"

    # Check for nucleotides (simplified: phosphate + sugar-like structure)
    phosphate = Chem.MolFromSmarts("[PX4]")
    sugar = Chem.MolFromSmarts("C1OCC(O)C1")  # Rough sugar pattern
    if mol.HasSubstructMatch(phosphate) and mol.HasSubstructMatch(sugar):
        return False, "Nucleotide detected"

    # Check if any nitrogen is in a ring (heterocyclic)
    ring_info = mol.GetRingInfo()
    for n in nitrogens:
        if ring_info.NumAtomRings(n.GetIdx()) > 0:
            return True, "Nitrogen atom in a heterocyclic ring"

    # Check for complex structures with rings and nitrogen (even if not in ring)
    if ring_info.NumRings() > 0:
        return True, "Complex structure with rings and nitrogen"

    # Exclude small amines (primary/secondary with low molecular weight)
    primary_secondary = Chem.MolFromSmarts("[NH2,NH1]")
    if mol.HasSubstructMatch(primary_secondary):
        mol_wt = Descriptors.ExactMolWt(mol)
        if mol_wt < 200:  # Arbitrary threshold for simple amines
            return False, "Simple low molecular weight amine"

    return True, "Nitrogen-containing compound meeting alkaloid criteria"