"""
Classifies: CHEBI:28892 ganglioside
"""
"""
Classifies: CHEBI:18287 ganglioside
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_ganglioside(smiles: str):
    """
    Determines if a molecule is a ganglioside based on its SMILES string.
    A ganglioside is a glycosphingolipid (ceramide + oligosaccharide) with one or more sialic acids.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a ganglioside, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"

    # Check for ceramide structure: sphingosine (long chain with amine) + fatty acid (amide)
    # Pattern for amide-linked long aliphatic chain
    ceramide_amide = Chem.MolFromSmarts("[NX3][C](=[O])[CX4][CX4,CX3]=[CX3,CX2]")
    if not mol.HasSubstructMatch(ceramide_amide):
        return False, "No ceramide amide group found"

    # Check for oligosaccharide (multiple connected sugars)
    # Look for at least two connected sugar units (e.g., hexoses with multiple OH and ether links)
    sugar_pattern = Chem.MolFromSmarts("[C;!$(C=O)]1O[C@H]([C@@H](O)[C@H](O)[C@H]1O)CO")
    sugar_matches = mol.GetSubstructMatches(sugar_pattern)
    if len(sugar_matches) < 2:
        return False, "Insufficient sugar units for oligosaccharide"

    # Detect sialic acid (Neu5Ac) - core structure: carboxylic acid + N-acetyl
    sialic_acid_pattern = Chem.MolFromSmarts("[C@@]1(O[C@H]([C@H](NC(=O)C)[C@H](O1)CO)[C@H](O)CO)C(=O)O")
    if not mol.HasSubstructMatch(sialic_acid_pattern):
        return False, "No sialic acid residue detected"

    # Ensure sialic acid is attached to the oligosaccharide
    # (Approximate check: sialic acid connected to a sugar oxygen)
    sialic_acid_atoms = mol.GetSubstructMatch(sialic_acid_pattern)
    for atom_idx in sialic_acid_atoms:
        atom = mol.GetAtomWithIdx(atom_idx)
        for neighbor in atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 8 and neighbor.GetDegree() >= 2:  # Oxygen in glycosidic bond
                return True, "Contains ceramide, oligosaccharide, and sialic acid(s)"

    return False, "Sialic acid not properly linked to oligosaccharide"