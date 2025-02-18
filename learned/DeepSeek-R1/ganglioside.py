"""
Classifies: CHEBI:28892 ganglioside
"""
"""
Classifies: CHEBI:28800 ganglioside
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
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for ceramide (amide-linked sphingosine and fatty acid)
    # Amide group: CON with adjacent chains
    amide_pattern = Chem.MolFromSmarts("[CX3](=O)[NX3H][CX4]")
    if not mol.HasSubstructMatch(amide_pattern):
        return False, "No ceramide amide group found"

    # Check for sphingosine (long chain with amino and hydroxyl groups)
    # Simplified: look for adjacent amino and hydroxyl groups in a chain
    sphingo_pattern = Chem.MolFromSmarts("[NH2][CX4H][CX4H]([OH])[CX4H2]")
    if not mol.HasSubstructMatch(sphingo_pattern):
        return False, "No sphingosine backbone detected"

    # Check for oligosaccharide (multiple sugar rings)
    # Look for at least two pyranose rings (approximate)
    pyranose_pattern = Chem.MolFromSmarts("[O]1~C~C~C~C~C~1")
    pyranose_matches = mol.GetSubstructMatches(pyranose_pattern)
    if len(pyranose_matches) < 2:
        return False, f"Only {len(pyranose_matches)} pyranose rings, need at least 2"

    # Check for sialic acid (Neu5Ac: has N-acetyl and carboxylate groups)
    # SMARTS pattern for Neu5Ac core structure
    sialic_acid_pattern = Chem.MolFromSmarts(
        "[C@@H]1([C@H](O)[C@H](O[C@H]2[C@H](O)[C@H](O[C@H]([C@H]2O)CO)NC(=O)C)[C@@H](C(=O)O)O1)"
    )
    if not mol.HasSubstructMatch(sialic_acid_pattern):
        return False, "No sialic acid (Neu5Ac) detected"

    return True, "Contains ceramide, oligosaccharide, and sialic acid(s)"