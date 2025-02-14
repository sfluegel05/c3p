"""
Classifies: CHEBI:55465 7-hydroxyisoflavones
"""
"""
Classifies: CHEBI:50009 7-hydroxyisoflavone
"A hydroxyisoflavone compound having a hydroxy group at the 7-position."
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_7_hydroxyisoflavones(smiles: str):
    """
    Determines if a molecule is a 7-hydroxyisoflavone based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 7-hydroxyisoflavone, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for isoflavone backbone (2 fused rings, one pyrone and one benzene)
    isoflavone_patterns = [
        Chem.MolFromSmarts("c1c(-c2ccc(O)cc2)oc2ccccc2c1=O"),
        Chem.MolFromSmarts("c1c(-c2ccc(O)cc2)oc2ccccn2c1=O"),
        Chem.MolFromSmarts("c1c(-c2ccc(O)cc2)oc2ncccc2c1=O"),
        # Add more patterns as needed
    ]
    if not any(mol.HasSubstructMatch(pattern) for pattern in isoflavone_patterns):
        return False, "No isoflavone backbone found"

    # Look for hydroxy group at 7-position
    hydroxy_7_patterns = [
        Chem.MolFromSmarts("c1c(-c2ccc(O)cc2)oc2cc(O)ccc2c1=O"),
        Chem.MolFromSmarts("c1c(-c2ccc(O)cc2)oc2cc(O)ccn2c1=O"),
        Chem.MolFromSmarts("c1c(-c2ccc(O)cc2)oc2cc(O)ncc2c1=O"),
        # Add more patterns as needed
    ]
    if not any(mol.HasSubstructMatch(pattern) for pattern in hydroxy_7_patterns):
        return False, "No hydroxy group at 7-position"

    # Check that there is only one hydroxy group at 7-position
    hydroxy_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8 and atom.GetTotalNumHs() == 1)
    if hydroxy_count != 1:
        return False, f"Found {hydroxy_count} hydroxy groups, should be exactly 1"

    # Enumerate tautomers and ionization states
    tautomers = AllChem.EnumerateTautomers(mol)
    for tautomer in tautomers:
        if any(tautomer.HasSubstructMatch(pattern) for pattern in hydroxy_7_patterns):
            return True, "Molecule is a 7-hydroxyisoflavone (after tautomer enumeration)"

    return True, "Molecule is a 7-hydroxyisoflavone"