"""
Classifies: CHEBI:87658 decanoate ester
"""
"""
Classifies: CHEBI:32792 decanoate ester
A fatty acid ester resulting from the formal condensation of the carboxy group of decanoic acid (capric acid) with the hydroxy group of an alcohol or phenol.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_decanoate_ester(smiles: str):
    """
    Determines if a molecule is a decanoate ester based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a decanoate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for decanoate group (CCCCCCCCCCC(=O)O)
    decanoate_pattern = Chem.MolFromSmarts("CCCCCCCCCCC(=O)O")
    decanoate_matches = mol.GetSubstructMatches(decanoate_pattern)
    
    if not decanoate_matches:
        return False, "No decanoate group found"

    # Look for alcohol or phenol group (-O-H or -O-C)
    alcohol_pattern = Chem.MolFromSmarts("[O;H1]")
    phenol_pattern = Chem.MolFromSmarts("[O;X2]C")
    alcohol_matches = mol.GetSubstructMatches(alcohol_pattern)
    phenol_matches = mol.GetSubstructMatches(phenol_pattern)
    
    if not alcohol_matches and not phenol_matches:
        return False, "No alcohol or phenol group found"

    # Check if decanoate and alcohol/phenol groups are connected by an ester bond
    for decanoate_idx in decanoate_matches:
        for alcohol_idx in alcohol_matches:
            if mol.GetBondBetweenAtoms(decanoate_idx[1], alcohol_idx[0]).GetBondType() == Chem.BondType.ESTER:
                return True, "Contains a decanoate group connected to an alcohol or phenol via an ester bond"
        for phenol_idx in phenol_matches:
            if mol.GetBondBetweenAtoms(decanoate_idx[1], phenol_idx[0]).GetBondType() == Chem.BondType.ESTER:
                return True, "Contains a decanoate group connected to an alcohol or phenol via an ester bond"

    return False, "No ester bond found between decanoate and alcohol/phenol groups"