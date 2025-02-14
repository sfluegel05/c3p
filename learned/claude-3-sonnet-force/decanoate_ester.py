"""
Classifies: CHEBI:87658 decanoate ester
"""
"""
Classifies: CHEBI:34892 decanoate ester
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

    # Look for decanoate group (-C(=O)CCCCCCCCC)
    decanoate_pattern = Chem.MolFromSmarts("C(=O)CCCCCCCCC")
    decanoate_matches = mol.GetSubstructMatches(decanoate_pattern)
    if len(decanoate_matches) == 0:
        return False, "No decanoate group found"

    # Look for alcohol/phenol group (-O-[!#6])
    alcohol_pattern = Chem.MolFromSmarts("O[!#6]")
    alcohol_matches = mol.GetSubstructMatches(alcohol_pattern)
    if len(alcohol_matches) == 0:
        return False, "No alcohol or phenol group found"

    # Check if the decanoate group is connected to an alcohol/phenol group via an ester bond
    for decanoate_idx in decanoate_matches:
        for alcohol_idx in alcohol_matches:
            if any(mol.GetBondBetweenAtoms(decanoate_idx[1], alcohol_idx[0]).GetBondType() == Chem.BondType.ESTER):
                return True, "Contains a decanoate group linked to an alcohol or phenol via an ester bond"

    return False, "Decanoate group not linked to an alcohol or phenol via an ester bond"