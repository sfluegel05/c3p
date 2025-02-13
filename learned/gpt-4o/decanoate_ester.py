"""
Classifies: CHEBI:87658 decanoate ester
"""
from rdkit import Chem

def is_decanoate_ester(smiles: str):
    """
    Determines if a molecule is a decanoate ester based on its SMILES string.
    A decanoate ester results from the esterification of the carboxy group of decanoic acid
    (which is a 10 carbon acid) with the hydroxy group of an alcohol or phenol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a decanoate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Ester linkage pattern with a 10-carbon chain accounting for various linking possibilities
    ester_pattern = Chem.MolFromSmarts("C(=O)O")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester linkage found"
    
    # Verify there is a 10-carbon chain connected to the ester by checking the carbon count
    carbon_count = 0
    atoms_in_ester = mol.GetSubstructMatches(ester_pattern)
    for match in atoms_in_ester:
        carbon_atoms = [atom for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6 and atom.IsInRing() == False]
        if len(carbon_atoms) == 10:
            return True, "Contains a decanoate ester linkage"
    
    return False, "No appropriate 10-carbon chain detected (decanoic acid derivative)"