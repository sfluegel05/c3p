"""
Classifies: CHEBI:28494 cardiolipin
"""
from rdkit import Chem

def is_cardiolipin(smiles: str):
    """
    Determines if a molecule is cardiolipin based on its SMILES string.
    A cardiolipin consists of two phosphatidic acids linked to a glycerol backbone,
    forming a symmetrical molecule with four ester-linked fatty acid chains.
    
    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a cardiolipin, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string into a molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # SMARTS for cardiolipin structure: central glycerol with two phosphate groups and ester linkages
    cardiolipin_pattern = Chem.MolFromSmarts(
        "[C@@H]1([O][P](=O)([O])[O][C@@H]2[C@H](O[P](=O)([O])[O2])[C@H](CO)O2)[C@H](CO[CX3](=O)[O][CX4H1,2])O1"
    )

    # Check for cardiolipin core structure
    if not mol.HasSubstructMatch(cardiolipin_pattern):
        return False, "Cardiolipin core structure not found"

    # Check for four ester-linked fatty acid chains (long alkyl chains connected to ester linkage)
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)[O][C]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 4:
        return False, f"Found {len(ester_matches)} ester linkages, need 4 for cardiolipin"

    return True, "Molecule matches cardiolipin structure"