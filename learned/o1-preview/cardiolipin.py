"""
Classifies: CHEBI:28494 cardiolipin
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_cardiolipin(smiles: str):
    """
    Determines if a molecule is a cardiolipin based on its SMILES string.
    A cardiolipin consists of a central glycerol connected via phosphodiester bonds
    to two phosphatidic acid units, each with two fatty acid chains.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a cardiolipin, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for cardiolipin
    cardiolipin_smarts = """
    [C@H]([O][P](=O)([O])[O][C@H](COC(=O)[CX4,CX3][CX4,CX3])OC(=O)[CX4,CX3][CX4,CX3])CO[P](=O)([O])[O][C@H](COC(=O)[CX4,CX3][CX4,CX3])OC(=O)[CX4,CX3][CX4,CX3]
    """
    cardiolipin_pattern = Chem.MolFromSmarts(cardiolipin_smarts)
    if cardiolipin_pattern is None:
        return False, "Invalid cardiolipin SMARTS pattern"

    # Check if molecule matches the cardiolipin pattern
    if not mol.HasSubstructMatch(cardiolipin_pattern):
        return False, "Molecule does not match cardiolipin structural pattern"

    # Additional checks for fatty acid chains length (optional)
    fatty_acid_chains = []
    ester_bonds = Chem.MolFromSmarts("C(=O)O[C@H]C")
    for match in mol.GetSubstructMatches(ester_bonds):
        fatty_acid_chain = []
        carbon = mol.GetAtomWithIdx(match[0])
        # Traverse the fatty acid chain
        while True:
            neighbors = [a for a in carbon.GetNeighbors() if a.GetAtomicNum() == 6]
            if not neighbors or len(neighbors) == 0:
                break
            carbon = neighbors[0]
            fatty_acid_chain.append(carbon.GetIdx())
            # Break if chain is unreasonably long
            if len(fatty_acid_chain) > 30:
                break
        fatty_acid_chains.append(len(fatty_acid_chain))

    if len(fatty_acid_chains) != 4:
        return False, f"Expected 4 fatty acid chains, found {len(fatty_acid_chains)}"

    # Check that there are two phosphate groups connected via phosphodiester bonds
    phosphodiester_pattern = Chem.MolFromSmarts("P(=O)(O)[O][C]")
    phosphodiester_matches = mol.GetSubstructMatches(phosphodiester_pattern)
    if len(phosphodiester_matches) < 2:
        return False, "Phosphodiester bonds not found or insufficient"

    return True, "Molecule matches cardiolipin structural features"