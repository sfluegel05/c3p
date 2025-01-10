"""
Classifies: CHEBI:10283 2-hydroxy fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_2_hydroxy_fatty_acid(smiles: str):
    """
    Classify a molecule as a 2-hydroxy fatty acid based on the SMILES string.
    A 2-hydroxy fatty acid contains a carboxylic acid group and a hydroxyl group on the second carbon of
    a hydrocarbon chain.

    Args:
        smiles (str): SMILES string

    Returns:
        bool: True if a 2-hydroxy fatty acid, False otherwise
        str: Reason for the decision
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Locate a carboxylic acid group at a terminal position
    carboxylic_acid_pattern = Chem.MolFromSmarts("C(=O)O")
    carboxylic_matches = mol.GetSubstructMatches(carboxylic_acid_pattern)
    
    # Accept cases where the carboxylic group is at terminal carbon
    candidate_atoms = [match[0] for match in carboxylic_matches if mol.GetAtomWithIdx(match[0]).GetDegree() == 1]
    
    if not candidate_atoms:
        return False, "No terminal carboxylic acid group found"

    # For each candidate, verify hydroxyl on second carbon
    for terminal_carbon in candidate_atoms:
        neighbors = mol.GetAtomWithIdx(terminal_carbon).GetNeighbors()
        for carbon in neighbors:
            if carbon.GetSymbol() != "C" or carbon.GetIdx() == terminal_carbon:
                continue
            
            # Confirm there is a hydroxyl on adjacent carbon after terminal carbon
            hydroxyl_pattern = Chem.MolFromSmarts("[OX2H]")
            if any(mol.HasSubstructMatch(hydroxyl_pattern, useChirality=True, maxMatches=1, atomIndices=[carbon.GetIdx()])):
                return True, "Carboxylic acid with adjacent hydroxyl set on the second carbon specifically matches the 2-hydroxy fatty acid definition"
    
    return False, "No proper hydroxyl and carboxylic group pairing on a chain to constitute a 2-hydroxy fatty acid"