"""
Classifies: CHEBI:59845 3-hydroxy fatty acid
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_3_hydroxy_fatty_acid(smiles: str):
    """
    Determines if a molecule is a 3-hydroxy fatty acid based on its SMILES string.
    A 3-hydroxy fatty acid has a carboxylic acid group and a hydroxyl group on the 3rd carbon from the acid.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a 3-hydroxy fatty acid, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES"
    
    # Find carboxylic acid groups [CX3](=O)[OH]
    carboxylic_acid = Chem.MolFromSmarts('[CX3](=O)[OX2H1]')
    ca_matches = mol.GetSubstructMatches(carboxylic_acid)
    if not ca_matches:
        return False, "No carboxylic acid group found"
    
    # Collect carboxylic acid oxygen atoms to exclude from hydroxyl search
    ca_oxygens = {match[2] for match in ca_matches}  # Index of hydroxyl O in COOH
    
    # Find all non-carboxylic hydroxyl groups
    hydroxyl_matches = mol.GetSubstructMatches(Chem.MolFromSmarts('[!C;OX2H1]'))
    for oh_match in hydroxyl_matches:
        o_idx = oh_match[0]
        if o_idx in ca_oxygens:
            continue
        
        # Get attached carbon and check distance from each carboxylic acid
        oh_carbon = mol.GetAtomWithIdx(o_idx).GetNeighbors()[0].GetIdx()
        
        for ca_match in ca_matches:
            carbonyl_c = ca_match[0]  # Carboxylic acid carbonyl carbon index
            try:
                path = Chem.GetShortestPath(mol, carbonyl_c, oh_carbon)
            except:
                continue  # No path exists
            if len(path) == 4:  # 3 bonds between carbonyl and hydroxyl carbons
                return True, "Hydroxyl group found on 3rd carbon from carboxylic acid"
    
    return False, "No 3-hydroxy group detected relative to carboxylic acid"