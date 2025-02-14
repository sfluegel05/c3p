"""
Classifies: CHEBI:25106 macrolide
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_macrolide(smiles: str):
    """
    Determines if a molecule is a macrolide based on its SMILES string.
    Macrolides are characterized by a macrocyclic lactone ring (12+ atoms) derived from polyketides.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a macrolide, False otherwise
        str: Reason for classification
    """
    
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Find lactone group
    lactone_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2]")
    lactone_matches = mol.GetSubstructMatches(lactone_pattern)

    if not lactone_matches:
          return False, "No lactone group found"
    
    for match in lactone_matches:
        lactone_carbon_index = match[0]
        
        # Find rings containing the lactone carbon
        rings = mol.GetAtomWithIdx(lactone_carbon_index).GetOwningRings()

        for ring in rings:
            if len(ring) >= 12:
                return True, f"Macrocyclic lactone ring (size: {len(ring)}) found."
            
    return False, "No macrocyclic lactone ring (12+ atoms) found."