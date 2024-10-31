from rdkit import Chem
from rdkit.Chem import AllChem

def is_sesterterpene(smiles: str):
    """
    Determines if a molecule is a sesterterpene (C25 terpene).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a sesterterpene, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Check carbon count
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C')
    if carbon_count != 25:
        return False, f"Not a C25 compound (contains {carbon_count} carbons)"
        
    # Calculate implicit and explicit H count
    total_h_count = 0
    for atom in mol.GetAtoms():
        total_h_count += atom.GetTotalNumHs()
        
    # Most sesterterpenes should have a formula close to C25H40
    # But we allow some variation due to oxidation, unsaturation etc.
    if total_h_count < 30 or total_h_count > 50:
        return False, f"H count ({total_h_count}) outside typical range for sesterterpenes"
        
    # Check for presence of common non-C/H atoms that might indicate this isn't a terpene
    for atom in mol.GetAtoms():
        if atom.GetSymbol() not in ['C', 'H', 'O']:  # Allow O for some oxidized forms
            return False, f"Contains non-typical atoms for sesterterpenes ({atom.GetSymbol()})"
            
    # If we get here, it's likely a sesterterpene
    return True, f"C25 compound with {total_h_count} hydrogens"
# Pr=1.0
# Recall=0.9591836734693877