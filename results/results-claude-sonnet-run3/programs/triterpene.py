from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_triterpene(smiles: str):
    """
    Determines if a molecule is a triterpene (C30 terpene).

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a triterpene, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"

    # Check carbon count = 30
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == 'C')
    if carbon_count != 30:
        return False, f"Carbon count is {carbon_count}, not 30"

    # Check for terpene characteristics:
    # - Mostly hydrocarbons (C,H)
    # - May contain some O atoms
    # - No other elements typically present
    allowed_elements = {'C', 'H', 'O'}
    for atom in mol.GetAtoms():
        if atom.GetSymbol() not in allowed_elements:
            return False, f"Contains non-terpene element: {atom.GetSymbol()}"

    # Count double bonds
    double_bond_count = len([bond for bond in mol.GetBonds() if bond.GetBondType() == Chem.BondType.DOUBLE])

    # Look for rings which are common in triterpenes
    ring_info = mol.GetRingInfo()
    num_rings = ring_info.NumRings()
    
    if num_rings > 0:
        return True, f"C30 molecule with {num_rings} rings and {double_bond_count} double bond(s), consistent with triterpene structure"
    elif double_bond_count > 0:
        return True, f"C30 molecule with {double_bond_count} double bond(s), consistent with triterpene structure"
    else:
        return True, "Saturated C30 molecule consistent with triterpene structure"
# Pr=1.0
# Recall=1.0