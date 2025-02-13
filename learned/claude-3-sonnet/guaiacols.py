"""
Classifies: CHEBI:134251 guaiacols
"""
"""
Classifies: CHEBI:27090 guaiacols
Any phenol carrying an additional methoxy substituent at the ortho-position.
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_guaiacols(smiles: str):
    """
    Determines if a molecule is a guaiacol based on its SMILES string.
    A guaiacol is any phenol carrying an additional methoxy substituent at the ortho-position.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a guaiacol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for guaiacol substructure pattern
    guaiacol_pattern = Chem.MolFromSmarts("Oc1ccc(OC)cc1")
    if not mol.HasSubstructMatch(guaiacol_pattern):
        return False, "No guaiacol substructure found"
    
    # Count methoxy and hydroxyl groups
    methoxy_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == "O" and sum(bond.GetBondTypeAsDouble() for bond in atom.GetBonds()) == 1)
    hydroxyl_count = sum(1 for atom in mol.GetAtoms() if atom.GetSymbol() == "O" and sum(bond.GetBondTypeAsDouble() for bond in atom.GetBonds()) == 2)
    if methoxy_count != 1 or hydroxyl_count != 1:
        return False, "Incorrect number of methoxy or hydroxyl groups"
    
    # Check molecular weight range
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 100 or mol_wt > 500:
        return False, "Molecular weight outside typical range for guaiacols"
    
    return True, "Contains a phenol ring with a methoxy substituent at the ortho position"