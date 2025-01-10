"""
Classifies: CHEBI:25029 leukotriene
"""
"""
Classifies: CHEBI:36044 leukotriene
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_leukotriene(smiles: str):
    """
    Determines if a molecule is a leukotriene based on its SMILES string.
    A leukotriene is a C20 polyunsaturated fatty acid with four double bonds, three of which are conjugated.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a leukotriene, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for approximately 20 carbons (C20, but allow for some variation)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if not (18 <= c_count <= 30):  # Relaxed range to accommodate derivatives
        return False, f"Expected around 20 carbons, found {c_count}"

    # Check for at least 4 double bonds (allow for more in derivatives)
    double_bonds = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.BondType.DOUBLE)
    if double_bonds < 4:
        return False, f"Expected at least 4 double bonds, found {double_bonds}"

    # Check for at least 3 conjugated double bonds in a specific arrangement
    conjugated_double_bonds = 0
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.DOUBLE and bond.GetIsConjugated():
            conjugated_double_bonds += 1
    if conjugated_double_bonds < 3:
        return False, f"Expected at least 3 conjugated double bonds, found {conjugated_double_bonds}"

    # Check for carboxylic acid group or its derivatives (e.g., esters, amides, etc.)
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H0OX1H0]")
    ester_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2][CX4]")
    amide_pattern = Chem.MolFromSmarts("[CX3](=O)[NX3]")
    
    if not (mol.HasSubstructMatch(carboxylic_acid_pattern) or 
            mol.HasSubstructMatch(ester_pattern) or 
            mol.HasSubstructMatch(amide_pattern)):
        return False, "No carboxylic acid group or derivative found"

    # Check molecular weight (leukotrienes typically have a MW around 300-500 Da, allowing for derivatives)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 250 or mol_wt > 800:  # Relaxed range to accommodate derivatives
        return False, f"Molecular weight {mol_wt:.2f} Da is outside the expected range for leukotrienes"

    return True, "C20 polyunsaturated fatty acid with at least 4 double bonds (3 conjugated) and a carboxylic acid group or derivative"