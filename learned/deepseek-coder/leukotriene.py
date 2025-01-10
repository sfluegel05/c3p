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

    # Check for 20 carbons (C20)
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if c_count != 20:
        return False, f"Expected 20 carbons, found {c_count}"

    # Check for 4 double bonds
    double_bonds = sum(1 for bond in mol.GetBonds() if bond.GetBondType() == Chem.BondType.DOUBLE)
    if double_bonds != 4:
        return False, f"Expected 4 double bonds, found {double_bonds}"

    # Check for conjugated system (at least 3 conjugated double bonds)
    conjugated_system = False
    for bond in mol.GetBonds():
        if bond.GetBondType() == Chem.BondType.DOUBLE:
            # Check if the double bond is part of a conjugated system
            if bond.GetBeginAtom().GetIsConjugated() and bond.GetEndAtom().GetIsConjugated():
                conjugated_system = True
                break
    if not conjugated_system:
        return False, "No conjugated system found"

    # Check for carboxylic acid group (-COOH)
    carboxylic_acid_pattern = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
    if not mol.HasSubstructMatch(carboxylic_acid_pattern):
        return False, "No carboxylic acid group found"

    # Check molecular weight (leukotrienes typically have a MW around 300-400 Da)
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300 or mol_wt > 400:
        return False, f"Molecular weight {mol_wt:.2f} Da is outside the expected range for leukotrienes"

    return True, "C20 polyunsaturated fatty acid with 4 double bonds (3 conjugated) and a carboxylic acid group"