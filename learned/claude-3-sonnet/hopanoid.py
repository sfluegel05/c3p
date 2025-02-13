"""
Classifies: CHEBI:51963 hopanoid
"""
"""
Classifies: CHEBI:35778 Hopanoid
A triterpenoid based on a hopane skeleton.
"""
from rdkit import Chem
from rdkit.Chem import rdFMCS

def is_hopanoid(smiles: str):
    """
    Determines if a molecule is a hopanoid based on its SMILES string.
    A hopanoid is a triterpenoid based on a hopane skeleton, which consists of
    a pentacyclic ring system with the following structure:
        https://en.wikipedia.org/wiki/Hopanoid#/media/File:Hopane_skeleton.svg

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a hopanoid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for pentacyclic ring system characteristic of hopanoids
    hopane_skeleton = Chem.MolFromSmiles("C1CCC2C(C3CCC4C5(CCC(C5CC4C3C2)C)C)C1")
    if not mol.HasSubstructMatch(hopane_skeleton):
        return False, "No pentacyclic hopane skeleton found"

    # Check for common substituents and stereochemistry of hopanoids
    substituted_hopanoids = [
        Chem.MolFromSmarts("[C&r6,r5]([C,O])[C&r6,r5]1[C&r6,r5]2[C&r6,r5][C&r6,r5]3[C&r6,r5]4[C&r6,r5]5[C&r6,r5]([C&r6,r5][C&r6,r5][C&r6,r5]5[C&r6,r5][C&r6,r5][C&r6,r5]4[C&r6,r5][C&r6,r5][C&r6,r5]3[C&r6,r5]2[C&r6,r5][C&r6,r5][C&r6,r5]1[C,O])[C,O]"),
        Chem.MolFromSmarts("[C&r6,r5]1[C&r6,r5][C&r6,r5]2[C&r6,r5][C&r6,r5]3[C&r6,r5]4[C&r6,r5]5[C&r6,r5]([C&r6,r5][C&r6,r5][C&r6,r5]5[C&r6,r5][C&r6,r5][C&r6,r5]4[C&r6,r5][C&r6,r5][C&r6,r5]3[C&r6,r5]2[C&r6,r5][C&r6,r5][C&r6,r5]1[C@@H]([C,O])[C,O])[C,O]"),
        Chem.MolFromSmarts("[C&r6,r5]1[C&r6,r5][C&r6,r5]2[C&r6,r5][C&r6,r5]3[C&r6,r5]4[C&r6,r5]5[C&r6,r5]([C&r6,r5][C&r6,r5][C&r6,r5]5[C&r6,r5][C&r6,r5][C&r6,r5]4[C&r6,r5][C&r6,r5][C&r6,r5]3[C&r6,r5]2[C&r6,r5][C&r6,r5][C&r6,r5]1[C@@H]([C@@H]([C,O])[C,O])[C,O])[C,O]"),
    ]
    if not any(mol.HasSubstructMatch(pattern) for pattern in substituted_hopanoids):
        return False, "Substituent pattern not typical of hopanoids"

    # Check molecular weight - hopanoids typically 400-800 Da
    mol_wt = Chem.rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 400 or mol_wt > 800:
        return False, "Molecular weight outside typical range for hopanoids"

    # Passed all checks
    return True, "Contains pentacyclic hopane skeleton with substituents typical of hopanoids"