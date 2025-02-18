"""
Classifies: CHEBI:61384 sulfolipid
"""
"""
Classifies: CHEBI:73404 sulfolipid
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_sulfolipid(smiles: str):
    """
    Determines if a molecule is a sulfolipid based on its SMILES string.
    A sulfolipid contains a sulfate group (O-SO3) attached to a lipid structure.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sulfolipid, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for sulfate group (O-SO3)
    sulfate_pattern = Chem.MolFromSmarts("[OX2]S(=O)(=O)[OX2]")
    if not mol.HasSubstructMatch(sulfate_pattern):
        return False, "No sulfate group found"

    # Check lipid characteristics: long chains and ester/amide groups
    # Minimum 14 carbons in a chain (common in fatty acids)
    long_chain = Chem.MolFromSmarts("[CH2]~[CH2]~[CH2]~[CH2]~[CH2]~[CH2]~[CH2]~[CH2]~[CH2]~[CH2]~[CH2]~[CH2]~[CH2]~[CH2]")
    if not mol.HasSubstructMatch(long_chain):
        return False, "No sufficiently long hydrocarbon chain (minimum 14 carbons)"

    # Check for ester or amide bonds characteristic of lipids
    ester = Chem.MolFromSmarts("[O][C]=O")
    amide = Chem.MolFromSmarts("[N][C]=O")
    if not mol.HasSubstructMatch(ester) and not mol.HasSubstructMatch(amide):
        return False, "No ester/amide groups found"

    # Verify significant lipid-like molecular weight (>500 Da typical)
    mol_wt = AllChem.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, f"Molecular weight too low ({mol_wt:.1f} Da) for lipid"

    return True, "Contains sulfate group attached to lipid structure with ester/amide bonds"