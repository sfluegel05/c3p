"""
Classifies: CHEBI:4986 fatty acid methyl ester
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_fatty_acid_methyl_ester(smiles: str):
    """
    Determines if a molecule is a fatty acid methyl ester based on its SMILES string.
    A fatty acid methyl ester is characterized by having a methyl ester subgroup and generally
    a sufficient hydrophobic carbon chain.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if the molecule is a fatty acid methyl ester, False otherwise
        str: Reason for classification
    """
    
    # Parse the SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return (False, "Invalid SMILES string")

    # Check for methyl ester group: C(=O)OC
    methyl_ester_pattern = Chem.MolFromSmarts("C(=O)OC")
    if not mol.HasSubstructMatch(methyl_ester_pattern):
        return False, "No methyl ester group found"
    
    # Check for the presence of a significant carbon chain
    # Here we allow for flexibility, capturing longer chains and branches
    carbon_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    
    # Consider that FAMEs typically have a significant carbon presence
    # (e.g. min of approx 8 carbons as a rough suggestion)
    if carbon_count < 8:
        return False, f"Insufficient carbon content for typical FAME; found {carbon_count} carbons"
    
    # Check for overall properties like hydrophobic chain characteristic
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 3:
        return False, f"Insufficient rotatable bonds, indicating limited chain length or flexibility"

    return True, "Contains methyl ester group with adequate carbon chain typical for a FAME"