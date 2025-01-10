"""
Classifies: CHEBI:36702 2-acyl-1-alkyl-sn-glycero-3-phosphocholine
"""
from rdkit import Chem

def is_2_acyl_1_alkyl_sn_glycero_3_phosphocholine(smiles: str):
    """
    Determines if a molecule is a 2-acyl-1-alkyl-sn-glycero-3-phosphocholine based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is classified as 2-acyl-1-alkyl-sn-glycero-3-phosphocholine, False otherwise
        str: Reason for classification
    """
    # Parse the SMILES string
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define enhanced substructure SMARTS for phosphocholine group
    phosphocholine_pattern = Chem.MolFromSmarts("[O-]P(=O)(O)OCC[N+](C)(C)C")
    if not mol.HasSubstructMatch(phosphocholine_pattern):
        return False, "No phosphocholine group found"

    # Define enhanced substructure SMARTS for sn-1 ether linkage
    ether_chain_pattern = Chem.MolFromSmarts("[C@H](COC)O")
    if not mol.HasSubstructMatch(ether_chain_pattern):
        return False, "No ether linkage with alkyl chain found (sn-1 position)"

    # Define enhanced substructure SMARTS for sn-2 ester linkage
    acyl_group_pattern = Chem.MolFromSmarts("OC(=O)[C,CX4]")
    if not mol.HasSubstructMatch(acyl_group_pattern):
        return False, "No acyl ester group found (sn-2 position)"

    # Optionally, we could also verify the presence of long carbon chains
    carbon_chain_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    if carbon_chain_count < 30:
        return False, "Insufficient carbon content typical for this chemical class"

    return True, "Molecule matches the 2-acyl-1-alkyl-sn-glycero-3-phosphocholine structure"