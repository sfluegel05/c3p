"""
Classifies: CHEBI:61384 sulfolipid
"""
"""
Classifies: CHEBI:35967 sulfolipid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_sulfolipid(smiles: str):
    """
    Determines if a molecule is a sulfolipid based on its SMILES string.
    A sulfolipid is defined as a compound containing a sulfonic acid residue
    joined by a carbon-sulfur bond to a lipid.

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
    
    # Look for sulfonic acid group (-S(=O)(=O)O)
    sulfonic_acid_pattern = Chem.MolFromSmarts("S(=O)(=O)O")
    sulfonic_acid_matches = mol.GetSubstructMatches(sulfonic_acid_pattern)
    if not sulfonic_acid_matches:
        return False, "No sulfonic acid group found"
    
    # Check if sulfonic acid group is attached to lipid chain
    for s_idx in sulfonic_acid_matches:
        sulfur_atom = mol.GetAtomWithIdx(s_idx)
        for bond in sulfur_atom.GetBonds():
            if bond.GetOtherAtom(sulfur_atom).GetAtomicNum() == 6:  # Carbon
                carbon_atom = bond.GetOtherAtom(sulfur_atom)
                lipid_chain_pattern = Chem.MolFromSmarts(f"[CX4]~[CX4]~[CX4]~[CX4]~{[carbon_atom.GetIdx()]}")
                lipid_chain_matches = mol.GetSubstructMatches(lipid_chain_pattern)
                if lipid_chain_matches:
                    return True, "Contains sulfonic acid group joined by a carbon-sulfur bond to a lipid chain"
    
    return False, "Sulfonic acid group not connected to lipid chain via carbon-sulfur bond"