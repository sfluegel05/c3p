"""
Classifies: CHEBI:37739 glycerophospholipid
"""
"""
Classifies: glycerophospholipid
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_glycerophospholipid(smiles: str):
    """
    Determines if a molecule is a glycerophospholipid based on its SMILES string.
    A glycerophospholipid is defined as any glycerolipid having a phosphate group
    ester-linked to a terminal carbon of the glycerol backbone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycerophospholipid, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define glycerol backbone pattern with potential stereochemistry
    # Pattern: Three connected carbon atoms (could be chiral), each not double-bonded
    glycerol_pattern = Chem.MolFromSmarts("[C;!$(C=O)]-[C;!$(C=O)]-[C;!$(C=O)]")
    backbone_matches = mol.GetSubstructMatches(glycerol_pattern)
    if not backbone_matches:
        return False, "No glycerol backbone found"

    # Define phosphate group patterns including different protonation states
    phosphate_patterns = [
        Chem.MolFromSmarts("P(=O)(O)(O)O"),
        Chem.MolFromSmarts("P(=O)(O)(O)[O-]"),
        Chem.MolFromSmarts("P(=O)(O)([O-])[O-]"),
        Chem.MolFromSmarts("P(=O)([O-])([O-])[O-]"),
        Chem.MolFromSmarts("P(=O)(O)(O)OC"),  # Phosphate connected to carbon (e.g., head groups)
    ]

    phosphate_found = False
    for p_pattern in phosphate_patterns:
        if mol.HasSubstructMatch(p_pattern):
            phosphate_found = True
            break
    if not phosphate_found:
        return False, "No phosphate group found"

    # Define patterns for ester and ether linkages at sn-1 and sn-2 positions
    ester_pattern = Chem.MolFromSmarts("C(=O)O[C;!$(C=O)]")
    ether_pattern = Chem.MolFromSmarts("[C;!$(C=O)]O[C;!$(C=O)]")

    ester_matches = mol.GetSubstructMatches(ester_pattern)
    ether_matches = mol.GetSubstructMatches(ether_pattern)

    # Count total fatty acid chains (esters and ethers)
    total_fatty_acid_chains = len(ester_matches) + len(ether_matches)
    if total_fatty_acid_chains < 2:
        return False, f"Found {total_fatty_acid_chains} fatty acid chains, need at least 2"

    # Optional: Check that fatty acid chains are attached to glycerol backbone carbons
    # Map the matches to the backbone carbons
    backbone_atom_indices = set()
    for match in backbone_matches:
        backbone_atom_indices.update(match)
    fatty_acid_linked = 0
    for match in ester_matches + ether_matches:
        for idx in match:
            if idx in backbone_atom_indices:
                fatty_acid_linked +=1
                break  # Only count each chain once
    if fatty_acid_linked < 2:
        return False, "Fatty acid chains are not properly attached to glycerol backbone"

    # Check that phosphate group is attached to a terminal carbon of the glycerol backbone
    phosphate_pattern = Chem.MolFromSmarts("[C;!$(C=O)]OP(=O)(O)O")
    if not mol.HasSubstructMatch(phosphate_pattern):
        # Check for deprotonated phosphate
        phosphate_pattern = Chem.MolFromSmarts("[C;!$(C=O)]OP(=O)(O)[O-]")
        if not mol.HasSubstructMatch(phosphate_pattern):
            return False, "Phosphate group not attached to terminal carbon"

    return True, "Contains glycerol backbone with fatty acid chains and a phosphate group attached"