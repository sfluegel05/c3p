"""
Classifies: CHEBI:36315 glycerophosphoinositol
"""
"""
Classifies: glycerophosphoinositol
"""
from rdkit import Chem

def is_glycerophosphoinositol(smiles: str):
    """
    Determines if a molecule is a glycerophosphoinositol based on its SMILES string.
    A glycerophosphoinositol is any glycerophospholipid having the polar alcohol inositol
    esterified to the phosphate group at the sn-3 position of the glycerol backbone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycerophosphoinositol, False otherwise
        str: Reason for classification
    """

    # Parse the SMILES string into an RDKit molecule
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for glycerol backbone with phosphate at sn-3 position
    glycerol_phosphate_smarts = "[CH2][CH](O[PX4](=O)(O)O[!H])[CH2]O"  # Glycerol backbone with phosphate at sn-3
    glycerol_phosphate = Chem.MolFromSmarts(glycerol_phosphate_smarts)
    if not mol.HasSubstructMatch(glycerol_phosphate):
        return False, "No glycerol backbone with phosphate at sn-3 position found"

    # Define SMARTS pattern for inositol ring attached to phosphate
    # Inositol is a cyclohexane ring with hydroxyl groups
    inositol_smarts = "O[PX4](=O)(O)OC1[C](O)C(O)C(O)C(O)C1O"  # Phosphate connected to inositol
    inositol = Chem.MolFromSmarts(inositol_smarts)
    if not mol.HasSubstructMatch(inositol):
        # Try a more general pattern without specifying all hydroxyls
        inositol_general_smarts = "O[PX4](=O)(O)OC1CCCCC1"  # Phosphate connected to cyclohexane ring
        inositol_general = Chem.MolFromSmarts(inositol_general_smarts)
        if not mol.HasSubstructMatch(inositol_general):
            return False, "No inositol group attached to phosphate found"

    # Check for fatty acid chains at sn-1 and sn-2 positions (ester or ether linkages)
    # General ester linkage pattern attached to glycerol carbons
    ester_linkage_smarts = "[C](=O)[O][CH][CH][O][CH2]"  # Ester linkage at sn-1 or sn-2
    ester_linkage = Chem.MolFromSmarts(ester_linkage_smarts)
    ester_matches = mol.GetSubstructMatches(ester_linkage)
    if len(ester_matches) < 2:
        # Try to find ether linkages as well
        ether_linkage_smarts = "[O][CH][CH][O][CH2]"  # Ether linkage at sn-1 or sn-2
        ether_linkage = Chem.MolFromSmarts(ether_linkage_smarts)
        ether_matches = mol.GetSubstructMatches(ether_linkage)
        total_chains = len(ester_matches) + len(ether_matches)
        if total_chains < 2:
            return False, f"Found {total_chains} fatty acid chains, expected 2"

    return True, "Molecule is a glycerophosphoinositol"

__metadata__ = {
    'chemical_class': {
        'name': 'glycerophosphoinositol',
        'definition': 'Any glycerophospholipid having the polar alcohol inositol esterified to the phosphate group at the sn-3 position of the glycerol backbone.'
    }
}