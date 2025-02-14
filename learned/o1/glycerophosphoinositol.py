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
    # Remove stereochemistry to make it more general
    glycerol_phosphate_smarts = "C(CO[PX4](=O)(O)O)CO"
    glycerol_phosphate = Chem.MolFromSmarts(glycerol_phosphate_smarts)
    if not mol.HasSubstructMatch(glycerol_phosphate):
        return False, "No glycerol backbone with phosphate at sn-3 position found"

    # Define SMARTS pattern for inositol ring attached to phosphate
    # Inositol is a cyclohexane ring with six hydroxyl groups connected via the phosphate
    inositol_smarts = "O[PX4](=O)(O)-O[C@@H]1[C@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@@H]1O"
    inositol = Chem.MolFromSmarts(inositol_smarts)
    if not mol.HasSubstructMatch(inositol):
        # Try a general pattern without stereochemistry
        inositol_smarts_general = "O[P](=O)(O)-OC1CCCC(C1O)O"
        inositol_general = Chem.MolFromSmarts(inositol_smarts_general)
        if not mol.HasSubstructMatch(inositol_general):
            return False, "No inositol group attached to phosphate found"

    # Check for fatty acid chains at sn-1 and sn-2 positions (ester or ether linkages)
    # General ester linkage pattern
    ester_sn1_sn2_smarts = "[#6][CX4H]([OX2H])[CX4H2][OX2H]"
    ester_sn1_sn2 = Chem.MolFromSmarts(ester_sn1_sn2_smarts)
    matches = mol.GetSubstructMatches(ester_sn1_sn2)
    if len(matches) < 2:
        return False, f"Found {len(matches)} fatty acid chains, expected 2"

    return True, "Molecule is a glycerophosphoinositol"

__metadata__ = {
    'chemical_class': {
        'name': 'glycerophosphoinositol',
        'definition': 'Any glycerophospholipid having the polar alcohol inositol esterified to the phosphate group at the sn-3 position of the glycerol backbone.'
    }
}