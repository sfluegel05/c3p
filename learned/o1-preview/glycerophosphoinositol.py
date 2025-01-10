"""
Classifies: CHEBI:36315 glycerophosphoinositol
"""
"""
Classifies: glycerophosphoinositol
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_glycerophosphoinositol(smiles: str):
    """
    Determines if a molecule is a glycerophosphoinositol based on its SMILES string.
    A glycerophosphoinositol is any glycerophospholipid having the polar alcohol inositol esterified
    to the phosphate group at the sn-3 position of the glycerol backbone.

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
    glycerol_phosphate_smarts = "[C@@H](COC(=O)[#6])[C@H](O)COP(=O)(O)O"
    glycerol_phosphate = Chem.MolFromSmarts(glycerol_phosphate_smarts)
    if not mol.HasSubstructMatch(glycerol_phosphate):
        return False, "No glycerol backbone with phosphate at sn-3 position found"

    # Define SMARTS pattern for inositol ring attached to phosphate
    # Inositol is a cyclohexane ring with six hydroxyl groups
    inositol_smarts = "OP(=O)(O[C@H]1[C@@H](O)[C@H](O)[C@@H](O)[C@H](O)[C@@H]1O)O"
    inositol = Chem.MolFromSmarts(inositol_smarts)
    if not mol.HasSubstructMatch(inositol):
        return False, "No inositol group attached to phosphate found"

    # Check for fatty acid chain at sn-1 position (ester or ether linkage)
    ester_sn1_smarts = "OC(=O)C[C@@H](O)CO"
    ester_sn1 = Chem.MolFromSmarts(ester_sn1_smarts)
    ether_sn1_smarts = "OCC[C@@H](O)CO"
    if not (mol.HasSubstructMatch(ester_sn1) or mol.HasSubstructMatch(ether_sn1_smarts)):
        return False, "No fatty acid chain at sn-1 position found"

    # Check for fatty acid chain at sn-2 position (ester linkage)
    ester_sn2_smarts = "[C@@H](OC(=O)[#6])[CH2]O"
    ester_sn2 = Chem.MolFromSmarts(ester_sn2_smarts)
    if not mol.HasSubstructMatch(ester_sn2):
        return False, "No fatty acid chain at sn-2 position found"

    return True, "Molecule is a glycerophosphoinositol"

__metadata__ = {
    'chemical_class': {
        'name': 'glycerophosphoinositol',
        'definition': 'Any glycerophospholipid having the polar alcohol inositol esterified to the phosphate group at the sn-3 position of the glycerol backbone.'
    }
}