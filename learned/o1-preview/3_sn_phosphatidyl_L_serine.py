"""
Classifies: CHEBI:11750 3-sn-phosphatidyl-L-serine
"""
"""
Classifies: CHEBI:8339 3-sn-phosphatidyl-L-serine
"""

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_3_sn_phosphatidyl_L_serine(smiles: str):
    """
    Determines if a molecule is a 3-sn-phosphatidyl-L-serine based on its SMILES string.
    A 3-sn-phosphatidyl-L-serine is a glycerophosphoserine compound having acyl substituents 
    at the 1- and 2-hydroxy positions, with correct stereochemistry.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 3-sn-phosphatidyl-L-serine, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Add hydrogens to consider stereochemistry
    mol = Chem.AddHs(mol)

    # Check for glycerol backbone with sn stereochemistry
    # The glycerol backbone with sn stereochemistry can be represented as [C@H](O)[C@@H](O)CO
    glycerol_sn_smarts = '[C@H](O)[C@@H](O)CO'
    glycerol_sn_pattern = Chem.MolFromSmarts(glycerol_sn_smarts)
    if not mol.HasSubstructMatch(glycerol_sn_pattern):
        return False, "No glycerol backbone with sn stereochemistry found"

    # Check for two acyl chains esterified at sn-1 and sn-2 positions
    # Acyl ester pattern: O-C(=O)-C (ester linkage)
    acyl_chain_smarts = 'O[C@H](CO[C@H](CO[P](=O)(O)OC[C@H](N)C(=O)O)OC=O)COC=O'
    acyl_chain_pattern = Chem.MolFromSmarts(acyl_chain_smarts)
    acyl_matches = mol.GetSubstructMatches(acyl_chain_pattern)
    if len(acyl_matches) < 1:
        return False, "Acyl chains esterified at sn-1 and sn-2 positions not found"

    # Check for phosphate group at sn-3 position
    phosphate_smarts = 'COP(=O)(O)OC[C@H](N)C(=O)O'
    phosphate_pattern = Chem.MolFromSmarts(phosphate_smarts)
    if not mol.HasSubstructMatch(phosphate_pattern):
        return False, "Phosphate group at sn-3 position not found"

    # Check for L-serine moiety connected to phosphate group
    serine_smarts = 'OC[C@H](N)C(=O)O'
    serine_pattern = Chem.MolFromSmarts(serine_smarts)
    if not mol.HasSubstructMatch(serine_pattern):
        return False, "L-serine moiety connected to phosphate group not found"

    # Check for overall phosphatidylserine structure
    phosphatidylserine_smarts = '[C@H](COC(=O)[#6])[C@@H](COC(=O)[#6])COP(=O)(O)OC[C@H](N)C(=O)O'
    phosphatidylserine_pattern = Chem.MolFromSmarts(phosphatidylserine_smarts)
    if not mol.HasSubstructMatch(phosphatidylserine_pattern):
        return False, "Overall phosphatidylserine structure not found"

    return True, "Molecule is a 3-sn-phosphatidyl-L-serine"