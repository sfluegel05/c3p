"""
Classifies: CHEBI:36315 glycerophosphoinositol
"""
"""
Classifies: CHEBI:27595 glycerophosphoinositol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_glycerophosphoinositol(smiles: str):
    """
    Determines if a molecule is a glycerophosphoinositol based on its SMILES string.
    A glycerophosphoinositol is a glycerophospholipid with an inositol group attached
    to the phosphate group at the sn-3 position of the glycerol backbone.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycerophosphoinositol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES and handle counterions
    mol = Chem.MolFromSmiles(smiles, removeHs=False)
    if mol is None:
        return False, "Invalid SMILES string"
    mol = Chem.RemoveHs(mol)
    
    # Look for glycerol backbone pattern (C-C-C with 2 oxygen attachments)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
    
    # Look for phosphate group (-O-P(=O)(-O)-O-)
    phosphate_pattern = Chem.MolFromSmarts("O=P(-O)(-O)-O")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    if not phosphate_matches:
        return False, "No phosphate group found"
    
    # Look for inositol ring attached to phosphate (-O-P-O-C1COCC(O)C1)
    inositol_patterns = [
        Chem.MolFromSmarts("[PX4](-[OX2]1[CH2X4][CH1X4]([OH1X3])[CH1X4]([OH1X3])[CH1X4]([OH1X3])[CH1X4]([OH1X3])[CH1X4]([OH1X3])O1)-[OX2]"),
        Chem.MolFromSmarts("[PX4](-[OX2]1[CH1X4]([OH1X3])[CH1X4]([OH1X3])[CH1X4]([OH1X3])[CH1X4]([OH1X3])[CH1X4]([OH1X3])O[CH2X4]1)-[OX2]")
    ]
    inositol_matches = False
    for pattern in inositol_patterns:
        inositol_matches = mol.GetSubstructMatches(pattern)
        if inositol_matches:
            break
    if not inositol_matches:
        return False, "No inositol group found"
    
    # Check for acyl/alkyl chains (long carbon chains attached to glycerol)
    chain_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    chain_matches = mol.GetSubstructMatches(chain_pattern)
    
    # Count rotatable bonds to verify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 4:
        return False, "Chains too short for glycerophosphoinositol"
    
    # Check for correct attachment of chains to glycerol backbone
    glycerol_atoms = set(mol.GetSubstructMatch(glycerol_pattern))
    chain_atoms = set()
    for match in chain_matches:
        chain_atoms.update(match)
    if not glycerol_atoms.intersection(chain_atoms):
        return False, "Acyl/alkyl chains not correctly attached to glycerol backbone"
    
    # Check molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 500:
        return False, "Molecular weight too low for glycerophosphoinositol"
    
    # Decide based on presence of chains
    if not chain_matches:
        return True, "Contains glycerol backbone with inositol group attached to phosphate"
    elif len(chain_matches) == 1:
        return True, "Contains glycerol backbone with inositol group attached to phosphate and one acyl/alkyl chain"
    elif len(chain_matches) == 2:
        return True, "Contains glycerol backbone with inositol group attached to phosphate and two acyl/alkyl chains"
    else:
        return False, "Contains more than two acyl/alkyl chains"