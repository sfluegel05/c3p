"""
Classifies: CHEBI:37739 glycerophospholipid
"""
Here's an improved version of the code to classify glycerophospholipids:

"""
Classifies: CHEBI:18059 glycerophospholipid
Any glycerolipid having a phosphate group ester-linked to a terminal carbon 
of the glycerol backbone.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_glycerophospholipid(smiles: str):
    """
    Determines if a molecule is a glycerophospholipid based on its SMILES string.

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
    
    # Look for glycerol backbone pattern
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][C@H]([OX2])[C@H]([OX2])[CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
    
    # Look for phosphate group connected to glycerol backbone
    phosphate_pattern = Chem.MolFromSmarts("P(~O)(~O)(~O)")
    phosphate_matches = mol.GetSubstructMatches(phosphate_pattern)
    
    for match in phosphate_matches:
        phosphate_atom = mol.GetAtomWithIdx(match[0])
        for bond in phosphate_atom.GetBonds():
            neighbor = bond.GetOtherAtom(phosphate_atom)
            if neighbor.IsInRing() and neighbor.GetRingSize() == 3:  # Check for glycerol ring
                break
        else:
            continue  # Phosphate not connected to glycerol backbone
        break
    else:
        return False, "Phosphate group not connected to glycerol backbone"
    
    # Look for ester groups (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) < 1:
        return False, "No ester groups found"
    
    # Check for fatty acid chains (long carbon chains attached to esters)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 1:
        return False, "Missing fatty acid chains"
    
    # Count rotatable bonds to verify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 10:
        return False, "Chains too short to be fatty acids"
    
    return True, "Contains glycerol backbone with phosphate and fatty acid chains"