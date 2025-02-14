"""
Classifies: CHEBI:35759 1-monoglyceride
"""
"""
Classifies: CHEBI:38636 1-monoglyceride
A monoglyceride in which the acyl substituent is located at position 1.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_1_monoglyceride(smiles: str):
    """
    Determines if a molecule is a 1-monoglyceride based on its SMILES string.
    A 1-monoglyceride is a glycerol backbone with a single fatty acid chain attached via an ester bond at position 1.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a 1-monoglyceride, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone pattern (C-C-C with a secondary alcohol group)
    glycerol_pattern = Chem.MolFromSmarts("[CH2][C]([OH])[CH2]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
        
    # Look for 1 ester group (-O-C(=O)-)
    ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) != 1:
        return False, f"Found {len(ester_matches)} ester groups, need exactly 1"

    # Check that the ester group is attached at position 1 of the glycerol backbone
    ester_atom_idx = ester_matches[0][0]
    ester_atom = mol.GetAtomWithIdx(ester_atom_idx)
    if ester_atom.GetAtomMapNum() != 1:
        return False, "Ester group not attached at position 1"

    # Check for fatty acid chain (long aliphatic chain with no other functional groups)
    fatty_acid_pattern = Chem.MolFromSmarts("[CX4]~[CX4]~[CX4]~[CX4]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
    if len(fatty_acid_matches) < 1:
        return False, "Missing fatty acid chain"

    # Count rotatable bonds to verify long chain
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 5:
        return False, "Chain too short to be fatty acid"

    # Check molecular weight - monoglycerides typically 200-400 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 200 or mol_wt > 400:
        return False, "Molecular weight outside typical range for 1-monoglyceride"

    # Count carbons and oxygens
    c_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 6)
    o_count = sum(1 for atom in mol.GetAtoms() if atom.GetAtomicNum() == 8)
    
    if c_count < 10:
        return False, "Too few carbons for 1-monoglyceride"
    if o_count != 4:
        return False, "Must have exactly 4 oxygens (1 ester group, 2 hydroxy groups)"

    return True, "Contains glycerol backbone with single fatty acid chain attached via ester bond at position 1"