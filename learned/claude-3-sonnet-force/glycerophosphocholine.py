"""
Classifies: CHEBI:36313 glycerophosphocholine
"""
"""
Classifies: CHEBI:17931 glycerophosphocholine
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_glycerophosphocholine(smiles: str):
    """
    Determines if a molecule is a glycerophosphocholine based on its SMILES string.
    A glycerophosphocholine is a glycerol phosphate ester of a phosphocholine.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a glycerophosphocholine, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone pattern (C-C-C with 2 oxygens attached)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    glycerol_matches = mol.GetSubstructMatches(glycerol_pattern)
    if not glycerol_matches:
        return False, "No glycerol backbone found"
        
    # Look for phosphate group (-O-P(=O)(-O)-O-) and choline group (-CC[N+](C)(C)C) connected to the glycerol backbone
    for glycerol_match in glycerol_matches:
        # Get the atoms in the glycerol backbone
        glycerol_atoms = [mol.GetAtomWithIdx(idx) for idx in glycerol_match]
        
        # Check if the phosphate group is connected to the glycerol backbone
        phosphate_pattern = Chem.MolFromSmarts("[OX2]P([OX1-,OX2-,OX3])(=[OX1])([OX1-,OX2-,OX3])[OX2-,OX1-]"
)
        phosphate_match = mol.GetSubstructMatch(phosphate_pattern)
        if phosphate_match:
            phosphate_atoms = [mol.GetAtomWithIdx(idx) for idx in phosphate_match]
            phosphate_connected = any(atom.GetIdx() in glycerol_match for atom in phosphate_atoms)
        else:
            phosphate_connected = False
        
        # Check if the choline group is connected to the glycerol backbone
        choline_pattern = Chem.MolFromSmarts("CC[N+](C)(C)C")
        choline_match = mol.GetSubstructMatch(choline_pattern)
        if choline_match:
            choline_atoms = [mol.GetAtomWithIdx(idx) for idx in choline_match]
            choline_connected = any(atom.GetIdx() in glycerol_match for atom in choline_atoms)
        else:
            choline_connected = False
        
        # If both phosphate and choline are connected to the glycerol backbone, it's a glycerophosphocholine
        if phosphate_connected and choline_connected:
            # Check for ester group (-O-C(=O)-) to connect glycerol and fatty acid chains
            ester_pattern = Chem.MolFromSmarts("[OX2][CX3](=[OX1])")
            ester_matches = mol.GetSubstructMatches(ester_pattern)
            if len(ester_matches) < 1:
                return False, "No ester groups found"
            
            # Look for fatty acid chains (long carbon chains)
            fatty_acid_pattern = Chem.MolFromSmarts("[CX4,CX3]~[CX4,CX3]~[CX4,CX3]~[CX4,CX3]")
            fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)
            if len(fatty_acid_matches) < 1:
                return False, "Missing fatty acid chains"
            
            # Count rotatable bonds to verify long chains
            n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
            if n_rotatable < 5:
                return False, "Chains too short to be fatty acids"

            return True, "Contains glycerol backbone with phosphocholine and fatty acid chains"

    return False, "Phosphate and choline groups not connected to the same glycerol backbone"