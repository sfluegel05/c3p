"""
Classifies: CHEBI:64611 ether lipid
"""
"""
Classifies: CHEBI:35748 ether lipid
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_ether_lipid(smiles: str):
    """
    Determines if a molecule is an ether lipid based on its SMILES string.
    An ether lipid is a lipid with a glycerol backbone where at least one carbon atom is bonded to an alkyl chain via an ether linkage.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an ether lipid, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for glycerol backbone pattern (C-C-C with at least 2 oxygens attached)
    glycerol_pattern = Chem.MolFromSmarts("[CH2X4][CHX4][CH2X4]")
    if not mol.HasSubstructMatch(glycerol_pattern):
        return False, "No glycerol backbone found"
        
    # Look for ether linkages specifically connected to alkyl chains
    # Pattern: [C]~[O]~[C] where at least one C is part of an alkyl chain
    ether_pattern = Chem.MolFromSmarts("[CX4][OX2][CX4]")
    ether_matches = mol.GetSubstructMatches(ether_pattern)
    if len(ether_matches) == 0:
        return False, "No ether linkages found"

    # Check for at least one ether linkage connected to the glycerol backbone
    glycerol_atoms = mol.GetSubstructMatch(glycerol_pattern)
    ether_connected_to_glycerol = False
    for match in ether_matches:
        if any(atom in glycerol_atoms for atom in match):
            # Verify the other end of the ether is part of an alkyl chain
            other_atom = match[0] if match[1] in glycerol_atoms else match[2]
            # Check if the other atom is part of a chain of at least 3 carbons
            chain_pattern = Chem.MolFromSmarts(f"[#{other_atom+1}]~[CX4]~[CX4]")
            if mol.HasSubstructMatch(chain_pattern):
                ether_connected_to_glycerol = True
                break
    if not ether_connected_to_glycerol:
        return False, "No ether linkage connected to both glycerol and alkyl chain"

    # Count rotatable bonds to verify long chains
    n_rotatable = rdMolDescriptors.CalcNumRotatableBonds(mol)
    if n_rotatable < 3:
        return False, "Chains too short to be alkyl chains"

    # Check molecular weight - ether lipids typically >300 Da
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    if mol_wt < 300:
        return False, "Molecular weight too low for ether lipid"

    return True, "Contains glycerol backbone with at least one ether linkage to an alkyl chain"