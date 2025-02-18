"""
Classifies: CHEBI:76579 triradylglycerol
"""
"""
Classifies: CHEBI:17855 triradylglycerol
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_triradylglycerol(smiles: str):
    """
    Determines if a molecule is a triradylglycerol based on its SMILES string.
    A triradylglycerol has a glycerol backbone with three substituents, each being
    acyl (ester), alkyl (ether), or alk-1-enyl (vinyl ether) groups.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a triradylglycerol, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Check for glycerol backbone with three O-linked substituents
    backbone_pattern = Chem.MolFromSmarts('[CH2]-[CH](-O-!@*)-[CH2]')
    if not mol.HasSubstructMatch(backbone_pattern):
        return False, "No glycerol backbone with three O-linked groups"

    # Get O atoms from the first backbone match
    matches = mol.GetSubstructMatches(backbone_pattern)
    if not matches:
        return False, "No backbone match found"
    
    # Extract O indices (positions 2, 4, 6 in the SMARTS [CH2:1]-[CH:2](-O:3-!@*)-[CH2:4])
    try:
        o_indices = [match[2], match[5], match[8]]  # Adjust indices based on SMARTS structure
    except IndexError:
        return False, "Invalid backbone structure"

    substituent_types = []
    for o_idx in o_indices:
        o_atom = mol.GetAtomWithIdx(o_idx)
        # Verify oxygen has single bond to carbon
        if o_atom.GetDegree() != 1:
            return False, f"Oxygen {o_idx} has multiple bonds"
        neighbor = o_atom.GetNeighbors()[0]
        if neighbor.GetAtomicNum() != 6:
            return False, f"Substituent at O {o_idx} is not carbon-based"

        # Check for ester (O-C=O)
        ester = False
        for bond in neighbor.GetBonds():
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                other_atom = bond.GetOtherAtom(neighbor)
                if other_atom.GetAtomicNum() == 8:  # Connected to oxygen
                    ester = True
                    break
        if ester:
            substituent_types.append('acyl')
            continue

        # Check for vinyl ether (O-C=C)
        vinyl = False
        for bond in neighbor.GetBonds():
            if bond.GetBondType() == Chem.BondType.DOUBLE and bond.GetOtherAtom(neighbor).GetAtomicNum() == 6:
                vinyl = True
                break
        if vinyl:
            substituent_types.append('alkenyl')
            continue

        # Default to alkyl if no other matches
        substituent_types.append('alkyl')

    # Verify all substituents are valid types
    if all(st in ['acyl', 'alkyl', 'alkenyl'] for st in substituent_types):
        return True, f"Triradylglycerol with substituents: {substituent_types}"
    else:
        return False, f"Invalid substituent types: {substituent_types}"