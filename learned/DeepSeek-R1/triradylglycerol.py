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
    A triradylglycerol has a glycerol backbone with three substituents (acyl, alkyl, or alkenyl).

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
    backbone_pattern = Chem.MolFromSmarts('[CH2]([O][!H])-[CH]([O][!H])-[CH2]([O][!H])')
    matches = mol.GetSubstructMatches(backbone_pattern)
    if not matches:
        return False, "No glycerol backbone with three O-linked groups"

    # Extract O indices from the first match (positions 1,3,5 in SMARTS match tuple)
    try:
        match = matches[0]  # Use first match
        o_indices = [match[1], match[3], match[5]]
    except IndexError:
        return False, "Invalid backbone structure"

    substituent_types = []
    for o_idx in o_indices:
        o_atom = mol.GetAtomWithIdx(o_idx)
        if o_atom.GetDegree() != 1:
            return False, f"Oxygen {o_idx} has multiple bonds"
        
        neighbor = o_atom.GetNeighbors()[0]  # Carbon attached to O
        if neighbor.GetAtomicNum() != 6:
            return False, f"Non-carbon substituent at O {o_idx}"

        # Check for acyl (ester: O-C=O)
        acyl = False
        for bond in neighbor.GetBonds():
            if bond.GetBondType() == Chem.BondType.DOUBLE:
                other_atom = bond.GetOtherAtom(neighbor)
                if other_atom.GetAtomicNum() == 8:  # Double bond to oxygen
                    acyl = True
                    break
        if acyl:
            substituent_types.append('acyl')
            continue

        # Check for alkenyl (vinyl ether: O-C=C)
        alkenyl = False
        for bond in neighbor.GetBonds():
            if bond.GetBondType() == Chem.BondType.DOUBLE and bond.GetOtherAtom(neighbor).GetAtomicNum() == 6:
                alkenyl = True
                break
        if alkenyl:
            substituent_types.append('alkenyl')
            continue

        # Default to alkyl (ether) if no double bonds
        substituent_types.append('alkyl')

    # Validate all substituent types
    valid_types = {'acyl', 'alkyl', 'alkenyl'}
    if all(st in valid_types for st in substituent_types):
        return True, f"Triradylglycerol with substituents: {substituent_types}"
    else:
        return False, f"Invalid substituents: {substituent_types}"