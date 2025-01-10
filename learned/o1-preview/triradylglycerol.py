"""
Classifies: CHEBI:76579 triradylglycerol
"""
"""
Classifies: CHEBI:<ID> triradylglycerol
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_triradylglycerol(smiles: str):
    """
    Determines if a molecule is a triradylglycerol based on its SMILES string.
    A triradylglycerol is a glycerol compound where all three hydroxyl groups are substituted
    with acyl, alkyl, or alk-1-enyl groups at positions sn-1, sn-2, and sn-3.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a triradylglycerol, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS pattern for glycerol backbone
    # Glycerol backbone with carbons connected to oxygens
    glycerol_pattern = Chem.MolFromSmarts("C(O*)C(O*)C(O*)")
    matches = mol.GetSubstructMatches(glycerol_pattern)
    if not matches:
        return False, "No glycerol backbone with three substituted positions found"

    # For each match, check substituents
    for match in matches:
        c1_idx, c2_idx, c3_idx, o1_idx, o2_idx, o3_idx = match
        c_atoms = [mol.GetAtomWithIdx(idx) for idx in [c1_idx, c2_idx, c3_idx]]
        o_atoms = [mol.GetAtomWithIdx(idx) for idx in [o1_idx, o2_idx, o3_idx]]

        # Check that each oxygen is connected to a substituent (not hydrogen)
        for o_atom in o_atoms:
            connected_atoms = [nbr for nbr in o_atom.GetNeighbors() if nbr.GetAtomicNum() > 1 and nbr.GetIdx() not in match]
            if not connected_atoms:
                break  # Oxygen is not connected to any substituent
            substituent = connected_atoms[0]

            # Check if substituent is acyl, alkyl, or alk-1-enyl
            # Define SMARTS patterns for acyl, alkyl, alk-1-enyl groups
            acyl_pattern = Chem.MolFromSmarts("C(=O)[#6]")
            alkyl_pattern = Chem.MolFromSmarts("[#6][#6]")
            alk1enyl_pattern = Chem.MolFromSmarts("C=C[#6]")

            frag = Chem.PathToSubmol(mol, Chem.FindAtomEnvironmentOfRadiusN(mol, 1, substituent.GetIdx()))
            smiles_frag = Chem.MolToSmiles(frag)
            is_acyl = frag.HasSubstructMatch(acyl_pattern)
            is_alkyl = frag.HasSubstructMatch(alkyl_pattern)
            is_alk1enyl = frag.HasSubstructMatch(alk1enyl_pattern)
            if not (is_acyl or is_alkyl or is_alk1enyl):
                break  # Substituent is not acyl, alkyl, or alk-1-enyl
        else:
            # All oxygens are connected to valid substituents
            return True, "Contains glycerol backbone fully substituted with acyl, alkyl, or alk-1-enyl groups"

    return False, "Does not match triradylglycerol structure"

__metadata__ = {   'chemical_class': {   'id': '<CHEBI ID>',
                          'name': 'triradylglycerol',
                          'definition': 'A glycerol compound having one of three possible substituent '
                                        'groups - either acyl, alkyl, or alk-1-enyl - at each of the '
                                        'three possible positions sn-1, sn-2 or sn-3. Has functional parent '
                                        'glycerol (CHEBI:17754), children: triglyceride (CHEBI:17855). Parent: '
                                        'is_a glycerolipid (CHEBI:35741)'},
    }