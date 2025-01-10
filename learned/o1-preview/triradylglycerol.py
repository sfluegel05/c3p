"""
Classifies: CHEBI:76579 triradylglycerol
"""
"""
Classifies: CHEBI:<ID> triradylglycerol
"""
from rdkit import Chem

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

    # Define SMARTS pattern for glycerol backbone with substituted oxygens
    glycerol_pattern = Chem.MolFromSmarts("[CH2](O[*])[CH](O[*])[CH2](O[*])")
    
    matches = mol.GetSubstructMatches(glycerol_pattern)
    if not matches:
        return False, "No glycerol backbone with substituted oxygens found"

    # For each match, check substituents
    for match in matches:
        # The pattern has 6 atoms: C1-O1-C2-O2-C3-O3
        c1_idx, o1_idx, c2_idx, o2_idx, c3_idx, o3_idx = match
        c_indices = [c1_idx, c2_idx, c3_idx]
        o_indices = [o1_idx, o2_idx, o3_idx]

        # For each oxygen, check that it is connected to a substituent
        for o_idx in o_indices:
            o_atom = mol.GetAtomWithIdx(o_idx)
            substituent_found = False
            for nbr in o_atom.GetNeighbors():
                nbr_idx = nbr.GetIdx()
                if nbr_idx not in c_indices:
                    # Found substituent connected to oxygen
                    substituent_found = True
                    
                    # Get the fragment connected via oxygen excluding the glycerol backbone
                    fragment = Chem.PathToSubmol(mol, Chem.FindAtomEnvironmentOfRadiusN(mol, 5, nbr_idx))
                    
                    # Define SMARTS patterns
                    acyl_pattern = Chem.MolFromSmarts("C(=O)[#6]")
                    alkyl_pattern = Chem.MolFromSmarts("[#6][#6]")
                    alk1enyl_pattern = Chem.MolFromSmarts("C=C[#6]")
                    
                    is_acyl = fragment.HasSubstructMatch(acyl_pattern)
                    is_alkyl = fragment.HasSubstructMatch(alkyl_pattern)
                    is_alk1enyl = fragment.HasSubstructMatch(alk1enyl_pattern)

                    if not (is_acyl or is_alkyl or is_alk1enyl):
                        return False, "Substituent attached to oxygen is not acyl, alkyl, or alk-1-enyl group"
                    break
            if not substituent_found:
                return False, "Oxygen is not connected to any substituent"
        # All oxygens have valid substituents
        return True, "Contains glycerol backbone fully substituted with acyl, alkyl, or alk-1-enyl groups"

    # If no matches passed the checks
    return False, "Does not match triradylglycerol structure"

__metadata__ = {   'chemical_class': {   'id': '<CHEBI ID>',
                          'name': 'triradylglycerol',
                          'definition': 'A glycerol compound having one of three possible substituent '
                                        'groups - either acyl, alkyl, or alk-1-enyl - at each of the '
                                        'three possible positions sn-1, sn-2 or sn-3. Has functional parent '
                                        'glycerol (CHEBI:17754), children: triglyceride (CHEBI:17855). Parent: '
                                        'is_a glycerolipid (CHEBI:35741)'},
    }