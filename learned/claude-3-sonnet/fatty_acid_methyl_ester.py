"""
Classifies: CHEBI:4986 fatty acid methyl ester
"""
"""
Classifies: CHEBI:35368 Fatty acid methyl ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_fatty_acid_methyl_ester(smiles: str):
    """
    Determines if a molecule is a fatty acid methyl ester based on its SMILES string.
    A fatty acid methyl ester is a fatty acid ester obtained by the formal condensation 
    of a fatty acid with methanol.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a fatty acid methyl ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for ester group (-C(=O)O-) with a methyl group (-CH3) attached
    ester_pattern = Chem.MolFromSmarts("[C:1](=[O:2])O[CH3:3]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)

    # No ester group or methyl group found
    if not ester_matches:
        return False, "No ester group with methyl group found"

    # Check if the ester and methyl groups are connected
    for match in ester_matches:
        if mol.GetBondBetweenAtoms(match[2], match[0]).GetIdx() >= 0:
            break  # Found a valid ester-methyl pair
    else:
        return False, "Ester group and methyl group are not connected"

    # Look for fatty acid chain (long carbon chain with optional double bonds)
    fatty_acid_pattern = Chem.MolFromSmarts("[C;H3]([C;H2])([C;H2])[C;H2]~[C;H2]~[C;H2]~[C;H2]")
    fatty_acid_matches = mol.GetSubstructMatches(fatty_acid_pattern)

    # No fatty acid chain found
    if not fatty_acid_matches:
        return False, "No fatty acid chain found"

    # Check if the fatty acid chain and ester group are connected
    for match in fatty_acid_matches:
        for atom_idx in match:
            if mol.GetBondBetweenAtoms(atom_idx, match[1]).GetIdx() >= 0:
                break  # Found a connection between fatty acid chain and ester
        else:
            continue
        break
    else:
        return False, "Fatty acid chain and ester group are not connected"

    return True, "Contains a fatty acid chain connected to a methyl ester group"