"""
Classifies: CHEBI:87657 octanoate ester
"""
"""
Classifies: CHEBI:38145 octanoate ester
Any fatty acid ester in which the carboxylic acid component is octanoic acid (caprylic acid).
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_octanoate_ester(smiles: str):
    """
    Determines if a molecule is an octanoate ester based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an octanoate ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Find all carboxylic acid groups
    acid_pattern = Chem.MolFromSmarts("C(=O)O")
    acid_matches = mol.GetSubstructMatches(acid_pattern)
    
    # Check if any carboxylic acid group has an 8-carbon chain attached
    for match in acid_matches:
        acid_atom = mol.GetAtomWithIdx(list(match)[0])
        neighbor_atoms = [mol.GetAtomWithIdx(nbr_idx) for nbr_idx in acid_atom.GetNeighbors()]
        
        # Check if any neighbor is part of an 8-carbon chain
        for nbr_atom in neighbor_atoms:
            if nbr_atom.GetAtomicNum() == 6:  # Carbon
                chain = Chem.FindAllPathsOfLengthN(mol, nbr_atom.GetIdx(), 8, useBonds=True)
                if chain:
                    # Check if the carboxylic acid is part of an ester linkage
                    for path in chain:
                        if any(mol.GetBondBetweenAtoms(path[i], path[i+1]).GetBondType() == Chem.BondType.DOUBLE for i in range(len(path)-1)):
                            continue  # Double bond present, not an ester
                        elif any(mol.GetAtomWithIdx(idx).GetAtomicNum() != 6 and mol.GetAtomWithIdx(idx).GetAtomicNum() != 8 for idx in path):
                            continue  # Non-carbon/oxygen atom in chain
                        else:
                            # Found an octanoate ester
                            return True, "Molecule contains an octanoic acid component attached via an ester linkage"
    
    return False, "No octanoate ester group found"