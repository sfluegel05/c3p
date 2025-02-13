"""
Classifies: CHEBI:51689 enone
"""
"""
Classifies: CHEBI:24768 enone
Definition: An alpha,beta-unsaturated ketone of general formula R(1)R(2)C=CR(3)-C(=O)R(4) (R(4) ≠ H) 
in which the C=O function is conjugated to a C=C double bond at the alpha,beta position.
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_enone(smiles: str):
    """
    Determines if a molecule is an enone based on its SMILES string.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an enone, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Look for alpha,beta-unsaturated ketone substructure
    enone_pattern = Chem.MolFromSmarts("[CD3]=[CD2][C@]=O")
    enone_matches = mol.GetSubstructMatches(enone_pattern)
    
    if not enone_matches:
        return False, "No alpha,beta-unsaturated ketone substructure found"
    
    # Check for conjugation
    is_conjugated = False
    for match in enone_matches:
        c_alpha, c_beta, c_carbonyl = match
        
        # Trace bonds from alpha and beta carbons to check for conjugation
        alpha_conjugated = _is_conjugated(mol, c_alpha, c_beta, c_carbonyl)
        beta_conjugated = _is_conjugated(mol, c_beta, c_alpha, c_carbonyl)
        
        if alpha_conjugated and beta_conjugated:
            is_conjugated = True
            break
    
    if not is_conjugated:
        return False, "Alpha,beta-unsaturated ketone not conjugated"
    
    # Check for non-hydrogen substituent at carbonyl carbon
    carbonyl_atom = mol.GetAtomWithIdx(c_carbonyl)
    if all(nbr.GetAtomicNum() == 1 for nbr in carbonyl_atom.GetNeighbors()):
        return False, "No non-hydrogen substituent at carbonyl carbon"
    
    return True, "Molecule contains an alpha,beta-unsaturated ketone with conjugation"

def _is_conjugated(mol, start_atom_idx, next_atom_idx, end_atom_idx):
    """
    Helper function to check if a series of bonds is conjugated.

    Args:
        mol (rdkit.Chem.rdchem.Mol): RDKit molecule object
        start_atom_idx (int): Index of the starting atom
        next_atom_idx (int): Index of the next atom in the conjugated system
        end_atom_idx (int): Index of the ending atom (carbonyl carbon)

    Returns:
        bool: True if the series of bonds is conjugated, False otherwise
    """
    start_atom = mol.GetAtomWithIdx(start_atom_idx)
    next_atom = mol.GetAtomWithIdx(next_atom_idx)
    end_atom = mol.GetAtomWithIdx(end_atom_idx)
    
    # Check if start atom is sp2 hybridized and bonded to next atom
    if start_atom.GetHybridization() != Chem.HybridizationType.SP2:
        return False
    bond = mol.GetBondBetweenAtoms(start_atom_idx, next_atom_idx)
    if bond is None or bond.GetBondType() != Chem.BondType.DOUBLE:
        return False
    
    # Recursively check remaining atoms in conjugated system
    if next_atom_idx == end_atom_idx:
        return True
    else:
        for nbr_idx in next_atom.GetNeighbors():
            nbr = mol.GetAtomWithIdx(nbr_idx)
            if nbr_idx != start_atom_idx and nbr.GetHybridization() == Chem.HybridizationType.SP2:
                if _is_conjugated(mol, nbr_idx, end_atom_idx, end_atom_idx):
                    return True
    
    return False