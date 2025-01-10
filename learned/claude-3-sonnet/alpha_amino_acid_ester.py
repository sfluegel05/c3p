"""
Classifies: CHEBI:46874 alpha-amino acid ester
"""
"""
Classifies: alpha-amino acid ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_alpha_amino_acid_ester(smiles: str):
    """
    Determines if a molecule is an alpha-amino acid ester based on its SMILES string.
    An alpha-amino acid ester has an amine and ester group attached to the same carbon.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is an alpha-amino acid ester, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Look for ester pattern (-C(=O)OR)
    ester_pattern = Chem.MolFromSmarts("[C;!$(C(=O)N)](=O)O[C,H]")
    if not mol.HasSubstructMatch(ester_pattern):
        return False, "No ester group found"

    # Look for basic alpha-amino acid ester pattern
    # Carbon with both amine and ester attached
    basic_pattern = Chem.MolFromSmarts("[NX3;!$(NC=O)][CH1,CH2][C](=O)O[C,H]")
    
    # Alternative pattern for N-substituted amino acid esters
    substituted_pattern = Chem.MolFromSmarts("[#7;!$(NC=O);!$(N=*);!$(N#*)][CH1,CH2][C](=O)O[C,H]")
    
    # Check if either pattern matches
    if not (mol.HasSubstructMatch(basic_pattern) or mol.HasSubstructMatch(substituted_pattern)):
        return False, "No alpha-amino acid ester pattern found"
    
    # Get all atoms that match either pattern
    matches = []
    for pattern in [basic_pattern, substituted_pattern]:
        if mol.HasSubstructMatch(pattern):
            matches.extend(mol.GetSubstructMatches(pattern))
    
    # Verify the matches
    for match in matches:
        n_atom = mol.GetAtomWithIdx(match[0])
        c_atom = mol.GetAtomWithIdx(match[1])
        
        # Check nitrogen hybridization (should be sp3)
        if n_atom.GetHybridization() != Chem.HybridizationType.SP3:
            continue
            
        # Check if nitrogen is part of an amide (exclude)
        amide_pattern = Chem.MolFromSmarts("[NX3][C](=O)")
        is_amide = False
        for neighbor in n_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 6:  # Carbon
                for bond in neighbor.GetBonds():
                    other_atom = bond.GetOtherAtom(neighbor)
                    if other_atom.GetAtomicNum() == 8 and bond.GetBondType() == Chem.BondType.DOUBLE:
                        is_amide = True
                        break
        if is_amide:
            continue
            
        # Verify carbon is really alpha to both N and C(=O)
        has_n = False
        has_ester = False
        for neighbor in c_atom.GetNeighbors():
            if neighbor.GetAtomicNum() == 7:  # Nitrogen
                has_n = True
            elif neighbor.GetAtomicNum() == 6:  # Carbon
                for bond in neighbor.GetBonds():
                    other_atom = bond.GetOtherAtom(neighbor)
                    if other_atom.GetAtomicNum() == 8 and bond.GetBondType() == Chem.BondType.DOUBLE:
                        has_ester = True
        
        if has_n and has_ester:
            return True, "Contains alpha-amino acid ester pattern"
            
    return False, "No valid alpha-amino acid ester pattern found"