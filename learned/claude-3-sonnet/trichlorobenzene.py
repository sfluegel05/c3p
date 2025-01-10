"""
Classifies: CHEBI:27096 trichlorobenzene
"""
"""
Classifies: trichlorobenzene
Definition: Any member of the class of chlorobenzenes carrying three chloro substituents at unspecified positions.
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_trichlorobenzene(smiles: str):
    """
    Determines if a molecule contains a trichlorobenzene substructure
    (benzene ring with exactly 3 chlorine substituents at any positions).
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule contains a trichlorobenzene substructure, False otherwise
        str: Reason for classification
    """
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Define SMARTS patterns for different trichlorobenzene substitution patterns
    # Note: [cH0] means a non-hydrogen-substituted aromatic carbon
    patterns = [
        # Generic pattern for benzene with exactly 3 chlorines
        "[cH0]1([Cl])[cH1,cH0][cH1,cH0][cH0]([Cl])[cH1,cH0][cH0]1[Cl]",
        # 1,2,3-trichlorobenzene pattern
        "c1c(Cl)c(Cl)c(Cl)cc1",
        # 1,2,4-trichlorobenzene pattern
        "c1c(Cl)cc(Cl)c(Cl)c1",
        # 1,3,5-trichlorobenzene pattern
        "c1c(Cl)cc(Cl)cc1Cl"
    ]
    
    # Check each pattern
    for pattern in patterns:
        pattern_mol = Chem.MolFromSmarts(pattern)
        if pattern_mol is not None and mol.HasSubstructMatch(pattern_mol):
            # Count chlorines in each matching substructure
            matches = mol.GetSubstructMatches(pattern_mol)
            for match in matches:
                # Get atoms in the match
                match_atoms = set(match)
                # Count chlorines connected to the matched benzene ring
                chlorine_count = 0
                for atom_idx in match_atoms:
                    atom = mol.GetAtomWithIdx(atom_idx)
                    if atom.GetSymbol() == "C":  # Only check carbons
                        for neighbor in atom.GetNeighbors():
                            if neighbor.GetSymbol() == "Cl" and neighbor.GetIdx() not in match_atoms:
                                chlorine_count += 1
                
                if chlorine_count == 3:
                    positions = []
                    # Determine substitution pattern
                    if mol.HasSubstructMatch(Chem.MolFromSmarts("c1c(Cl)c(Cl)c(Cl)cc1")):
                        positions = "1,2,3"
                    elif mol.HasSubstructMatch(Chem.MolFromSmarts("c1c(Cl)cc(Cl)c(Cl)c1")):
                        positions = "1,2,4"
                    elif mol.HasSubstructMatch(Chem.MolFromSmarts("c1c(Cl)cc(Cl)cc1Cl")):
                        positions = "1,3,5"
                    else:
                        positions = "unspecified"
                        
                    return True, f"Contains {positions}-trichlorobenzene substructure"
    
    return False, "No trichlorobenzene substructure found"