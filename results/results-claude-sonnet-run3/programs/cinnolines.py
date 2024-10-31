from rdkit import Chem
from rdkit.Chem import AllChem

def is_cinnolines(smiles: str):
    """
    Determines if a molecule is a cinnoline (1,2-diaza analogue of naphthalene) or derivative.
    
    Args:
        smiles (str): SMILES string of the molecule
        
    Returns:
        bool: True if molecule is a cinnoline, False otherwise
        str: Reason for classification
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return None, "Invalid SMILES string"
        
    # Match cinnoline core structure using SMARTS pattern
    # Matches fused bicyclic system with N=N in specific position
    # The pattern describes:
    # - A fused bicyclic system where one ring is benzene
    # - Two adjacent nitrogen atoms (N=N) in positions 1,2 of the second ring
    cinnoline_pattern = Chem.MolFromSmarts('[c]1[c][c][c][c]2[c]1[c][n][n][c]2')
    
    if mol.HasSubstructMatch(cinnoline_pattern):
        matches = mol.GetSubstructMatches(cinnoline_pattern)
        for match in matches:
            # Get the atoms in the matched pattern
            matched_atoms = [mol.GetAtomWithIdx(i) for i in match]
            
            # Verify we have a benzene ring fused to a heterocycle with N=N
            benzene_atoms = matched_atoms[:6]
            if all(atom.GetIsAromatic() for atom in benzene_atoms):
                # Check for the N=N pattern
                n_atoms = [atom for atom in matched_atoms[6:] if atom.GetSymbol() == 'N']
                if len(n_atoms) == 2:
                    # Verify the nitrogens are adjacent
                    n1_idx = mol.GetAtomWithIdx(match[7])
                    n2_idx = mol.GetAtomWithIdx(match[8])
                    bond = mol.GetBondBetweenAtoms(match[7], match[8])
                    if bond is not None and bond.GetBondTypeAsDouble() >= 1.5:
                        return True, "Contains cinnoline core structure"
                            
    return False, "Does not contain cinnoline core structure"
# Pr=None
# Recall=0.0