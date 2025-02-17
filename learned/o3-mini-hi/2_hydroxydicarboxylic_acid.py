"""
Classifies: CHEBI:50263 2-hydroxydicarboxylic acid
"""
"""
Classifies: 2-hydroxydicarboxylic acid
Definition: Any dicarboxylic acid carrying a hydroxy group on the carbon atom 
at position alpha to a carboxy group.
"""

from rdkit import Chem

def is_2_hydroxydicarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a 2-hydroxydicarboxylic acid based on its SMILES string.
    A 2-hydroxydicarboxylic acid must have exactly two carboxyl groups (–COOH) and at least one
    of the carboxyl groups is adjacent (alpha carbon) to a carbon that carries a hydroxyl (-OH) group.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is a 2-hydroxydicarboxylic acid, False otherwise.
        str : Reason for classification.
    """
    # Parse the SMILES string.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Add explicit hydrogens for reliable substructure matching.
    mol_h = Chem.AddHs(mol)
    
    # Define a SMARTS pattern for a carboxylic acid group.
    # We use a loose pattern "C(=O)O" and then check that the 'O' has at least one H.
    carboxyl_smarts = "C(=O)O"
    carboxyl_pattern = Chem.MolFromSmarts(carboxyl_smarts)
    carboxyl_matches = mol_h.GetSubstructMatches(carboxyl_pattern)
    
    # Filter matches to ensure the hydroxyl oxygen actually carries an explicit hydrogen.
    valid_carboxyl_indices = []
    for match in carboxyl_matches:
        # match: (carboxyl carbon, carbonyl oxygen, hydroxyl oxygen)
        acid_carbon_idx, carbonyl_O_idx, hydroxyl_O_idx = match
        hydroxyl_O = mol_h.GetAtomWithIdx(hydroxyl_O_idx)
        # Check if this oxygen has any hydrogen neighbor.
        has_H = any(neighbor.GetAtomicNum() == 1 for neighbor in hydroxyl_O.GetNeighbors())
        if has_H:
            valid_carboxyl_indices.append(acid_carbon_idx)
    
    # Remove duplicate carboxyl groups if matching overlaps (using set)
    valid_carboxyl_indices = list(set(valid_carboxyl_indices))
    
    # For a dicarboxylic acid, we need exactly two distinct carboxyl groups.
    if len(valid_carboxyl_indices) != 2:
        return False, f"Molecule has {len(valid_carboxyl_indices)} valid carboxylic acid groups (exactly 2 required)"
    
    # Function to check if a given carbon atom (candidate alpha carbon) has an OH substituent.
    def has_OH(candidate):
        for nbr in candidate.GetNeighbors():
            if nbr.GetAtomicNum() == 8:
                # Check that the bond is single and that oxygen has at least one hydrogen.
                bond = mol_h.GetBondBetweenAtoms(candidate.GetIdx(), nbr.GetIdx())
                if bond is not None and bond.GetBondType()==Chem.BondType.SINGLE:
                    if any(x.GetAtomicNum() == 1 for x in nbr.GetNeighbors()):
                        return True
        return False
    
    # Flag to indicate if at least one carboxyl group meets the condition.
    alpha_OH_found = False
    
    # For each valid carboxyl group, consider carbons adjacent to it that are NOT the carbonyl or hydroxyl oxygens.
    for acid_carbon_idx in valid_carboxyl_indices:
        acid_carbon = mol_h.GetAtomWithIdx(acid_carbon_idx)
        for neighbor in acid_carbon.GetNeighbors():
            # We are only interested in carbon neighbors.
            if neighbor.GetAtomicNum() != 6:
                continue
            # Candidate alpha carbon found; check if it has an –OH group.
            if has_OH(neighbor):
                alpha_OH_found = True
                break
        if alpha_OH_found:
            break

    if not alpha_OH_found:
        return False, "No alpha carbon adjacent to a carboxyl group carries a hydroxyl group"
    
    return True, "Molecule is a 2-hydroxydicarboxylic acid with an alpha hydroxy substituent"

# Example test cases from the user-provided examples:
if __name__ == "__main__":
    test_smiles = [
        "C[C@H](C(O)=O)[C@@](C)(O)C(O)=O",            # (2R,3S)-2,3-dimethylmalic acid
        "O[C@@H](CCC(O)=O)C(O)=O",                      # (S)-2-hydroxyglutaric acid
        "O[C@H](CC(O)=O)C(O)=O",                        # (R)-malic acid
        "CCC(C(O)C(O)=O)C(O)=O",                        # 3-ethylmalic acid; should fail (OH is not alpha)
        "CC(C)[C@@H]([C@@H](O)C(O)=O)C(O)=O",           # (2R,3S)-3-isopropylmalic acid
        "OC(=O)C(\\O)=C(/O)C(O)=O"                      # dihydroxyfumaric acid; alpha not sp3 so should fail
    ]
    
    for smi in test_smiles:
        is_class, reason = is_2_hydroxydicarboxylic_acid(smi)
        print(f"SMILES: {smi} --> {is_class}, Reason: {reason}")