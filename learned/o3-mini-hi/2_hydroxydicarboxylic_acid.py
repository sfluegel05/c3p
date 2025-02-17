"""
Classifies: CHEBI:50263 2-hydroxydicarboxylic acid
"""
"""
Classifies: 2-hydroxydicarboxylic acid
Definition: Any dicarboxylic acid carrying a hydroxy group on the carbon atom 
at position alpha to the carboxy group.
"""

from rdkit import Chem
from rdkit.Chem import rdchem

def is_2_hydroxydicarboxylic_acid(smiles: str):
    """
    Determines if a molecule is a 2-hydroxydicarboxylic acid based on its SMILES string.
    A 2-hydroxydicarboxylic acid is defined as any molecule with exactly two carboxy groups
    where at least one of the carboxyl carbons is connected to an sp3-hybridized carbon bearing an -OH substituent.
    
    Args:
        smiles (str): SMILES string of the molecule.
        
    Returns:
        bool: True if molecule is a 2-hydroxydicarboxylic acid, False otherwise.
        str : Reason for classification.
    """
    # Parse SMILES string and add explicit hydrogens for reliable substructure matching.
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    mol_h = Chem.AddHs(mol)
    
    # Define SMARTS for a carboxyl group: [CX3](=O)[OX2H]
    carboxyl_smarts = "[CX3](=O)[OX2H]"
    carboxyl_pattern = Chem.MolFromSmarts(carboxyl_smarts)
    carboxyl_matches = mol_h.GetSubstructMatches(carboxyl_pattern)
    
    # Expect exactly two carboxyl groups for a dicarboxylic acid.
    if len(carboxyl_matches) != 2:
        return False, f"Molecule has {len(carboxyl_matches)} carboxylic acid groups (exactly 2 required)"
    
    # Flag to indicate if any alpha carbon carries an -OH substituent
    alpha_hydroxy_found = False
    
    # Iterate over each carboxyl group match.
    # Each match is a tuple: (acid_carbon, carbonyl oxygen, carboxyl hydroxyl oxygen)
    for match in carboxyl_matches:
        acid_carbon_idx = match[0]
        acid_carbon = mol_h.GetAtomWithIdx(acid_carbon_idx)
        
        # For each neighbor of the acid carbon, consider it as a candidate alpha carbon.
        for neighbor in acid_carbon.GetNeighbors():
            # Only consider carbon atoms that are not part of the carboxyl group (exclude O's that belong to COOH)
            if neighbor.GetAtomicNum() != 6:
                continue
            # Ensure candidate alpha carbon is sp3 hybridized (to rule out double bonds etc.)
            if neighbor.GetHybridization() != rdchem.HybridizationType.SP3:
                continue

            # Search for an -OH group attached to the candidate alpha carbon.
            # We look for an oxygen neighbor that has at least one hydrogen attached.
            for nb in neighbor.GetNeighbors():
                if nb.GetAtomicNum() == 8:
                    # Exclude oxygen atoms that are directly in the carboxyl group by checking if nb is bonded to the acid carbon.
                    # For an -OH group on the candidate, the oxygen should be attached to the candidate carbon and carry at least one H.
                    if nb.GetTotalNumHs() >= 1:
                        alpha_hydroxy_found = True
                        break
            if alpha_hydroxy_found:
                break
        if alpha_hydroxy_found:
            break
            
    if not alpha_hydroxy_found:
        return False, "No alpha carbon adjacent to a carboxyl group carries a hydroxyl group"

    return True, "Molecule is a 2-hydroxydicarboxylic acid (dicarboxylic acid with an alpha hydroxy substituent)"

# Example test cases:
if __name__ == "__main__":
    test_smiles = [
        "C[C@H](C(O)=O)[C@@](C)(O)C(O)=O",  # (2R,3S)-2,3-dimethylmalic acid
        "O[C@@H](CCC(O)=O)C(O)=O",            # (S)-2-hydroxyglutaric acid
        "O[C@H](CC(O)=O)C(O)=O",              # (R)-malic acid
    ]
    
    for smi in test_smiles:
        is_class, reason = is_2_hydroxydicarboxylic_acid(smi)
        print(f"SMILES: {smi} --> {is_class}, Reason: {reason}")