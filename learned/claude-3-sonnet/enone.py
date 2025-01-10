"""
Classifies: CHEBI:51689 enone
"""
"""
Classifies: enone (alpha,beta-unsaturated ketone)
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_enone(smiles: str):
    """
    Determines if a molecule contains an enone group based on its SMILES string.
    An enone is an alpha,beta-unsaturated ketone with general formula R(1)R(2)C=CR(3)-C(=O)R(4)
    where R(4) ≠ H, and the C=O is conjugated to a C=C double bond.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule contains an enone group, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for enone groups
    # Pattern 1: Linear enone pattern
    # [C;!$(C=O)] = carbon not part of a carbonyl
    # [CD1,CD2,CD3,CD4] = carbon with 1-4 connections (covers both chain and ring)
    patterns = [
        # Basic enone pattern
        "[C;!$(C=O)]=[C;!$(C=O)]-C(=[OD1])-[#6]",
        # Cyclic enone pattern
        "[C;R]=[C;R]-C(=[OD1])-[#6]",
        # Alternative pattern for complex rings
        "[C;!$(C=O)]=[C]-C(=O)-[#6;!$([#6][OD1])]"
    ]

    for pattern in patterns:
        substructure = Chem.MolFromSmarts(pattern)
        if substructure is None:
            continue
            
        matches = mol.GetSubstructMatches(substructure)
        
        if matches:
            for match in matches:
                # Verify conjugation and ketone nature
                try:
                    # Get the relevant atoms
                    c_beta = mol.GetAtomWithIdx(match[0])  # Beta carbon
                    c_alpha = mol.GetAtomWithIdx(match[1]) # Alpha carbon
                    c_carbonyl = mol.GetAtomWithIdx(match[2]) # Carbonyl carbon
                    c_substituent = mol.GetAtomWithIdx(match[3]) # R group on carbonyl
                    
                    # Check bond types
                    double_bond = mol.GetBondBetweenAtoms(match[0], match[1])
                    single_bond = mol.GetBondBetweenAtoms(match[1], match[2])
                    
                    # Verify conjugation pattern
                    if (double_bond.GetBondType() == Chem.BondType.DOUBLE and 
                        single_bond.GetBondType() == Chem.BondType.SINGLE):
                        
                        # Verify it's a ketone (not an aldehyde)
                        if c_carbonyl.GetTotalNumHs() == 0:
                            # Additional check for carbonyl oxygen
                            for neighbor in c_carbonyl.GetNeighbors():
                                if neighbor.GetAtomicNum() == 8 and neighbor.GetTotalNumHs() == 0:
                                    return True, "Contains alpha,beta-unsaturated ketone (enone) group"
                
                except (IndexError, AttributeError):
                    continue

    return False, "No enone group found or incorrect conjugation pattern"