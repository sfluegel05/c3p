"""
Classifies: CHEBI:46874 alpha-amino acid ester
"""
"""
Classifies: CHEBI:85259 alpha-amino acid ester
"""
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors

def is_alpha_amino_acid_ester(smiles: str):
    """
    Determines if a molecule is an alpha-amino acid ester based on its SMILES string.
    An alpha-amino acid ester is formed by the condensation of an alpha-amino acid with an alcohol.

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

    # More flexible ester pattern to catch different configurations
    ester_pattern = Chem.MolFromSmarts("[CX3](=[OX1])[OX2][#6]")
    ester_matches = mol.GetSubstructMatches(ester_pattern)
    if len(ester_matches) == 0:
        return False, "No ester group found"

    # More flexible alpha-amino acid pattern
    # Accounts for substituted amines, different carboxylate forms, and stereochemistry
    amino_acid_pattern = Chem.MolFromSmarts("[NX3;H2,H1,H0][CX4;H1,H2][CX3](=[OX1])[OX2H0]")
    amino_acid_matches = mol.GetSubstructMatches(amino_acid_pattern)
    
    # Alternative pattern for cases where the carboxylate is part of a larger structure
    alt_amino_acid_pattern = Chem.MolFromSmarts("[NX3;H2,H1,H0][CX4;H1,H2][CX3](=[OX1])[OX2][#6]")
    alt_amino_acid_matches = mol.GetSubstructMatches(alt_amino_acid_pattern)
    
    # Combine matches
    all_aa_matches = amino_acid_matches + alt_amino_acid_matches
    if len(all_aa_matches) == 0:
        return False, "No alpha-amino acid moiety found"

    # Check if any ester group is connected to the carboxylate of any amino acid
    ester_connected = False
    for aa_match in all_aa_matches:
        # Get the carboxylate carbon and oxygen indices
        carboxyl_c = aa_match[2]
        carboxyl_o = aa_match[3]
        
        # Check if any ester group is connected to the carboxylate carbon or oxygen
        for ester_match in ester_matches:
            ester_c = ester_match[0]
            ester_o = ester_match[1]
            if ester_c == carboxyl_c or ester_o == carboxyl_o:
                ester_connected = True
                break
        if ester_connected:
            break

    if not ester_connected:
        return False, "Ester group not connected to the carboxylate of the amino acid"

    # Check if the amino acid has an alpha-carbon with at least one hydrogen
    # and that it's connected to both the nitrogen and carboxylate
    alpha_carbon_valid = False
    for aa_match in all_aa_matches:
        alpha_c = aa_match[1]
        atom = mol.GetAtomWithIdx(alpha_c)
        if atom.GetTotalNumHs() > 0:
            # Verify connections
            neighbors = [n.GetIdx() for n in atom.GetNeighbors()]
            if aa_match[0] in neighbors and aa_match[2] in neighbors:
                alpha_carbon_valid = True
                break

    if not alpha_carbon_valid:
        return False, "No valid alpha-carbon found in the amino acid moiety"

    return True, "Contains alpha-amino acid moiety with ester group connected to the carboxylate"