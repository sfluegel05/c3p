"""
Classifies: CHEBI:64583 sphingomyelin
"""
"""
Classifies: sphingomyelin
"""
from rdkit import Chem

def is_sphingomyelin(smiles: str):
    """
    Determines if a molecule is a sphingomyelin based on its SMILES string.
    A sphingomyelin is defined as a phospholipid where the amino group of a sphingoid base 
    is in amide linkage with a fatty acid, and the terminal hydroxy group of 
    the sphingoid base is esterified to phosphorylcholine.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a sphingomyelin, False otherwise
        str: Reason for classification
    """

    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns
    # Sphingoid base: long-chain amino alcohol with variable chain length and unsaturation
    sphingoid_pattern = Chem.MolFromSmarts("N[C;@H][C;@H](O)[C;!H0]CCCCCCCCCCCC")
    # Amide linkage: amide bond between amino group and fatty acid
    amide_pattern = Chem.MolFromSmarts("N[C;!H0]=O")
    # Phosphorylcholine group
    phosphocholine_pattern = Chem.MolFromSmarts("O[P](=O)(OC[C,N](C)(C)C)O")

    # Ensure patterns are valid
    if sphingoid_pattern is None or amide_pattern is None or phosphocholine_pattern is None:
        return False, "Invalid SMARTS pattern(s)"

    # Search for sphingoid base
    sphingoid_matches = mol.GetSubstructMatches(sphingoid_pattern)
    if not sphingoid_matches:
        return False, "No sphingoid base found"

    # Search for amide linkage
    amide_matches = mol.GetSubstructMatches(amide_pattern)
    if not amide_matches:
        return False, "No amide linkage found"

    # Search for phosphorylcholine group
    phosphocholine_matches = mol.GetSubstructMatches(phosphocholine_pattern)
    if not phosphocholine_matches:
        return False, "No phosphorylcholine group found"

    # Verify connectivity
    is_valid_sphingomyelin = False
    for sphingoid_match in sphingoid_matches:
        sphingoid_n_idx = sphingoid_match[0]
        sphingoid_c1_idx = sphingoid_match[1]
        sphingoid_o_idx = sphingoid_match[2]
        sphingoid_c_chain_idx = sphingoid_match[3]

        # Check if amide nitrogen is the same as sphingoid base nitrogen
        for amide_match in amide_matches:
            amide_n_idx = amide_match[0]
            amide_c_idx = amide_match[1]

            if sphingoid_n_idx == amide_n_idx:
                # Check if phosphorylcholine is connected to sphingoid terminal hydroxyl
                for phosphocholine_match in phosphocholine_matches:
                    phospho_o_idx = phosphocholine_match[0]  # Oxygen attached to phosphate
                    # Check for bond between sphingoid O and phospho O
                    if mol.GetBondBetweenAtoms(sphingoid_o_idx, phospho_o_idx):
                        is_valid_sphingomyelin = True
                        break
            if is_valid_sphingomyelin:
                break
        if is_valid_sphingomyelin:
            break

    if is_valid_sphingomyelin:
        return True, "Molecule is a sphingomyelin with sphingoid base, amide-linked fatty acid, and phosphorylcholine group"
    else:
        return False, "Molecule does not meet the criteria for sphingomyelin"