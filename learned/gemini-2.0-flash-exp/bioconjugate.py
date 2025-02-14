"""
Classifies: CHEBI:64985 bioconjugate
"""
from rdkit import Chem
from rdkit.Chem import AllChem

def is_bioconjugate(smiles: str):
    """
    Determines if a molecule is a bioconjugate based on its SMILES string.
    A bioconjugate is defined as at least two biological molecules covalently linked.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a bioconjugate, False otherwise
        str: Reason for classification.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"

    # Define SMARTS patterns for common biological building blocks
    amino_acid_pattern = Chem.MolFromSmarts("[NX3;H2,H1,H0][CX4][CX3](=[OX1])[OX2;H1,H0]")  # Covers both N- and C- termini
    sugar_pattern = Chem.MolFromSmarts("OC[C,c]1([O,N,S][C,c]([O,N,S])[C,c]([O,N,S])[C,c]([O,N,S])[C,c]1[O,N,S])")  # Covers many sugar types
    nucleotide_base_pattern = Chem.MolFromSmarts("c1[nc][nc][nc]1")  # covers purine and pyrimidine bases
    phosphate_pattern = Chem.MolFromSmarts("P(=O)([O])([O])[O]")  # Covers phosphates
    fatty_acid_pattern = Chem.MolFromSmarts("C(=O)(O)C[C,c]([C,c])")
    
    patterns = [amino_acid_pattern, sugar_pattern, nucleotide_base_pattern, phosphate_pattern, fatty_acid_pattern]

    # Find matches for each pattern
    matches = []
    for pattern in patterns:
      if mol.HasSubstructMatch(pattern):
        matches.append(pattern)

    # Check if we have at least 2 different classes of substructures.
    if len(matches) < 2:
      return False, "Does not contain at least two different bio-molecule substructures"


    # Define SMARTS patterns for common bioconjugate linkages
    peptide_bond_pattern = Chem.MolFromSmarts("C(=O)N")
    glycosidic_bond_pattern = Chem.MolFromSmarts("C-O-[C,c]1([O,N,S][C,c]([O,N,S])[C,c]([O,N,S])[C,c]([O,N,S])[C,c]1[O,N,S])")  # C-O-C in a sugar context
    thioester_bond_pattern = Chem.MolFromSmarts("C(=O)S") # For CoA linkages
    phosphodiester_bond_pattern = Chem.MolFromSmarts("P(=O)(-O)-O") # For nucleotide linkages

    linkage_patterns = [peptide_bond_pattern, glycosidic_bond_pattern, thioester_bond_pattern, phosphodiester_bond_pattern]

    # Check if any linkage is present and connects different substructure types
    linked = False
    
    for linkage_pattern in linkage_patterns:
      linkage_matches = mol.GetSubstructMatches(linkage_pattern)
      
      if linkage_matches:
        for match in linkage_matches:
            # Get the atoms of the linkage.
            link_atoms = set(match)

            # Check if linkage connects two different substructures.
            connected_substructures = 0
            for pattern in matches:
                match_indices = mol.GetSubstructMatches(pattern)
                if match_indices:
                    for match_idx in match_indices:
                      substructure_atoms = set(match_idx)
                      if len(link_atoms.intersection(substructure_atoms)) > 0:
                        connected_substructures += 1
            if connected_substructures > 1:
                linked = True
                break
      if linked:
        break

    if not linked:
        return False, "Contains at least two different bio-molecule substructures, but not covalently linked."

    return True, "Contains at least two different covalently linked bio-molecule substructures."