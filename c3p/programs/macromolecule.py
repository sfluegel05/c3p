"""
Classifies: CHEBI:33839 macromolecule
"""
from rdkit import Chem
from rdkit.Chem import rdMolDescriptors

def is_macromolecule(smiles: str):
    """
    Determines if a molecule is a macromolecule based on its SMILES string.
    Considers molecular weight, presence of typical recurring monomer linkages, 
    and other structural motifs to determine macromolecule classification.

    Args:
        smiles (str): SMILES string of the molecule

    Returns:
        bool: True if molecule is a macromolecule, False otherwise
        str: Reason for classification
    """
    
    # Parse SMILES
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return False, "Invalid SMILES string"
    
    # Calculate molecular weight
    mol_wt = rdMolDescriptors.CalcExactMolWt(mol)
    min_macromolecule_weight = 500  # Consider molecular weight threshold
    if mol_wt < min_macromolecule_weight:
        return False, f"Molecular weight {mol_wt:.2f} Da is less than the minimum macromolecule threshold {min_macromolecule_weight} Da"

    # Detect repeating linkages typical in macromolecules
    try:
        # Improved patterns for identifying macromolecular linkages
        peptide_bond_pattern = Chem.MolFromSmarts('C(=O)N[C@H]')
        glycosidic_bond_pattern = Chem.MolFromSmarts('O[C@@H]1C(O)C([C@@H]([C@H]1)O)O')
        nucleic_acid_pattern = Chem.MolFromSmarts('C1=NC=CC(=O)N1')
        
        # Expanded search for recurring patterns
        peptide_matches = mol.GetSubstructMatches(peptide_bond_pattern)
        glycosidic_matches = mol.GetSubstructMatches(glycosidic_bond_pattern)
        nucleic_matches = mol.GetSubstructMatches(nucleic_acid_pattern)

        # More vigorous pattern detection for typical macromolecules
        if len(peptide_matches) > 5 or len(glycosidic_matches) > 3 or len(nucleic_matches) > 2:
            if len(peptide_matches) + len(glycosidic_matches) + len(nucleic_matches) > 10:
                return True, "Considered a macromolecule based on recurrent macromolecular linkages."
        
    except Exception as e:
        return False, f"Failed to match structures due to error: {str(e)}"
    
    # Check for additional macromolecular characteristics
    if mol.GetNumAtoms() > 100:
        # Heuristic: Large structures often correlate with macromolecules
        return True, "Large molecular structure identified as potential macromolecule."

    return False, "No sufficient recurring motifs typical of macromolecules detected."